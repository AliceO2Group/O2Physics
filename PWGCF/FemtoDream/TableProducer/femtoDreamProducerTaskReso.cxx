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

#include "Common/Core/Zorro.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

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
  Configurable<bool> confIsActivateV0{"confIsActivateV0", true, "Activate filling of V0 (Lambdas) into femtodream tables"};
  Configurable<bool> confIsActivateV0K0S{"confIsActivateV0K0S", true, "Activate filling of V0 into femtodream tables"};
  Configurable<bool> confIsActivateKStar{"confIsActivateKStar", true, "Activate filling of KStars into femtodream tables"};
  Configurable<bool> confIsActivatePhi{"confIsActivatePhi", true, "Activates cuts on Phi's and fills tables"};
  Configurable<bool> confIsActivateXi{"confIsActivateXi", false, "Activate filling of Xis into femtodream tables"};
  Configurable<bool> confIsActivateOmega{"confIsActivateOmega", false, "Activate filling of Omegas into femtodream tables"};
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

  FemtoDreamV0Selection LambdaCuts;
  FemtoDreamV0Selection K0SCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confLambdaSign{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0Sign, "confLambda"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaPtMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMin, "confLambda"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaPtMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMax, "confLambda"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaEtaMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0etaMax, "confLambda"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaDCADaughMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DCADaughMax, "confLambda"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaCPAMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0CPAMin, "confLambda"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaTranRadMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMin, "confLambda"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaTranRadMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMax, "confLambda"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaDecVtxMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DecVtxMax, "confLambda"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<float> confLambdaInvMassLowLimit{"confLambdaInvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> confLambdaInvMassUpLimit{"confLambdaInvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};
    Configurable<bool> confLambdaRejectKaons{"confLambdaRejectKaons", false, "Switch to reject kaons"};
    Configurable<bool> confLambdaRejectLambdas{"confLambdaRejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
    Configurable<float> confLambdaInvKaonMassLowLimit{"confLambdaInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> confLambdaInvKaonMassUpLimit{"confLambdaInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    Configurable<std::vector<float>> confLambdaChildSign{"confLambdaChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> confLambdaChildEtaMax{"confLambdaChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> confLambdaChildTPCnClsMin{"confLambdaChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confLambdaChildDCAMin{"confLambdaChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confLambdaChildPIDnSigmaMax{"confLambdaChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confLambdaChildPIDspecies{"confLambdaChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

    // cuts and object for v0 2
    Configurable<std::vector<float>> confK0shortSign{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0Sign, "confK0short"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortPtMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMin, "confK0short"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortPtMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMax, "confK0short"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortEtaMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0etaMax, "confK0short"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortDCADaughMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DCADaughMax, "confK0short"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortCPAMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0CPAMin, "confK0short"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortTranRadMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMin, "confK0short"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortTranRadMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMax, "confK0short"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortDecVtxMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DecVtxMax, "confK0short"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<float> confK0shortInvMassLowLimit{"confK0shortInvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> confK0shortInvMassUpLimit{"confK0shortInvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};
    Configurable<bool> confK0shortRejectKaons{"confK0shortRejectKaons", false, "Switch to reject kaons"};
    Configurable<bool> confK0shortRejectLambdas{"confK0shortRejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
    Configurable<float> confK0shortInvKaonMassLowLimit{"confK0shortInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> confK0shortInvKaonMassUpLimit{"confK0shortInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    Configurable<std::vector<float>> confK0shortChildSign{"confK0shortChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> confK0shortChildEtaMax{"confK0shortChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> confK0shortChildTPCnClsMin{"confK0shortChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confK0shortChildDCAMin{"confK0shortChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confK0shortChildPIDnSigmaMax{"confK0shortChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confK0shortChildPIDspecies{"confK0shortChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};
  } V0Sel;

  FemtoDreamCascadeSelection xiCuts;
  FemtoDreamCascadeSelection omegaCuts;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = std::string("Cascade");
    // Xi Selection
    Configurable<float> confXiInvMassLowLimit{"confXiInvMassLowLimit", 1.2, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confXiInvMassUpLimit{"confXiInvMassUpLimit", 1.5, "Upper limit of the Cascade invariant mass"};
    // Cascade
    Configurable<std::vector<float>> confXiSign{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeSign, "confXi"), std::vector<float>{-1, 1}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeSign, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiPtMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMin, "confXi"), std::vector<float>{0.3f, 0.4f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiPtMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMax, "confXi"), std::vector<float>{5.5f, 6.0f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiEtaMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeEtaMax, "confXi"), std::vector<float>{0.8f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeEtaMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiDCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDCADaughMax, "confXi"), std::vector<float>{1.f, 1.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiCPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeCPAMin, "confXi"), std::vector<float>{0.99f, 0.95f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeCPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiTranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMin, "confXi"), std::vector<float>{0.2f, 0.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiTranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMax, "confXi"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confXiDecVtxMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDecVtxMax, "confXi"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDecVtxMax, "Cascade selection: ")};

    // Cascade v0 daughters
    Configurable<std::vector<float>> confXiV0DCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "confXi"), std::vector<float>{1.2f, 1.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confXiV0CPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0CPAMin, "confXi"), std::vector<float>{0.99f, 0.995f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0CPAMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confXiV0TranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "confXi"), std::vector<float>{0.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confXiV0TranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "confXi"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confXiV0DCAtoPVMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "confXi"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confXiV0DCAtoPVMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "confXi"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "CascV0 selection: ")};
    Configurable<float> confXiV0InvMassLowLimit{"confXiV0InvMassLowLimit", 1.011461, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confXiV0InvMassUpLimit{"confXiV0InvMassUpLimit", 1.027461, "Upper limit of the Cascade invariant mass"};
    // Cascade Daughter Tracks
    Configurable<std::vector<float>> confXiV0ChildSign{"confXiV0ChildSign", std::vector<float>{-1, 1}, "CascV0 Child sel: Charge"};
    Configurable<std::vector<float>> confXiV0ChildPtMin{"confXiV0ChildPtMin", std::vector<float>{0.8f}, "CascV0 Child sel: min pt"};
    Configurable<std::vector<float>> confXiV0ChildEtaMax{"confXiV0ChildEtaMax", std::vector<float>{0.8f}, "CascV0 Child sel: max eta"};
    Configurable<std::vector<float>> confXiV0ChildTPCnClsMin{"confXiV0ChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "CascV0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confXiV0ChildDCAMin{"confXiV0ChildDCAMin", std::vector<float>{0.05f, 0.06f}, "CascV0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confXiV0ChildPIDnSigmaMax{"confXiV0ChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "CascV0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confXiV0ChildPIDspecies{"confXiV0ChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "CascV0 Child sel: Particles species for PID"};
    // Cascade Bachelor Track
    Configurable<std::vector<float>> confXiBachelorSign{"confXiBachelorSign", std::vector<float>{-1, 1}, "Cascade Bachelor sel: Charge"};
    Configurable<std::vector<float>> confXiBachelorPtMin{"confXiBachelorPtMin", std::vector<float>{0.8f}, "Cascade Bachelor sel: min pt"};
    Configurable<std::vector<float>> confXiBachelorEtaMax{"confXiBachelorEtaMax", std::vector<float>{0.8f}, "Cascade Bachelor sel: max eta"};
    Configurable<std::vector<float>> confXiBachelorTPCnClsMin{"confXiBachelorTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Cascade Bachelor sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confXiBachelorDCAMin{"confXiBachelorDCAMin", std::vector<float>{0.05f, 0.06f}, "Cascade Bachelor sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confXiBachelorPIDnSigmaMax{"confXiBachelorPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "Cascade Bachelor sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confXiBachelorPIDspecies{"confXiBachelorPIDspecies", std::vector<int>{o2::track::PID::Pion}, "Cascade Bachelor sel: Particles species for PID"};

    Configurable<bool> confXiRejectCompetingMass{"confXiRejectCompetingMass", false, "Switch on to reject Omegas (for Xi) or Xis (for Omegas)"};
    Configurable<float> confXiInvCompetingMassLowLimit{"confXiInvCompetingMassLowLimit", 1.66, "Lower limit of the cascade invariant mass for competing mass rejection"};
    Configurable<float> confXiInvCompetingMassUpLimit{"confXiInvCompetingMassUpLimit", 1.68, "Upper limit of the cascade invariant mass for competing mass rejection"};

    // Omega selection
    Configurable<float> confOmegaInvMassLowLimit{"confOmegaInvMassLowLimit", 1.2, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confOmegaInvMassUpLimit{"confOmegaInvMassUpLimit", 1.5, "Upper limit of the Cascade invariant mass"};
    // Cascade
    Configurable<std::vector<float>> confOmegaSign{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeSign, "confOmega"), std::vector<float>{-1, 1}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeSign, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaPtMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMin, "confOmega"), std::vector<float>{0.3f, 0.4f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaPtMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMax, "confOmega"), std::vector<float>{5.5f, 6.0f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaEtaMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeEtaMax, "confOmega"), std::vector<float>{0.8f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeEtaMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaDCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDCADaughMax, "confOmega"), std::vector<float>{1.f, 1.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaCPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeCPAMin, "confOmega"), std::vector<float>{0.99f, 0.95f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeCPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaTranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMin, "confOmega"), std::vector<float>{0.2f, 0.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaTranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMax, "confOmega"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confOmegaDecVtxMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDecVtxMax, "confOmega"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDecVtxMax, "Cascade selection: ")};

    // Cascade v0 daughters
    Configurable<std::vector<float>> confOmegaV0DCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "confOmega"), std::vector<float>{1.2f, 1.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confOmegaV0CPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0CPAMin, "confOmega"), std::vector<float>{0.99f, 0.995f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0CPAMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confOmegaV0TranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "confOmega"), std::vector<float>{0.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confOmegaV0TranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "confOmega"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confOmegaV0DCAtoPVMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "confOmega"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confOmegaV0DCAtoPVMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "confOmega"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "CascV0 selection: ")};
    Configurable<float> confOmegaV0InvMassLowLimit{"confOmegaV0InvMassLowLimit", 1.011461, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confOmegaV0InvMassUpLimit{"confOmegaV0InvMassUpLimit", 1.027461, "Upper limit of the Cascade invariant mass"};
    // Cascade Daughter Tracks
    Configurable<std::vector<float>> confOmegaV0ChildSign{"confOmegaV0ChildSign", std::vector<float>{-1, 1}, "CascV0 Child sel: Charge"};
    Configurable<std::vector<float>> confOmegaV0ChildPtMin{"confOmegaV0ChildPtMin", std::vector<float>{0.8f}, "CascV0 Child sel: min pt"};
    Configurable<std::vector<float>> confOmegaV0ChildEtaMax{"confOmegaV0ChildEtaMax", std::vector<float>{0.8f}, "CascV0 Child sel: max eta"};
    Configurable<std::vector<float>> confOmegaV0ChildTPCnClsMin{"confOmegaV0ChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "CascV0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confOmegaV0ChildDCAMin{"confOmegaV0ChildDCAMin", std::vector<float>{0.05f, 0.06f}, "CascV0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confOmegaV0ChildPIDnSigmaMax{"confOmegaV0ChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "CascV0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confOmegaV0ChildPIDspecies{"confOmegaV0ChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "CascV0 Child sel: Particles species for PID"};
    // Cascade Bachelor Track
    Configurable<std::vector<float>> confOmegaBachelorSign{"confOmegaBachelorSign", std::vector<float>{-1, 1}, "Cascade Bachelor sel: Charge"};
    Configurable<std::vector<float>> confOmegaBachelorPtMin{"confOmegaBachelorPtMin", std::vector<float>{0.8f}, "Cascade Bachelor sel: min pt"};
    Configurable<std::vector<float>> confOmegaBachelorEtaMax{"confOmegaBachelorEtaMax", std::vector<float>{0.8f}, "Cascade Bachelor sel: max eta"};
    Configurable<std::vector<float>> confOmegaBachelorTPCnClsMin{"confOmegaBachelorTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Cascade Bachelor sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confOmegaBachelorDCAMin{"confOmegaBachelorDCAMin", std::vector<float>{0.05f, 0.06f}, "Cascade Bachelor sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confOmegaBachelorPIDnSigmaMax{"confOmegaBachelorPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "Cascade Bachelor sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confOmegaBachelorPIDspecies{"confOmegaBachelorPIDspecies", std::vector<int>{o2::track::PID::Pion}, "Cascade Bachelor sel: Particles species for PID"};

    Configurable<bool> confOmegaRejectCompetingMass{"confOmegaRejectCompetingMass", false, "Switch on to reject Omegas (for Xi) or Xis (for Omegas)"};
    Configurable<float> confOmegaInvCompetingMassLowLimit{"confOmegaInvCompetingMassLowLimit", 1.66, "Lower limit of the cascade invariant mass for competing mass rejection"};
    Configurable<float> confOmegaInvCompetingMassUpLimit{"confOmegaInvCompetingMassUpLimit", 1.68, "Upper limit of the cascade invariant mass for competing mass rejection"};
  } confCascadeSel;

  // Resonances
  FemtoDreamResoSelection resoCuts;      // Phi
  FemtoDreamResoSelection resoCutsKStar; // KStar
  struct : ConfigurableGroup {
    std::string prefix = std::string("Resonance");

    // PhiCuts
    Configurable<float> confPhiInvMassLowLimit{"confPhiInvMassLowLimit", 0.9, "Lower limit of the Reso invariant mass"}; // 1.011461
    Configurable<float> confPhiInvMassUpLimit{"confPhiInvMassUpLimit", 1.15, "Upper limit of the Reso invariant mass"};  // 1.027461
    Configurable<bool> confDoLikeSignPhi{"confDoLikeSignPhi", false, "(De)activates likeSign histo filling"};
    Configurable<int> confMassQAPhiPart2PID{"confMassQAPhiPart2PID", o2::track::PID::Pion, "Daughter PID for massQAPart2 (bothPID2) histograms"};
    Configurable<std::vector<float>> confPhiSign{"confPhiSign", std::vector<float>{-1., 1.}, "Reso Sign selection"};

    Configurable<std::vector<float>> confPhiDaughterCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "confPhiDaughter"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterPtMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "confPhiDaughter"), std::vector<float>{0.1, 0.15, 0.2}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterPtMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "confPhiDaughter"), std::vector<float>{5.0, 4.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterEtaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "confPhiDaughter"), std::vector<float>{0.8, 0.85, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterTPCnClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "confPhiDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterTPCfClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "confPhiDaughter"), std::vector<float>{0.7, 0.8, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "confPhiDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "confPhiDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "confPhiDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Reso selection: ")};
    Configurable<std::vector<float>> confPhiDaughterPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "confPhiDaughter"), std::vector<float>{3.0, 2.5, 2.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Reso selection: ")};
    Configurable<std::vector<int>> confPhiDaughterPIDspecies{"confPhiDaughterPIDspecies", std::vector<int>{o2::track::PID::Kaon, o2::track::PID::Kaon}, "Reso Daughter sel: Particles species for PID"};
    Configurable<std::vector<float>> confPhiDaughterPTPCThr{"confPhiDaughterPTPCThr", std::vector<float>{0.5, 0.5}, "p_T (GeV/c) & pid dependent Thresholds for case distinction between TPC/TPCTOF"};

    // KStarCuts
    Configurable<float> confKstarInvMassLowLimit{"confKstarInvMassLowLimit", 0.9, "Lower limit of the Reso invariant mass"}; // 1.011461
    Configurable<float> confKstarInvMassUpLimit{"confKstarInvMassUpLimit", 1.15, "Upper limit of the Reso invariant mass"};  // 1.027461
    Configurable<bool> confDoLikeSignKstar{"confDoLikeSignKstar", false, "(De)activates likeSign histo filling"};
    Configurable<int> confMassQAKstarPart2PID{"confMassQAKstarPart2PID", o2::track::PID::Pion, "Daughter PID for massQAPart2 (bothPID2) histograms"};
    Configurable<std::vector<float>> confKstarSign{"confKstarSign", std::vector<float>{-1., 1.}, "Reso Sign selection"};

    Configurable<std::vector<float>> confKstarDaughterCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "confKstarDaughter"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterPtMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "confKstarDaughter"), std::vector<float>{0.1, 0.15, 0.2}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterPtMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "confKstarDaughter"), std::vector<float>{5.0, 4.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterEtaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "confKstarDaughter"), std::vector<float>{0.8, 0.85, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterTPCnClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "confKstarDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterTPCfClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "confKstarDaughter"), std::vector<float>{0.7, 0.8, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "confKstarDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "confKstarDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "confKstarDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Reso selection: ")};
    Configurable<std::vector<float>> confKstarDaughterPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "confKstarDaughter"), std::vector<float>{3.0, 2.5, 2.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Reso selection: ")};
    Configurable<std::vector<int>> confKstarDaughterPIDspecies{"confKstarDaughterPIDspecies", std::vector<int>{o2::track::PID::Kaon, o2::track::PID::Kaon}, "Reso Daughter sel: Particles species for PID"};
    Configurable<std::vector<float>> confKstarDaughterPTPCThr{"confKstarDaughterPTPCThr", std::vector<float>{0.5, 0.5}, "p_T (GeV/c) & pid dependent Thresholds for case distinction between TPC/TPCTOF"};
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
  HistogramRegistry qaRegistryV0{"QAHistosV0", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistryCascade{"QAHistosCascade", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistryReso{"QAHistosReso", {}, OutputObjHandlingPolicy::AnalysisObject};
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
    resoRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});
    cascadeRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});

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
      LambdaCuts.setSelection(V0Sel.confLambdaSign, femto_dream_v0_selection::kV0Sign, femtoDreamSelection::kEqual);
      LambdaCuts.setSelection(V0Sel.confLambdaPtMin, femto_dream_v0_selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaPtMax, femto_dream_v0_selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaEtaMax, femto_dream_v0_selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaDCADaughMax, femto_dream_v0_selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaCPAMin, femto_dream_v0_selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaTranRadMin, femto_dream_v0_selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaTranRadMax, femto_dream_v0_selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      LambdaCuts.setSelection(V0Sel.confLambdaDecVtxMax, femto_dream_v0_selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      LambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      LambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      LambdaCuts.setChildPIDSpecies(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildPIDspecies);
      LambdaCuts.setChildPIDSpecies(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildPIDspecies);
      LambdaCuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistryV0, &v0Registry);
      LambdaCuts.setInvMassLimits(V0Sel.confLambdaInvMassLowLimit, V0Sel.confLambdaInvMassUpLimit);
      LambdaCuts.setIsMother(true);

      LambdaCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kPosTrack, confTrkRejectNotPropagated);
      LambdaCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kNegTrack, confTrkRejectNotPropagated);

      LambdaCuts.setnSigmaPIDOffsetTPC(Track.confTrkPIDnSigmaOffsetTPC);
      LambdaCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kPosTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);
      LambdaCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kNegTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);

      if (V0Sel.confLambdaRejectKaons) {
        LambdaCuts.setKaonInvMassLimits(V0Sel.confLambdaInvKaonMassLowLimit, V0Sel.confLambdaInvKaonMassUpLimit);
      }
    }

    if (confIsActivateV0K0S) {
      K0SCuts.setSelection(V0Sel.confK0shortSign, femto_dream_v0_selection::kV0Sign, femtoDreamSelection::kEqual);
      K0SCuts.setSelection(V0Sel.confK0shortPtMin, femto_dream_v0_selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      K0SCuts.setSelection(V0Sel.confK0shortPtMax, femto_dream_v0_selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      K0SCuts.setSelection(V0Sel.confK0shortEtaMax, femto_dream_v0_selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      K0SCuts.setSelection(V0Sel.confK0shortDCADaughMax, femto_dream_v0_selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      K0SCuts.setSelection(V0Sel.confK0shortCPAMin, femto_dream_v0_selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      K0SCuts.setSelection(V0Sel.confK0shortTranRadMin, femto_dream_v0_selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      K0SCuts.setSelection(V0Sel.confK0shortTranRadMax, femto_dream_v0_selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      K0SCuts.setSelection(V0Sel.confK0shortDecVtxMax, femto_dream_v0_selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      K0SCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      K0SCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      K0SCuts.setChildPIDSpecies(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildPIDspecies);
      K0SCuts.setChildPIDSpecies(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildPIDspecies);
      K0SCuts.init<aod::femtodreamparticle::ParticleType::kV0K0Short, aod::femtodreamparticle::ParticleType::kV0K0ShortChild, aod::femtodreamparticle::cutContainerType>(&qaRegistryV0, &v0Registry);
      K0SCuts.setInvMassLimits(V0Sel.confK0shortInvMassLowLimit, V0Sel.confK0shortInvMassUpLimit);
      K0SCuts.setIsMother(false);

      K0SCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kPosTrack, confTrkRejectNotPropagated);
      K0SCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kNegTrack, confTrkRejectNotPropagated);

      K0SCuts.setnSigmaPIDOffsetTPC(Track.confTrkPIDnSigmaOffsetTPC);
      K0SCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kPosTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);
      K0SCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kNegTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);

      K0SCuts.setRejectLambda(V0Sel.confK0shortRejectLambdas);
      K0SCuts.setKaonInvMassLimits(V0Sel.confK0shortInvKaonMassLowLimit, V0Sel.confK0shortInvKaonMassUpLimit);
    }

    if (confIsActivateXi) {
      // Cascades
      xiCuts.setSelection(confCascadeSel.confXiSign, femtoDreamCascadeSelection::kCascadeSign, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeSign));
      xiCuts.setSelection(confCascadeSel.confXiPtMin, femtoDreamCascadeSelection::kCascadePtMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMin));
      xiCuts.setSelection(confCascadeSel.confXiPtMax, femtoDreamCascadeSelection::kCascadePtMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMax));
      xiCuts.setSelection(confCascadeSel.confXiEtaMax, femtoDreamCascadeSelection::kCascadeEtaMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeEtaMax));
      xiCuts.setSelection(confCascadeSel.confXiDCADaughMax, femtoDreamCascadeSelection::kCascadeDCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDCADaughMax));
      xiCuts.setSelection(confCascadeSel.confXiCPAMin, femtoDreamCascadeSelection::kCascadeCPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeCPAMin));
      xiCuts.setSelection(confCascadeSel.confXiTranRadMin, femtoDreamCascadeSelection::kCascadeTranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMin));
      xiCuts.setSelection(confCascadeSel.confXiTranRadMax, femtoDreamCascadeSelection::kCascadeTranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMax));
      xiCuts.setSelection(confCascadeSel.confXiDecVtxMax, femtoDreamCascadeSelection::kCascadeDecVtxMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDecVtxMax));
      // Cascade v0
      xiCuts.setSelection(confCascadeSel.confXiV0DCADaughMax, femtoDreamCascadeSelection::kCascadeV0DCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCADaughMax));
      xiCuts.setSelection(confCascadeSel.confXiV0CPAMin, femtoDreamCascadeSelection::kCascadeV0CPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0CPAMin));
      xiCuts.setSelection(confCascadeSel.confXiV0TranRadMin, femtoDreamCascadeSelection::kCascadeV0TranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMin));
      xiCuts.setSelection(confCascadeSel.confXiV0TranRadMax, femtoDreamCascadeSelection::kCascadeV0TranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMax));
      xiCuts.setSelection(confCascadeSel.confXiV0DCAtoPVMin, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin));
      xiCuts.setSelection(confCascadeSel.confXiV0DCAtoPVMax, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax));

      // Cascade Daughter Tracks
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confXiV0ChildPIDspecies);

      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confXiV0ChildPIDspecies);

      // Cascade Bachelor Track
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      xiCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      xiCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confXiBachelorPIDspecies);

      xiCuts.init<aod::femtodreamparticle::ParticleType::kCascade, aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeBachelor, aod::femtodreamparticle::cutContainerType>(&qaRegistryCascade, &cascadeRegistry, false);
      xiCuts.setInvMassLimits(confCascadeSel.confXiInvMassLowLimit, confCascadeSel.confXiInvMassUpLimit);
      xiCuts.setV0InvMassLimits(confCascadeSel.confXiV0InvMassLowLimit, confCascadeSel.confXiV0InvMassUpLimit);
      if (confCascadeSel.confXiRejectCompetingMass) {
        xiCuts.setCompetingInvMassLimits(confCascadeSel.confXiInvCompetingMassLowLimit, confCascadeSel.confXiInvCompetingMassUpLimit);
      }
    }

    if (confIsActivateOmega) {
      // Cascades
      omegaCuts.setSelection(confCascadeSel.confOmegaSign, femtoDreamCascadeSelection::kCascadeSign, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeSign));
      omegaCuts.setSelection(confCascadeSel.confOmegaPtMin, femtoDreamCascadeSelection::kCascadePtMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaPtMax, femtoDreamCascadeSelection::kCascadePtMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaEtaMax, femtoDreamCascadeSelection::kCascadeEtaMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeEtaMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaDCADaughMax, femtoDreamCascadeSelection::kCascadeDCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDCADaughMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaCPAMin, femtoDreamCascadeSelection::kCascadeCPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeCPAMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaTranRadMin, femtoDreamCascadeSelection::kCascadeTranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaTranRadMax, femtoDreamCascadeSelection::kCascadeTranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaDecVtxMax, femtoDreamCascadeSelection::kCascadeDecVtxMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDecVtxMax));
      // Cascade v0
      omegaCuts.setSelection(confCascadeSel.confOmegaV0DCADaughMax, femtoDreamCascadeSelection::kCascadeV0DCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCADaughMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaV0CPAMin, femtoDreamCascadeSelection::kCascadeV0CPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0CPAMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaV0TranRadMin, femtoDreamCascadeSelection::kCascadeV0TranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaV0TranRadMax, femtoDreamCascadeSelection::kCascadeV0TranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMax));
      omegaCuts.setSelection(confCascadeSel.confOmegaV0DCAtoPVMin, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin));
      omegaCuts.setSelection(confCascadeSel.confOmegaV0DCAtoPVMax, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax));

      // Cascade Daughter Tracks
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kPosTrack, confCascadeSel.confOmegaV0ChildPIDspecies);

      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kNegTrack, confCascadeSel.confOmegaV0ChildPIDspecies);

      // Cascade Bachelor Track
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      omegaCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      omegaCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kBachTrack, confCascadeSel.confOmegaBachelorPIDspecies);

      omegaCuts.init<aod::femtodreamparticle::ParticleType::kOmega, aod::femtodreamparticle::ParticleType::kOmegaV0Child, aod::femtodreamparticle::ParticleType::kOmegaBachelor, aod::femtodreamparticle::cutContainerType>(&qaRegistryCascade, &cascadeRegistry, true);
      omegaCuts.setInvMassLimits(confCascadeSel.confOmegaInvMassLowLimit, confCascadeSel.confOmegaInvMassUpLimit);
      omegaCuts.setV0InvMassLimits(confCascadeSel.confOmegaV0InvMassLowLimit, confCascadeSel.confOmegaV0InvMassUpLimit);
      if (confCascadeSel.confOmegaRejectCompetingMass) {
        omegaCuts.setCompetingInvMassLimits(confCascadeSel.confOmegaInvCompetingMassLowLimit, confCascadeSel.confOmegaInvCompetingMassUpLimit);
      }
    }

    if (confIsActivatePhi.value) {
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCuts.init<aod::femtodreamparticle::ParticleType::kReso,
                    aod::femtodreamparticle::ParticleType::kResoChild>(&qaRegistryReso, &resoRegistry);

      resoCuts.assign(Resonance.confPhiDaughterPTPCThr); // assigns Configurable value to class member
      resoCuts.setDaughterPIDSpecies(femto_dream_reso_selection::kPosdaugh, Resonance.confPhiDaughterPIDspecies);
      resoCuts.setDaughterPIDSpecies(femto_dream_reso_selection::kNegdaugh, Resonance.confPhiDaughterPIDspecies);
      resoCuts.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kPosdaugh, 0.f, 0.f);
      resoCuts.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kNegdaugh, 0.f, 0.f);

      resoCuts.setSelection(Resonance.confPhiSign, femto_dream_reso_selection::kResoSign, femtoDreamSelection::kEqual);
    }

    if (confIsActivateKStar.value) {
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCutsKStar.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCutsKStar.init<aod::femtodreamparticle::ParticleType::kResoKStar, /// chose this particle type, since there is no kResoKstar-Particle type (implementing it would only serve the naming of the producer histos)
                         aod::femtodreamparticle::ParticleType::kResoKStarChild>(&qaRegistryReso, &resoRegistry);

      resoCutsKStar.assign(Resonance.confKstarDaughterPTPCThr); // assigns Configurable value to class member
      resoCutsKStar.setDaughterPIDSpecies(femto_dream_reso_selection::kPosdaugh, Resonance.confKstarDaughterPIDspecies);
      resoCutsKStar.setDaughterPIDSpecies(femto_dream_reso_selection::kNegdaugh, Resonance.confKstarDaughterPIDspecies);
      resoCutsKStar.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kPosdaugh, 0.f, 0.f);
      resoCutsKStar.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kNegdaugh, 0.f, 0.f);

      resoCutsKStar.setSelection(Resonance.confKstarSign, femto_dream_reso_selection::kResoSign, femtoDreamSelection::kEqual);
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

  template <aod::femtodreamparticle::ParticleType part, typename P, typename V, typename E>
  void fillLikeSign(P const& sliceDaughters, FemtoDreamResoSelection& resoSelectionObject, V const& pidVector, E const& extraPID)
  {
    for (const auto& track1 : sliceDaughters) {
      if (!resoSelectionObject.daughterSelectionPos(track1) || !resoSelectionObject.isSelectedMinimalPIDPos(track1, pidVector)) {
        continue;
      }
      for (const auto& track2 : sliceDaughters) {
        if (!resoSelectionObject.daughterSelectionPos(track2) || !resoSelectionObject.isSelectedMinimalPIDPos(track2, pidVector)) {
          continue;
        }
        bool resoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
        if ((pidVector.size() > 1) && (pidVector[0] != pidVector[1])) {
          auto [isNormal, WrongCombination] = resoSelectionObject.checkCombination(track1, track2, pidVector);
          if (WrongCombination) {
            continue;
          }
          resoIsNotAnti = isNormal;
        }
        /// Resos, where both daughters have the same PID are defaulted to sign 1. and resoIsNotAnti = true
        resoSelectionObject.fillLikeSignHistos<part>(track1, track2, pidVector, extraPID, resoIsNotAnti);
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

    // minimal selection
    bool bEmptyColl = true;
    if (!colCuts.isEmptyCollision(col, tracks, trackCuts)) {
      bEmptyColl = false;
    }
    if (bEmptyColl && confIsActivateXi.value) {
      if (!colCuts.isCollisionWithoutTrkCasc(col, fullCascades, xiCuts, tracks)) {
        bEmptyColl = false;
      }
    }
    if (bEmptyColl && confIsActivateOmega.value) {
      if (!colCuts.isCollisionWithoutTrkCasc(col, fullCascades, omegaCuts, tracks)) {
        bEmptyColl = false;
      }
    }
    if (bEmptyColl && confIsActivateV0.value) {
      if (!colCuts.isEmptyCollision(col, fullV0s, LambdaCuts, tracks)) {
        bEmptyColl = false;
      }
    }
    if (bEmptyColl && confIsActivateV0K0S.value) {
      if (!colCuts.isEmptyCollision(col, fullV0s, K0SCuts, tracks)) {
        bEmptyColl = false;
      }
    }
    if (bEmptyColl) {
      return;
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

    // v01 (Lambdas)
    if (confIsActivateV0.value) {
      for (const auto& v0 : fullV0s) {

        auto postrack = v0.template posTrack_as<TrackType>();
        auto negtrack = v0.template negTrack_as<TrackType>();
        ///\tocheck funnily enough if we apply the filter the
        /// sign of Pos and Neg track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        LambdaCuts.fillLambdaQA<aod::femtodreamparticle::ParticleType::kV0>(col, v0, postrack, negtrack);

        if (!LambdaCuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }

        // if (confRejectITSHitandTOFMissing) {
        // Uncomment only when TOF timing is solved
        // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack,
        // TrackSelection::TrackCuts::kITSHits); bool itsHit =
        // o2PhysicsTrackSelection->IsSelected(negtrack,
        // TrackSelection::TrackCuts::kITSHits);
        // }

        LambdaCuts.fillQA<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 = LambdaCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);

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

        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    aod::femtodreamparticle::ParticleType::kV0,
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    v0.mLambda(),
                    v0.mAntiLambda());
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

    // v02 (K0S)
    if (confIsActivateV0K0S.value) {
      for (const auto& v0 : fullV0s) {

        auto postrack = v0.template posTrack_as<TrackType>();
        auto negtrack = v0.template negTrack_as<TrackType>();
        ///\tocheck funnily enough if we apply the filter the
        /// sign of Pos and Neg track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        K0SCuts.fillLambdaQA<aod::femtodreamparticle::ParticleType::kV0K0Short>(col, v0, postrack, negtrack);

        if (!K0SCuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }

        // if (confRejectITSHitandTOFMissing) {
        // Uncomment only when TOF timing is solved
        // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack,
        // TrackSelection::TrackCuts::kITSHits); bool itsHit =
        // o2PhysicsTrackSelection->IsSelected(negtrack,
        // TrackSelection::TrackCuts::kITSHits);
        // }

        K0SCuts.fillQA<aod::femtodreamparticle::ParticleType::kV0K0Short, aod::femtodreamparticle::ParticleType::kV0K0ShortChild>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 = K0SCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);

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

        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    aod::femtodreamparticle::ParticleType::kV0K0Short,
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    v0.mK0Short(),
                    v0.mK0Short());
        if (confIsDebug.value) {
          fillDebugParticle<true, false>(postrack); // QA for positive daughter
          fillDebugParticle<true, false>(negtrack); // QA for negative daughter
          fillDebugParticle<false, false>(v0);      // QA for v0
        }
        if constexpr (isMC) {
          fillMCParticle(col, v0, o2::aod::femtodreamparticle::ParticleType::kV0); /// change particle Type here as well?
        }
      }
    }

    // Cascades (Xi and Omega)
    if (confIsActivateXi.value || confIsActivateOmega.value) {
      for (const auto& casc : fullCascades) {
        // get the daughter tracks
        const auto& posTrackCasc = casc.template posTrack_as<TrackType>();
        const auto& negTrackCasc = casc.template negTrack_as<TrackType>();
        const auto& bachTrackCasc = casc.template bachelor_as<TrackType>();

        if (confIsActivateXi.value) {
          // xiCuts.fillQA<0, aod::femtodreamparticle::ParticleType::kCascade, aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeBachelor>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

          if (xiCuts.isSelectedMinimal(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc)) {

            // xiCuts.fillQA<1, aod::femtodreamparticle::ParticleType::kCascade, aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeBachelor>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);
            //  auto cutContainerCasc = xiCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, v0daugh, posTrackCasc, negTrackCasc, bachTrackCasc);
            auto cutContainerCasc = xiCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

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
                        casc.mXi(),
                        casc.mLambda());
            // TODO: include here MC filling
            //------

            if (confIsDebug.value) {
              fillDebugParticle<true, false>(posTrackCasc);  // QA for positive daughter
              fillDebugParticle<true, false>(negTrackCasc);  // QA for negative daughter
              fillDebugParticle<true, false>(bachTrackCasc); // QA for negative daughter
              fillDebugCascade(casc, col);                   // QA for Cascade
            }

            // continue;
          } // if xiCuts.isSelectedMinimal
        } // if confIsActivateXi
        /*
        if (confIsActivateOmega.value){
          omegaCuts.fillQA<0, aod::femtodreamparticle::ParticleType::kOmega, aod::femtodreamparticle::ParticleType::kOmegaV0Child, aod::femtodreamparticle::ParticleType::kOmegaBachelor>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

          if (omegaCuts.isSelectedMinimal(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc)) {

            omegaCuts.fillQA<1, aod::femtodreamparticle::ParticleType::kOmega, aod::femtodreamparticle::ParticleType::kOmegaV0Child, aod::femtodreamparticle::ParticleType::kOmegaBachelor>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);
            // auto cutContainerCasc = xiCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, v0daugh, posTrackCasc, negTrackCasc, bachTrackCasc);
            auto cutContainerCasc = omegaCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

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
                        aod::femtodreamparticle::ParticleType::kOmegaV0Child,
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
                        aod::femtodreamparticle::ParticleType::kOmegaV0Child,
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
                        aod::femtodreamparticle::ParticleType::kOmegaBachelor,
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
            std::vector<int> indexCascadeChildID = {rowOfPosCascadeTrack, rowOfNegCascadeTrack, rowOfBachelorCascadeTrack};
            outputParts(outputCollision.lastIndex(),
                        casc.pt(),
                        casc.eta(),
                        casc.phi(),
                        aod::femtodreamparticle::ParticleType::kOmega,
                        cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kCascade),
                        0,
                        casc.casccosPA(col.posX(), col.posY(), col.posZ()),
                        indexCascadeChildID,
                        casc.mOmega(),
                        casc.mLambda());
            // TODO: include here MC filling
            //------

            if (confIsDebug.value) {
              fillDebugParticle<true, false>(posTrackCasc);  // QA for positive daughter
              fillDebugParticle<true, false>(negTrackCasc);  // QA for negative daughter
              fillDebugParticle<true, false>(bachTrackCasc); // QA for negative daughter
              fillDebugCascade(casc, col);                   // QA for Cascade
            }


            //continue;
          } //if omegaCuts.isSelectedMinimal
        } //if confIsActivateOmega
        */
      } // loop over cascades
    } // at least one cascade active

    if (confIsActivatePhi.value) {

      resoCuts.updateSigmaPIDMax();

      auto slicePosdaugh = daughter1.sliceByCached(aod::track::collisionId, col.globalIndex(), cache); // o2::framework defined in AnalysisHelper.h
      auto sliceNegdaugh = daughter2.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);

      for (const auto& track1 : slicePosdaugh) {
        if (!resoCuts.daughterSelectionPos(track1) || !resoCuts.isSelectedMinimalPIDPos(track1, Resonance.confPhiDaughterPIDspecies.value)) {
          continue;
        }

        for (const auto& track2 : sliceNegdaugh) {
          if (!resoCuts.daughterSelectionNeg(track2) || !resoCuts.isSelectedMinimalPIDNeg(track2, Resonance.confPhiDaughterPIDspecies.value)) {
            continue; /// loosest cuts for track2
          }

          /// This only works for the case where the mass of opposite charged particles are the same (for example K+/K- have same mass)
          float massPart1 = o2::track::PID::getMass(Resonance.confPhiDaughterPIDspecies.value[0]);
          float massPart2 = massPart1;
          if (Resonance.confPhiDaughterPIDspecies->size() > 1)
            massPart2 = o2::track::PID::getMass(Resonance.confPhiDaughterPIDspecies.value[1]);

          /// Resonance
          ROOT::Math::PtEtaPhiMVector tempD1(track1.pt(), track1.eta(), track1.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempD2(track2.pt(), track2.eta(), track2.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempReso = tempD1 + tempD2;
          /// Anti-resonance
          ROOT::Math::PtEtaPhiMVector tempDA1(track1.pt(), track1.eta(), track1.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempDA2(track2.pt(), track2.eta(), track2.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempAntiReso = tempDA1 + tempDA2;

          bool resoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
          float resoSign = 1.;
          if ((Resonance.confPhiDaughterPIDspecies->size() > 1) && (Resonance.confPhiDaughterPIDspecies.value[0] != Resonance.confPhiDaughterPIDspecies.value[1])) {
            auto [isNormal, WrongCombination] = resoCuts.checkCombination(track1, track2, Resonance.confPhiDaughterPIDspecies.value);
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
            resoCuts.fillResoQA<aod::femtodreamparticle::ParticleType::kReso>(track1, track2, true, tempReso.M(), tempAntiReso.M(), Resonance.confPhiDaughterPIDspecies.value, Resonance.confMassQAPhiPart2PID.value);
            if (!(tempReso.M() > Resonance.confPhiInvMassLowLimit.value && tempReso.M() < Resonance.confPhiInvMassUpLimit.value))
              continue;
            resoCuts.fillMassSelectedQA<aod::femtodreamparticle::ParticleType::kReso>(tempReso.M(), true);
          } else {
            resoCuts.fillResoQA<aod::femtodreamparticle::ParticleType::kReso>(track1, track2, false, tempAntiReso.M(), tempReso.M(), Resonance.confPhiDaughterPIDspecies.value, Resonance.confMassQAPhiPart2PID.value);
            if (!(tempAntiReso.M() > Resonance.confPhiInvMassLowLimit.value && tempAntiReso.M() < Resonance.confPhiInvMassUpLimit.value))
              continue;
            resoCuts.fillMassSelectedQA<aod::femtodreamparticle::ParticleType::kReso>(tempAntiReso.M(), false);
          }

          resoCuts.fillQA<aod::femtodreamparticle::ParticleType::kResoChild,
                          aod::femtodreamparticle::TrackType::kPosChild,
                          aod::femtodreamparticle::TrackType::kNegChild>(track1, track2);

          auto type = resoCuts.getType<aod::femtodreamparticle::kReso>(track1, track2, resoIsNotAnti); //   kResoPosdaughTPC_NegdaughTPC
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

      if (Resonance.confDoLikeSignPhi.value) {
        fillLikeSign<aod::femtodreamparticle::ParticleType::kReso>(slicePosdaugh, resoCuts, Resonance.confPhiDaughterPIDspecies.value, Resonance.confMassQAPhiPart2PID.value);
        fillLikeSign<aod::femtodreamparticle::ParticleType::kReso>(sliceNegdaugh, resoCuts, Resonance.confPhiDaughterPIDspecies.value, Resonance.confMassQAPhiPart2PID.value);
      }
    } // if (confIsActivatePhi.value)

    if (confIsActivateKStar.value) {

      resoCutsKStar.updateSigmaPIDMax();

      auto slicePosdaugh = daughter1.sliceByCached(aod::track::collisionId, col.globalIndex(), cache); // o2::framework defined in AnalysisHelper.h
      auto sliceNegdaugh = daughter2.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);

      for (const auto& track1 : slicePosdaugh) {
        if (!resoCutsKStar.daughterSelectionPos(track1) || !resoCutsKStar.isSelectedMinimalPIDPos(track1, Resonance.confKstarDaughterPIDspecies.value)) {
          continue;
        }

        for (const auto& track2 : sliceNegdaugh) {
          if (!resoCutsKStar.daughterSelectionNeg(track2) || !resoCutsKStar.isSelectedMinimalPIDNeg(track2, Resonance.confKstarDaughterPIDspecies.value)) {
            continue; /// loosest cuts for track2
          }

          /// This only works for the case where the mass of opposite charged particles are the same (for example K+/K- have same mass)
          float massPart1 = o2::track::PID::getMass(Resonance.confKstarDaughterPIDspecies.value[0]);
          float massPart2 = massPart1;
          if (Resonance.confKstarDaughterPIDspecies->size() > 1)
            massPart2 = o2::track::PID::getMass(Resonance.confKstarDaughterPIDspecies.value[1]);

          /// Resonance
          ROOT::Math::PtEtaPhiMVector tempD1(track1.pt(), track1.eta(), track1.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempD2(track2.pt(), track2.eta(), track2.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempReso = tempD1 + tempD2;
          /// Anti-resonance
          ROOT::Math::PtEtaPhiMVector tempDA1(track1.pt(), track1.eta(), track1.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempDA2(track2.pt(), track2.eta(), track2.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempAntiReso = tempDA1 + tempDA2;

          bool resoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
          float resoSign = 1.;
          if ((Resonance.confKstarDaughterPIDspecies->size() > 1) && (Resonance.confKstarDaughterPIDspecies.value[0] != Resonance.confKstarDaughterPIDspecies.value[1])) {
            auto [isNormal, WrongCombination] = resoCutsKStar.checkCombination(track1, track2, Resonance.confKstarDaughterPIDspecies.value);
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
            resoCutsKStar.fillResoQA<aod::femtodreamparticle::ParticleType::kResoKStar>(track1, track2, true, tempReso.M(), tempAntiReso.M(), Resonance.confKstarDaughterPIDspecies.value, Resonance.confMassQAKstarPart2PID.value);
            if (!(tempReso.M() > Resonance.confKstarInvMassLowLimit.value && tempReso.M() < Resonance.confKstarInvMassUpLimit.value))
              continue;
            resoCutsKStar.fillMassSelectedQA<aod::femtodreamparticle::ParticleType::kResoKStar>(tempReso.M(), true);
          } else {
            resoCutsKStar.fillResoQA<aod::femtodreamparticle::ParticleType::kResoKStar>(track1, track2, false, tempAntiReso.M(), tempReso.M(), Resonance.confKstarDaughterPIDspecies.value, Resonance.confMassQAKstarPart2PID.value);
            if (!(tempAntiReso.M() > Resonance.confKstarInvMassLowLimit.value && tempAntiReso.M() < Resonance.confKstarInvMassUpLimit.value))
              continue;
            resoCutsKStar.fillMassSelectedQA<aod::femtodreamparticle::ParticleType::kResoKStar>(tempAntiReso.M(), false);
          }

          resoCutsKStar.fillQA<aod::femtodreamparticle::ParticleType::kResoKStarChild,
                               aod::femtodreamparticle::TrackType::kPosChild,
                               aod::femtodreamparticle::TrackType::kNegChild>(track1, track2);

          auto type = resoCutsKStar.getType<aod::femtodreamparticle::kResoKStar>(track1, track2, resoIsNotAnti); //   kResoPosdaughTPC_NegdaughTPC
                                                                                                                 //   kResoPosdaughTPC_NegdaughTOF
                                                                                                                 //   kResoPosdaughTPC_NegdaughTPC
                                                                                                                 //   kResoPosdaughTOF_NegdaughTOF as possible output

          auto bitmask = resoCutsKStar.getCutContainer<aod::femtodreamparticle::cutContainerType>(track1, track2, resoSign);

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

      if (Resonance.confDoLikeSignKstar.value) {
        fillLikeSign<aod::femtodreamparticle::ParticleType::kResoKStar>(slicePosdaugh, resoCutsKStar, Resonance.confKstarDaughterPIDspecies.value, Resonance.confMassQAKstarPart2PID.value);
        fillLikeSign<aod::femtodreamparticle::ParticleType::kResoKStar>(sliceNegdaugh, resoCutsKStar, Resonance.confKstarDaughterPIDspecies.value, Resonance.confMassQAKstarPart2PID.value);
      }
    } // if (confIsActivateKStar.value)
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
