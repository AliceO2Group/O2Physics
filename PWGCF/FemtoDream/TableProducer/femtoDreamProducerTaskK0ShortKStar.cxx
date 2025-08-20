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

/// \file femtoDreamProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include "FairLogger.h" // delete after debugging

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCascadeSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamResoSelectionKStar.h"
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
#include <Framework/Logger.h> // delete after debugging

#include "Math/Vector4D.h"
#include "TMath.h"

#include <fairlogger/Logger.h>

#include <fstream>  // delete after debugging
#include <iostream> // delete in the end
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
using FemtoFullCollision_noCent = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollision_CentPbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
using FemtoFullCollision_noCent_MC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullCollisionMC_CentPbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::McCollisionLabels>::iterator;
using FemtoFullMCgenCollisions = soa::Join<aod::McCollisions, MultsExtraMC>;
using FemtoFullMCgenCollision = FemtoFullMCgenCollisions::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
            aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
            aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;

} // namespace o2::aod

namespace softwareTriggers
{
static const int nTriggers = 6;
static const std::vector<std::string> triggerNames{"fPPP", "fPPL", "fPLL", "fLLL", "fPD", "fLD"};
static const float triggerSwitches[1][nTriggers]{
  {0, 0, 0, 0, 0, 0}};
} // namespace softwareTriggers

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

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  Configurable<bool> ConfUseItsPid{"ConfUseItsPid", false, "Enable Debug tables"};
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"}; // true?
  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};
  /// Event cuts
  FemtoDreamCollisionSelection colCuts;
  // Event cuts - Triggers
  Configurable<bool> ConfEnableTriggerSelection{"ConfEnableTriggerSelection", false, "Should the trigger selection be enabled for collisions?"};
  Configurable<LabeledArray<float>> ConfTriggerSwitches{
    "ConfTriggerSwitches",
    {softwareTriggers::triggerSwitches[0], 1, softwareTriggers::nTriggers, std::vector<std::string>{"Switch"}, softwareTriggers::triggerNames},
    "Turn on which trigger should be checked for recorded events to pass selection"};
  Configurable<std::string> ConfBaseCCDBPathForTriggers{"ConfBaseCCDBPathForTriggers", "Users/m/mpuccio/EventFiltering/OTS/Chunked/", "Provide ccdb path for trigger table; default - trigger coordination"};

  // Event cuts - usual selection criteria
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfEvtAddOfflineCheck{"ConfEvtAddOfflineCheck", false, "Evt sel: additional checks for offline selection (not part of sel8 yet)"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", true, "Activate filling of V0 into femtodream tables"};
  Configurable<bool> ConfIsActivateReso{"ConfIsActivateReso", true, "Activate filling of sl Resonances into femtodream tables"};
  Configurable<bool> ConfIsActivatePhi{"ConfIsActivatePhi", true, "Activates cuts on Phi's and fills tables"};

  Configurable<bool> ConfTrkRejectNotPropagated{"ConfTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> ConfRejectITSHitandTOFMissing{ "ConfRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};
  Configurable<int> ConfTrkPDGCode{"ConfTrkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  FemtoDreamTrackSelection trackCuts;
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track");
    Configurable<std::vector<float>> ConfTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.1f, 0.15f, 0.2f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{4.4f, 4.6f, 4.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.85f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 90.f, 100.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.2f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
    Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"}; // set to zero for run3 or so
    Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
    Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
    // missing DCA configurable?? because implemented in TrackSelection.h
  } Track;

  FemtoDreamV0Selection v0Cuts;
  Configurable<std::vector<float>> ConfV0Sign{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0Sign, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMin{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0pTMin, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0pTMin, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0pTMax, "ConfV0"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0pTMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0EtaMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0etaMax, "ConfV0"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0etaMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0DCADaughMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0DCADaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0DCADaughMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0CPAMin{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0CPAMin, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0CPAMin, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0TranRadMin{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0TranRadMin, "ConfV0"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0TranRadMin, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0TranRadMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0TranRadMax, "ConfV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0TranRadMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0DecVtxMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0DecVtxMax, "ConfV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0DecVtxMax, "V0 selection: ")};

  Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};
  Configurable<bool> ConfV0RejectKaons{"ConfV0RejectKaons", false, "Switch to reject kaons"};
  Configurable<bool> ConfV0RejectLambdas{"ConfV0RejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
  Configurable<float> ConfV0InvKaonMassLowLimit{"ConfV0InvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> ConfV0InvKaonMassUpLimit{"ConfV0InvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};
  Configurable<bool> ConfV0MotherIsLambda{"ConfV0MotherIsLambda", true, "True: Lambda, False: K0Short"};

  Configurable<std::vector<float>> ConfChildCharge{"ConfChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
  Configurable<std::vector<float>> ConfChildEtaMax{"ConfChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
  Configurable<std::vector<float>> ConfChildTPCnClsMin{"ConfChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
  Configurable<std::vector<float>> ConfChildDCAMin{"ConfChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> ConfChildPIDnSigmaMax{"ConfChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
  Configurable<std::vector<int>> ConfChildPIDspecies{"ConfChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

  // Resonances
  FemtoDreamResoSelection resoCuts;
  struct : ConfigurableGroup {
    std::string prefix = std::string("Resonance");
    Configurable<bool> ConfThreshold{"ConfThreshold", true, "Enable threshold selection between TPC/TPCTOF"};

    Configurable<float> ConfResoInvMassLowLimit{"ConfResoInvMassLowLimit", 0.9, "Lower limit of the Reso invariant mass"}; // 1.011461
    Configurable<float> ConfResoInvMassUpLimit{"ConfResoInvMassUpLimit", 1.15, "Upper limit of the Reso invariant mass"};  // 1.027461

    Configurable<std::vector<float>> ConfDaughterCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfDaughter"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterPtMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfDaughter"), std::vector<float>{0.1, 0.15, 0.2}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterPtMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "ConfDaughter"), std::vector<float>{5.0, 4.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterEtaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "ConfDaughter"), std::vector<float>{0.8, 0.85, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterTPCnClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "ConfDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterTPCfClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "ConfDaughter"), std::vector<float>{0.7, 0.8, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "ConfDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "ConfDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "ConfDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Reso selection: ")};
    Configurable<std::vector<float>> ConfDaughterPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "ConfDaughter"), std::vector<float>{3.0, 2.5, 2.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Reso selection: ")};
    // Configurable<float> ConfResoMassUp{"ConfResoMassUp", 0.52, "Upper limit for the mass selection of the daughters"};
    // Configurable<float> ConfResoMassLow{"ConfResoMassLow", 0.48, "Lower limit for the mass selection of the daughters"};
    Configurable<std::vector<int>> ConfDaughterPIDspecies{"ConfDaughterPIDspecies", std::vector<int>{o2::track::PID::Kaon}, "Reso Daughter sel: Particles species for PID"};
    Configurable<float> ConfResoPTPCThr{"ConfResoPTPCThr", 0.5, "p_T (GeV/c)Threshold for case distinction between TPC/TPCTOF"};
    Configurable<std::vector<float>> ConfResoSign{"ConfResoSign", std::vector<float>{-1., 1.}, "Reso Sign selection"};
    Configurable<int> ConfResoMotherID{"ConfResoMotherID", static_cast<int>(femtoDreamResoSelection::kPhi), "ResoID of Mother [0: Phi, 1: KStar]"};
  } Resonance;

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  /// General options
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfTrkMinChi2PerClusterTPC{"ConfTrkMinChi2PerClusterTPC", 0.f, "Lower limit for chi2 of TPC; currently for testing only"};
    Configurable<float> ConfTrkMaxChi2PerClusterTPC{"ConfTrkMaxChi2PerClusterTPC", 1000.f, "Upper limit for chi2 of TPC; currently for testing only"};
    Configurable<float> ConfTrkMaxChi2PerClusterITS{"ConfTrkMaxChi2PerClusterITS", 1000.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default
    Configurable<bool> ConfTrkTPCRefit{"ConfTrkTPCRefit", false, "True: require TPC refit"};
    Configurable<bool> ConfTrkITSRefit{"ConfTrkITSRefit", false, "True: require ITS refit"};

  } OptionTrackSpecialSelections;

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry TrackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry V0Registry{"V0", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry ResoRegistry{"Reso", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float mMagField;
  std::string zorroTriggerNames = "";
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  RCTFlagsChecker rctChecker;

  void init(InitContext&)
  {
    if (doprocessData == false && doprocessData_noCentrality == false && doprocessData_CentPbPb == false && doprocessMC == false && doprocessMC_noCentrality == false && doprocessMC_CentPbPb == false) {
      LOGF(fatal, "Neither processData nor processMC enabled. Please choose one.");
    }
    if ((doprocessData == true && doprocessMC == true) || (doprocessData == true && doprocessMC_noCentrality == true) || (doprocessMC == true && doprocessMC_noCentrality == true) || (doprocessData_noCentrality == true && doprocessData == true) || (doprocessData_noCentrality == true && doprocessMC == true) || (doprocessData_noCentrality == true && doprocessMC_noCentrality == true) || (doprocessData_CentPbPb == true && doprocessData == true) || (doprocessData_CentPbPb == true && doprocessData_noCentrality == true) || (doprocessData_CentPbPb == true && doprocessMC == true) || (doprocessData_CentPbPb == true && doprocessMC_noCentrality == true) || (doprocessData_CentPbPb == true && doprocessMC_CentPbPb == true)) {
      LOGF(fatal,
           "Cannot enable more than one process switch at the same time. "
           "Please choose one.");
    }

    int CutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    TrackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});
    TrackRegistry.add("AnalysisQA/Chi2ITSTPCperCluster", "; ITS_Chi2; TPC_Chi2", kTH2F, {{100, 0, 50}, {100, 0, 20}});
    TrackRegistry.add("AnalysisQA/RefitITSTPC", "; ITS_Refit; TPC_Refit", kTH2F, {{2, 0, 2}, {2, 0, 2}});
    TrackRegistry.add("AnalysisQA/getGenStatusCode", "; Bit; Entries", kTH1F, {{200, 0, 200}});
    TrackRegistry.add("AnalysisQA/getProcess", "; Bit; Entries", kTH1F, {{200, 0, 200}});
    TrackRegistry.add("AnalysisQA/Mother", "; Bit; Entries", kTH1F, {{4000, -4000, 4000}});
    TrackRegistry.add("AnalysisQA/Particle", "; Bit; Entries", kTH1F, {{4000, -4000, 4000}});
    V0Registry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});

    ResoRegistry.add("AnalysisQA/Reso/InvMass", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    ResoRegistry.add("AnalysisQA/Reso/InvMassAnti", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    ResoRegistry.add("AnalysisQA/Reso/InvMass_phi_selected", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    ResoRegistry.add("AnalysisQA/Reso/InvMassAnti_phi_selected", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});

    ResoRegistry.add("AnalysisQA/Reso/Daughter1/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter1/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter1/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter1/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/Daughter1/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/Daughter2/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter2/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter2/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    ResoRegistry.add("AnalysisQA/Reso/Daughter2/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/Daughter2/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here

    ResoRegistry.add("AnalysisQA/Reso/Pt_posdaughter_selected", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    ResoRegistry.add("AnalysisQA/Reso/Eta_posdaughter_selected", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    ResoRegistry.add("AnalysisQA/Reso/Phi_posdaughter_selected", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    ResoRegistry.add("AnalysisQA/Reso/DCAxy_posdaughter_selected", "dcaXY of all processed tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/DCAz_posdaughter_selected", "dcaZ of all processed tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/Eta_negdaughter_selected", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    ResoRegistry.add("AnalysisQA/Reso/Phi_negdaughter_selected", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    ResoRegistry.add("AnalysisQA/Reso/Pt_negdaughter_selected", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    ResoRegistry.add("AnalysisQA/Reso/DCAxy_negdaughter_selected", "dcaXY of all processed tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    ResoRegistry.add("AnalysisQA/Reso/DCAz_negdaughter_selected", "dcaZ of all processed tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here

    if (ConfEnableTriggerSelection) {
      for (const std::string& triggerName : softwareTriggers::triggerNames) {
        if (ConfTriggerSwitches->get("Switch", triggerName.c_str())) {
          zorroTriggerNames += triggerName + ",";
        }
      }
      zorroTriggerNames.pop_back();
    }

    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, false, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    colCuts.setCuts(ConfEvtZvtx.value, ConfEvtTriggerCheck.value, ConfEvtTriggerSel.value, ConfEvtOfflineCheck.value, ConfEvtAddOfflineCheck.value, ConfIsRun3.value);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(Track.ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(Track.ConfTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(Track.ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.ConfTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(Track.ConfTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.ConfTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.ConfTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.ConfTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(Track.ConfTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(Track.ConfTrkPIDnSigmaOffsetTPC, Track.ConfTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &TrackRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a
    /// different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0),
    // femtoDreamV0Selection::kDecVtxMax, femtoDreamSelection::kAbsUpperLimit);
    if (ConfIsActivateV0) {
      v0Cuts.setSelection(ConfV0Sign, femtoDreamV0Selection::kV0Sign, femtoDreamSelection::kEqual);
      v0Cuts.setSelection(ConfV0PtMin, femtoDreamV0Selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0PtMax, femtoDreamV0Selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0EtaMax, femtoDreamV0Selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setSelection(ConfV0DCADaughMax, femtoDreamV0Selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0CPAMin, femtoDreamV0Selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0TranRadMin, femtoDreamV0Selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0TranRadMax, femtoDreamV0Selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0DecVtxMax, femtoDreamV0Selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfChildCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfChildCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kPosTrack, ConfChildPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kNegTrack, ConfChildPIDspecies);
      v0Cuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &V0Registry);
      v0Cuts.setInvMassLimits(ConfV0InvMassLowLimit, ConfV0InvMassUpLimit);
      v0Cuts.setIsMother(ConfV0MotherIsLambda);

      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kPosTrack, ConfTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kNegTrack, ConfTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(Track.ConfTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kPosTrack, Track.ConfTrkPIDnSigmaOffsetTPC, Track.ConfTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kNegTrack, Track.ConfTrkPIDnSigmaOffsetTPC, Track.ConfTrkPIDnSigmaOffsetTOF);

      if (ConfV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfV0InvKaonMassLowLimit, ConfV0InvKaonMassUpLimit);
      }
      v0Cuts.setRejectLambda(ConfV0RejectLambdas);
    }

    if (ConfIsActivateReso.value) {
      // resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh,Resonance.ConfDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      // resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh,Resonance.ConfDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCuts.init<aod::femtodreamparticle::ParticleType::kReso,
                    aod::femtodreamparticle::ParticleType::kResoChild>(&qaRegistry, &V0Registry);

      resoCuts.assign(Resonance.ConfResoPTPCThr); // assigns configurable value to class member
      resoCuts.setDaughterPIDSpecies(femtoDreamResoSelection::kPosdaugh, Resonance.ConfDaughterPIDspecies);
      resoCuts.setDaughterPIDSpecies(femtoDreamResoSelection::kNegdaugh, Resonance.ConfDaughterPIDspecies);
      resoCuts.setDaughternSigmaPIDOffset(femtoDreamResoSelection::kPosdaugh, 0.f, 0.f);
      resoCuts.setDaughternSigmaPIDOffset(femtoDreamResoSelection::kNegdaugh, 0.f, 0.f);

      resoCuts.setSelection(Resonance.ConfResoSign, femtoDreamResoSelection::kResoSign, femtoDreamSelection::kEqual);

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
  void initCCDB_Mag_Trig(aod::BCsWithTimestamps::iterator bc)
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

    // Init for zorro to get trigger flags
    if (ConfEnableTriggerSelection) {
      zorro.setCCDBpath(ConfBaseCCDBPathForTriggers);
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

  template <typename CollisionType, typename ParticleType>
  void fillMCParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();
      TrackRegistry.fill(HIST("AnalysisQA/Particle"), pdgCode);
      int particleOrigin = 99;
      int pdgCodeMother = -1;
      // get list of mothers, but it could be empty (for example in case of injected light nuclei)
      auto motherparticlesMC = particleMC.template mothers_as<aod::McParticles>();
      // check pdg code
      TrackRegistry.fill(HIST("AnalysisQA/getGenStatusCode"), particleMC.getGenStatusCode());
      TrackRegistry.fill(HIST("AnalysisQA/getProcess"), particleMC.getProcess());
      // if this fails, the particle is a fake
      if (abs(pdgCode) == abs(ConfTrkPDGCode.value)) {
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
        } else if (particleMC.getProcess() == 4 && particleMC.getGenStatusCode() == -1 && !motherparticlesMC.empty()) {
          // get direct mother
          auto motherparticleMC = motherparticlesMC.front();
          pdgCodeMother = motherparticleMC.pdgCode();
          TrackRegistry.fill(HIST("AnalysisQA/Mother"), pdgCodeMother);
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
          // check if particle is material
          // particle is from inelastic hadronic interaction -> getProcess() == 23
          // particle is generated during transport -> getGenStatusCode() == -1
        } else if (particleMC.getProcess() == 23 && particleMC.getGenStatusCode() == -1) {
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
      if (ConfIsDebug) {
        outputPartsExtMCLabels(outputPartsMC.lastIndex());
        outputDebugPartsMC(pdgCodeMother);
      }
    } else {
      outputPartsMCLabels(-1);
      if (ConfIsDebug) {
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
  template <bool isMC, bool hasItsPid, bool useCentrality, bool analysePbPb, typename V0Type, typename TrackType, typename TrackTypeWithItsPid, typename CollisionType>
  void fillCollisionsAndTracksAndV0(CollisionType const& col, TrackType const& tracks, TrackTypeWithItsPid const& tracksWithItsPid, V0Type const& fullV0s)
  {
    // If triggering is enabled, select only events which were triggered wit our triggers
    if (ConfEnableTriggerSelection) {
      bool zorroSelected = zorro.isSelected(col.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// check if event was selected by triggers of interest
      if (!zorroSelected) {
        return;
      }
    }

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    float mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
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

    if (ConfIsActivateV0.value) {
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

    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracksWithItsPid) {

      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, 0>(track);

      if (track.tpcChi2NCl() < OptionTrackSpecialSelections.ConfTrkMinChi2PerClusterTPC || track.tpcChi2NCl() > OptionTrackSpecialSelections.ConfTrkMaxChi2PerClusterTPC) {
        continue;
      }
      if (track.itsChi2NCl() > OptionTrackSpecialSelections.ConfTrkMaxChi2PerClusterITS) {
        continue;
      }
      if ((OptionTrackSpecialSelections.ConfTrkTPCRefit && !track.hasTPC()) || (OptionTrackSpecialSelections.ConfTrkITSRefit && !track.hasITS())) {
        continue;
      }

      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      TrackRegistry.fill(HIST("AnalysisQA/Chi2ITSTPCperCluster"), track.itsChi2NCl(), track.tpcChi2NCl());
      TrackRegistry.fill(HIST("AnalysisQA/RefitITSTPC"), track.hasITS(), track.hasTPC());

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
      if (ConfIsDebug.value) {
        fillDebugParticle<true, hasItsPid>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }

    if (ConfIsActivateV0.value) {
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
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kPosCuts),
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kPosPID),
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
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegCuts),
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegPID),
                    negtrack.dcaXY(),
                    childIDs,
                    0,
                    0);
        const int rowOfNegTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(col, negtrack, o2::aod::femtodreamparticle::ParticleType::kV0Child);
        }
        std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};
        float FillMass = 0.;
        float FillMassAnti = 0.;
        if (ConfV0MotherIsLambda) {
          FillMass = v0.mLambda();
          FillMassAnti = v0.mAntiLambda();
        } else {
          FillMass = v0.mK0Short();
          FillMassAnti = FillMass;
        }
        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    aod::femtodreamparticle::ParticleType::kV0,
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    FillMass,
                    FillMassAnti);
        if (ConfIsDebug.value) {
          fillDebugParticle<true, false>(postrack); // QA for positive daughter
          fillDebugParticle<true, false>(negtrack); // QA for negative daughter
          fillDebugParticle<false, false>(v0);      // QA for v0
        }
        if constexpr (isMC) {
          fillMCParticle(col, v0, o2::aod::femtodreamparticle::ParticleType::kV0);
        }
      }
    }

    if (ConfIsActivatePhi.value) {
      resoCuts.updateThreshold();

      auto SlicePosdaugh = daughter1.sliceByCached(aod::track::collisionId, col.globalIndex(), cache); // o2::framework defined in AnalysisHelper.h
      auto SliceNegdaugh = daughter2.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);

      for (auto& track1 : SlicePosdaugh) {
        if (!resoCuts.DaughterSelectionPos(track1, Resonance.ConfThreshold.value))
          continue; /// loosest cuts for track1
        for (auto& track2 : SliceNegdaugh) {
          if (!resoCuts.DaughterSelectionNeg(track2, Resonance.ConfThreshold.value))
            continue; /// loosest cuts for track2

          bool ResoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
          float ResoSign = 1.;
          if ((Resonance.ConfDaughterPIDspecies->size() > 1) && (Resonance.ConfDaughterPIDspecies.value[0] != Resonance.ConfDaughterPIDspecies.value[1])) {
            auto [isNormal, WrongCombination] = resoCuts.checkCombination(track1, track2, static_cast<femtoDreamResoSelection::ResoMothers>(Resonance.ConfResoMotherID.value));
            if (WrongCombination) {
              continue;
            }
            if (!isNormal) {
              ResoSign = -1.;
            }
            ResoIsNotAnti = isNormal;
          }
          /// Resos, where both daughters have the same PID are defaulted to sign 1. and ResoIsNotAnti = true

          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter1/Pt"), track1.pt());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter2/Pt"), track2.pt());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter1/Eta"), track1.eta());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter2/Eta"), track2.eta());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter1/DcaXY"), track1.dcaXY());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter2/DcaXY"), track2.dcaXY());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter1/DcaZ"), track1.dcaZ());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter2/DcaZ"), track2.dcaZ());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter1/Phi"), track1.phi());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Daughter2/Phi"), track2.phi());

          /// Get masses for calculating invariant Mass
          /// This only works for the case where the mass of opposite charged particles are the same (for example K+/K- have same mass)
          auto [MassPart1, MassPart2] = resoCuts.getMassDaughters(static_cast<femtoDreamResoSelection::ResoMothers>(Resonance.ConfResoMotherID.value));

          /// Resonance
          ROOT::Math::PtEtaPhiMVector tempD1(track1.pt(), track1.eta(), track1.phi(), MassPart1);
          ROOT::Math::PtEtaPhiMVector tempD2(track2.pt(), track2.eta(), track2.phi(), MassPart2);
          ROOT::Math::PtEtaPhiMVector tempReso = tempD1 + tempD2;
          /// Anti-resonance
          ROOT::Math::PtEtaPhiMVector tempDA1(track1.pt(), track1.eta(), track1.phi(), MassPart2);
          ROOT::Math::PtEtaPhiMVector tempDA2(track2.pt(), track2.eta(), track2.phi(), MassPart1);
          ROOT::Math::PtEtaPhiMVector tempAntiReso = tempDA1 + tempDA2;

          if (ResoIsNotAnti) {
            ResoRegistry.fill(HIST("AnalysisQA/Reso/InvMass"), tempReso.M());
            if (!(tempReso.M() > Resonance.ConfResoInvMassLowLimit.value && tempReso.M() < Resonance.ConfResoInvMassUpLimit.value))
              continue; /// MassCut
            ResoRegistry.fill(HIST("AnalysisQA/Reso/InvMass_phi_selected"), tempReso.M());
          } else {
            ResoRegistry.fill(HIST("AnalysisQA/Reso/InvMassAnti"), tempAntiReso.M());
            if (!(tempAntiReso.M() > Resonance.ConfResoInvMassLowLimit.value && tempAntiReso.M() < Resonance.ConfResoInvMassUpLimit.value))
              continue; /// MassCut
            ResoRegistry.fill(HIST("AnalysisQA/Reso/InvMassAnti_phi_selected"), tempAntiReso.M());
          }

          resoCuts.fillQA<aod::femtodreamparticle::ParticleType::kResoChild,
                          aod::femtodreamparticle::TrackType::kPosChild,
                          aod::femtodreamparticle::TrackType::kNegChild>(track1, track2);

          ResoRegistry.fill(HIST("AnalysisQA/Reso/Pt_posdaughter_selected"), track1.pt());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Pt_negdaughter_selected"), track2.pt());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Eta_posdaughter_selected"), track1.eta());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Eta_negdaughter_selected"), track2.eta());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/DCAxy_posdaughter_selected"), track1.dcaXY());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/DCAxy_negdaughter_selected"), track2.dcaXY());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/DCAz_posdaughter_selected"), track1.dcaZ());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/DCAz_negdaughter_selected"), track2.dcaZ());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Phi_posdaughter_selected"), track1.phi());
          ResoRegistry.fill(HIST("AnalysisQA/Reso/Phi_negdaughter_selected"), track2.phi());

          auto Type = resoCuts.getType(track1, track2); //   kPhiPosdaughTPC_NegdaughTPC
                                                        //   kPhiPosdaughTPC_NegdaughTOF
                                                        //   kPhiPosdaughTPC_NegdaughTPC
                                                        //   kPhiPosdaughTOF_NegdaughTOF as possible output

          auto Bitmask = resoCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track1, track2, ResoSign);

          /// Get Variables for Output
          auto outputReso = tempReso;
          auto outputDaugh1 = tempD1;
          auto outputDaugh2 = tempD2;
          if (!ResoIsNotAnti) {
            outputReso = tempAntiReso;
            outputDaugh1 = tempDA1;
            outputDaugh2 = tempDA2;
          }

          // fill FDParticles
          int postrkId = track1.globalIndex();
          int RowOfPosTrack = -1;
          RowOfPosTrack = getRowDaughters(postrkId, tmpIDtrack);

          childIDs[0] = RowOfPosTrack; // should give me the row
          childIDs[1] = 0;
          outputParts(outputCollision.lastIndex(),
                      track1.pt(),
                      track1.eta(),
                      track1.phi(),
                      aod::femtodreamparticle::ParticleType::kResoChild,
                      Bitmask[1],
                      Bitmask[2],
                      track1.dcaXY(),
                      childIDs,
                      outputDaugh1.M(),
                      outputDaugh2.M()); // fill tempFitVar with dcaXY?
          const int rowPosTrk = outputParts.lastIndex();

          int negtrkId = track2.globalIndex();
          int RowOfNegTrack = -1;
          RowOfNegTrack = getRowDaughters(negtrkId, tmpIDtrack);

          childIDs[0] = 0;
          childIDs[1] = RowOfNegTrack;
          outputParts(outputCollision.lastIndex(),
                      track2.pt(),
                      track2.eta(),
                      track2.phi(),
                      aod::femtodreamparticle::ParticleType::kResoChild,
                      Bitmask[3],
                      Bitmask[4],
                      track2.dcaXY(),
                      childIDs,
                      outputDaugh2.M(),
                      outputDaugh1.M()); // maybe CPA instead of dcaXY()? as tempFitVar?
          const int rowNegTrk = outputParts.lastIndex();

          // Reso
          std::vector<int> IndexChildIds = {rowPosTrk, rowNegTrk};
          outputParts(outputCollision.lastIndex(),
                      outputReso.pt(),
                      outputReso.eta(),
                      outputReso.phi(),
                      Type,
                      0,
                      Bitmask[0],
                      -999.f,
                      IndexChildIds,
                      tempReso.M(),
                      tempAntiReso.M()); // no TempFitVar !!

          // needed?
          if (ConfIsDebug.value) {
            fillDebugParticle<true, false>(track1); // QA for positive daughter
            fillDebugParticle<true, false>(track2); // QA for negative daughter
            fillDebugParticle<false, false, ROOT::Math::PtEtaPhiMVector, true>(outputReso);
          }
        }
      }
    }
  }

  void
    processData(aod::FemtoFullCollision const& col,
                aod::BCsWithTimestamps const&,
                aod::FemtoFullTracks const& tracks,
                o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (ConfUseItsPid.value) {
      fillCollisionsAndTracksAndV0<false, true, true, false>(col, tracks, tracksWithItsPid, fullV0s);
    } else {
      fillCollisionsAndTracksAndV0<false, false, true, false>(col, tracks, tracks, fullV0s);
    }
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processData,
                 "Provide experimental data", true);

  void
    processData_noCentrality(aod::FemtoFullCollision_noCent const& col,
                             aod::BCsWithTimestamps const&,
                             aod::FemtoFullTracks const& tracks,
                             o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (ConfUseItsPid.value) {
      fillCollisionsAndTracksAndV0<false, true, false, false>(col, tracks, tracksWithItsPid, fullV0s);
    } else {
      fillCollisionsAndTracksAndV0<false, false, false, false>(col, tracks, tracks, fullV0s);
    }
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processData_noCentrality,
                 "Provide experimental data without centrality information", false);

  void processData_CentPbPb(aod::FemtoFullCollision_CentPbPb const& col,
                            aod::BCsWithTimestamps const&,
                            aod::FemtoFullTracks const& tracks,
                            o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (ConfUseItsPid.value) {
      fillCollisionsAndTracksAndV0<false, true, true, true>(col, tracks, tracksWithItsPid, fullV0s);
    } else {
      fillCollisionsAndTracksAndV0<false, false, true, true>(col, tracks, tracks, fullV0s);
    }
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processData_CentPbPb,
                 "Provide experimental data with centrality information for PbPb collisions", false);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::FemtoFullMCgenCollisions const&,
                 aod::McParticles const&,
                 soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<false, false, true, false>(col, tracks, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMC, "Provide MC data", false);

  void processMC_noCentrality(aod::FemtoFullCollision_noCent_MC const& col,
                              aod::BCsWithTimestamps const&,
                              soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                              aod::FemtoFullMCgenCollisions const&,
                              aod::McParticles const&,
                              soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<true, false, false, false>(col, tracks, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMC_noCentrality, "Provide MC data without requiring a centrality calibration", false);

  void processMC_CentPbPb(aod::FemtoFullCollisionMC_CentPbPb const& col,
                          aod::BCsWithTimestamps const&,
                          soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                          aod::FemtoFullMCgenCollisions const&,
                          aod::McParticles const&,
                          soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    initCCDB_Mag_Trig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<true, false, true, true>(col, tracks, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMC_CentPbPb, "Provide MC data with centrality information for PbPb collisions", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamProducerTask>(cfgc)};
  return workflow;
}
