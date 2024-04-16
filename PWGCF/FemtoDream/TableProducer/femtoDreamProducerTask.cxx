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

#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamV0Selection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Math/Vector4D.h"
#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
using FemtoFullCollision_noCent_MC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using FemtoHFTrack = FemtoHFTracks::iterator;
} // namespace o2::aod

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

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDHfCand> rowCandCharmHad;
  Produces<aod::FDHfCandMC> rowCandMCCharmHad;
  Produces<aod::FDHfCandMCGen> rowCandCharmHadGen;
  Produces<aod::FDParticlesIndex> outputPartsIndex;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
  Produces<aod::FDExtMCLabels> outputPartsExtMCLabels;

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  /// Event cuts
  FemtoDreamCollisionSelection colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfEvtAddOfflineCheck{"ConfEvtAddOfflineCheck", false, "Evt sel: additional checks for offline selection (not part of sel8 yet)"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", true, "Activate filling of V0 into femtodream tables"};

  /// Lc table
  Configurable<bool> useCent{"useCent", false, "Enable centrality for lc"};
  Configurable<bool> fillLcMCGenTable{"fillLcMCGenTable", false, "Enable mc gen table for lc"};
  Configurable<bool> keepOnlySignalMc{"keepOnlySignalMc", false, "Fill MC tree only with signal candidates"};
  Configurable<bool> keepOnlyBkg{"keepOnlyBkg", false, "Fill MC tree only with background candidates"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};

  Configurable<bool> ConfTrkRejectNotPropagated{"ConfTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> ConfRejectITSHitandTOFMissing{ "ConfRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};
  Configurable<int> ConfTrkPDGCode{"ConfTrkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  FemtoDreamTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

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
  Configurable<float> ConfV0InvKaonMassLowLimit{"ConfV0InvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> ConfV0InvKaonMassUpLimit{"ConfV0InvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

  Configurable<std::vector<float>> ConfChildCharge{"ConfChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
  Configurable<std::vector<float>> ConfChildEtaMax{"ConfChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
  Configurable<std::vector<float>> ConfChildTPCnClsMin{"ConfChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
  Configurable<std::vector<float>> ConfChildDCAMin{"ConfChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> ConfChildPIDnSigmaMax{"ConfChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
  Configurable<std::vector<int>> ConfChildPIDspecies{"ConfChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why
  using CandidateLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;
  using CandidateLcMC = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>;
  using GenratedMC = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry TrackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry V0Registry{"V0", {}, OutputObjHandlingPolicy::AnalysisObject};

  HfHelper hfHelper;
  int mRunNumber;
  float mMagField;
  //  std::array<float, 3> outputMlPKPi{-1., -1., -1.};
  //  std::array<float, 3> outputMlPiKP{-1., -1., -1.};
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    std::array<bool, 7> processes = {doprocessData, doprocessMC,
                                     doprocessMCNoCentrality, doprocessDataCharmHad, doprocessMCCharmHad, doprocessDataCharmHadWithML, doprocessMCCharmHadWithML};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    // if (doprocessData == false && doprocessMC == false && doprocessMCNoCentrality == false) {
    //   LOGF(fatal, "Neither processData nor processMC enabled. Please choose one.");
    // }
    // if ((doprocessData == true && doprocessMC == true) || (doprocessData == true && doprocessMCNoCentrality == true) || (doprocessMC == true && doprocessMCNoCentrality == true)) {
    //   LOGF(fatal,
    //        "Cannot enable more than one process switch at the same time. "
    //        "Please choose one.");
    // }

    int CutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    TrackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});
    V0Registry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});

    colCuts.setCuts(ConfEvtZvtx.value, ConfEvtTriggerCheck.value, ConfEvtTriggerSel.value, ConfEvtOfflineCheck.value, ConfEvtAddOfflineCheck.value, ConfIsRun3.value);
    colCuts.init(&qaRegistry);

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
    trackCuts.setPIDSpecies(ConfTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
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

      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kPosTrack, ConfTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kNegTrack, ConfTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(ConfTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kPosTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kNegTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);

      if (ConfV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfV0InvKaonMassLowLimit, ConfV0InvKaonMassUpLimit);
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

  template <bool isTrackOrV0, bool isHF = false, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
      if constexpr (isHF) {
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
                         particle.tpcNSigmaPi(),
                         particle.tpcNSigmaKa(),
                         particle.tpcNSigmaPr(),
                         particle.tofNSigmaPi(),
                         particle.tofNSigmaKa(),
                         particle.tofNSigmaPr(),
                         -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.);

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
                         particle.tofNSigmaEl(),
                         particle.tofNSigmaPi(),
                         particle.tofNSigmaKa(),
                         particle.tofNSigmaPr(),
                         particle.tofNSigmaDe(),
                         -999., -999., -999., -999., -999., -999.);
      }
    } else {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       particle.dcaV0daughters(),
                       particle.v0radius(),
                       particle.x(),
                       particle.y(),
                       particle.z(),
                       particle.mK0Short());
    }
  }

  template <typename CollisionType, typename ParticleType>
  void fillMCParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();
      int particleOrigin = 99;
      auto motherparticleMC = particleMC.template mothers_as<aod::McParticles>().front();

      if (abs(pdgCode) == abs(ConfTrkPDGCode.value)) {
        if ((col.has_mcCollision() && (particleMC.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
        } else if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (motherparticleMC.isPhysicalPrimary() && particleMC.getProcess() == 4) {
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
        } else if (particleMC.getGenStatusCode() == -1) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
      } else {
        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
      if (ConfIsDebug) {
        outputPartsExtMCLabels(outputPartsMC.lastIndex());
        outputDebugPartsMC(motherparticleMC.pdgCode());
      }
    } else {
      outputPartsMCLabels(-1);
      if (ConfIsDebug) {
        outputPartsExtMCLabels(-1);
      }
    }
  }

  template <bool isMC, bool useCentrality, typename V0Type, typename TrackType, typename CollisionType>
  void fillCollisionsAndTracksAndV0(CollisionType const& col, TrackType const& tracks, V0Type const& fullV0s)
  {
    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    float mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      if constexpr (useCentrality) {
        mult = col.centFT0M();
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

    outputCollision(vtxZ, mult, multNtr, spher, mMagField);

    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

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
        fillDebugParticle<true>(track);
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
        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    aod::femtodreamparticle::ParticleType::kV0,
                    cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    v0.mLambda(),
                    v0.mAntiLambda());
        if (ConfIsDebug) {
          fillDebugParticle<true>(postrack); // QA for positive daughter
          fillDebugParticle<true>(negtrack); // QA for negative daughter
          fillDebugParticle<false>(v0);      // QA for v0
        }
        if constexpr (isMC) {
          fillMCParticle(col, v0, o2::aod::femtodreamparticle::ParticleType::kV0);
        }
      }
    }
  }

  template <bool isMC, typename TrackType>
  bool fillTracksForCharmHadron(TrackType const& tracks, o2::aod::FemtoHFTrack const& prong0, o2::aod::FemtoHFTrack const& prong1, o2::aod::FemtoHFTrack const& prong2, int candSize)
  {

    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    // std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    bool fIsTrackFilled = false;

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      if ((candSize == 1) && (track.globalIndex() == prong0.globalIndex() || track.globalIndex() == prong1.globalIndex() || track.globalIndex() == prong2.globalIndex()))
        continue;

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, true>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

      // track global index
      outputPartsIndex(track.globalIndex());
      // now the table is filled

      outputParts(outputCollision.lastIndex() + 1,
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtodreamparticle::ParticleType::kTrack,
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);
      fIsTrackFilled = true;
      // tmpIDtrack.push_back(track.globalIndex());
      if (ConfIsDebug.value) {
        fillDebugParticle<true, true>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }
    return fIsTrackFilled;
  }

  template <bool isMC, bool useCharmMl, typename TrackType, typename CollisionType, typename CandType>
  void fillCharmHadronTable(CollisionType const& col, TrackType const& tracks, CandType const& candidates)
  {
    const auto vtxZ = col.posZ();
    const auto sizeCand = candidates.size();
    if (sizeCand == 0)
      return;

    const auto spher = colCuts.computeSphericity(col, tracks);
    float mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      if (useCent) {
        mult = col.centFT0M();
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

    if (colCuts.isEmptyCollision(col, tracks, trackCuts)) {
      return;
    }
    // Filling candidate properties
    rowCandCharmHad.reserve(sizeCand);
    bool isTrackFilled = false;
    for (const auto& candidate : candidates) {
      std::array<float, 3> outputMlPKPi{-1., -1., -1.};
      std::array<float, 3> outputMlPiKP{-1., -1., -1.};
      if constexpr (useCharmMl) {
        /// fill with ML information
        /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
        if (candidate.mlProbLcToPKPi().size() > 0) {
          outputMlPKPi.at(0) = candidate.mlProbLcToPKPi()[0]; /// bkg score
          outputMlPKPi.at(1) = candidate.mlProbLcToPKPi()[1]; /// prompt score
          outputMlPKPi.at(2) = candidate.mlProbLcToPKPi()[2]; /// non-prompt score
        }
        if (candidate.mlProbLcToPiKP().size() > 0) {
          outputMlPiKP.at(0) = candidate.mlProbLcToPiKP()[0]; /// bkg score
          outputMlPiKP.at(1) = candidate.mlProbLcToPiKP()[1]; /// prompt score
          outputMlPiKP.at(2) = candidate.mlProbLcToPiKP()[2]; /// non-prompt score
        }
      }
      auto trackPos1 = candidate.template prong0_as<o2::aod::FemtoHFTracks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<o2::aod::FemtoHFTracks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<o2::aod::FemtoHFTracks>(); // positive daughter (negative for the antiparticles)
      bool isMcCandidateSignal = false;

      if constexpr (isMC) {
        isMcCandidateSignal = std::abs(candidate.flagMcMatchRec()) == (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi);
      }

      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float BDTScoreBkg,
                           float BDTScorePrompt,
                           float BDTScoreFD,
                           float FunctionInvMass,
                           float FunctionY) {
        if (FunctionSelection >= 1 && (!isMC || (/*keep all*/ (!keepOnlySignalMc && !keepOnlyBkg) || /*keep only signal*/ (keepOnlySignalMc && isMcCandidateSignal) || /*keep only background and downsample it*/ (keepOnlyBkg && !isMcCandidateSignal)))){
        // Fill tracks if it is not filled for Lc Candidate in an event
            if (!isTrackFilled) {
                isTrackFilled = fillTracksForCharmHadron<false>(tracks, trackPos1, trackNeg, trackPos2, sizeCand);

                // If track filling was successful, fill the collision table
                if (isTrackFilled) {
                    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
                }
            }

            // fill collision table if track table is filled, i.e., there is at least one Lc-p pair
            if (isTrackFilled) {
                // Row for candidate charm hadron
                rowCandCharmHad(
                    outputCollision.lastIndex(),
                    trackPos1.sign() + trackNeg.sign() + trackPos2.sign(),
                    trackPos1.globalIndex(),
                    trackNeg.globalIndex(),
                    trackPos2.globalIndex(),
                    trackPos1.pt(),
                    trackNeg.pt(),
                    trackPos2.pt(),
                    trackPos1.eta(),
                    trackNeg.eta(),
                    trackPos2.eta(),
                    trackPos1.phi(),
                    trackNeg.phi(),
                    trackPos2.phi(),
                    1 << CandFlag,
                    BDTScoreBkg,
                    BDTScorePrompt,
                    BDTScoreFD,
                    FunctionInvMass,
                    candidate.pt(),
                    candidate.p(),
                    candidate.eta(),
                    candidate.phi(),
                    FunctionY);

                // Row for MC candidate charm hadron (if constexpr isMC)
                if constexpr (isMC) {
                    rowCandMCCharmHad(
                        candidate.flagMcMatchRec(),
                        candidate.originMcRec());
                }
            }
      } };

      fillTable(0, candidate.isSelLcToPKPi(), outputMlPKPi.at(0), outputMlPKPi.at(1), outputMlPKPi.at(2), hfHelper.invMassLcToPKPi(candidate), hfHelper.yLc(candidate));
      fillTable(1, candidate.isSelLcToPiKP(), outputMlPiKP.at(0), outputMlPiKP.at(1), outputMlPiKP.at(2), hfHelper.invMassLcToPiKP(candidate), hfHelper.yLc(candidate));
    }
  }

  template <typename CollisionType, typename TrackType, typename ParticleType>
  void fillCharmHadMCGen(CollisionType const& col, TrackType const& tracks, ParticleType particles)
  {

    if (!colCuts.isSelectedCollision(col)) {
      return;
    }

    if (colCuts.isEmptyCollision(col, tracks, trackCuts)) {
      return;
    }
    // Filling particle properties
    rowCandCharmHadGen.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        rowCandCharmHadGen(
          outputCollision.lastIndex(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassLambdaCPlus),
          particle.flagMcMatchGen(),
          particle.originMcGen());
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
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<false, true>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processData,
                 "Provide experimental data", false);
  void
    processDataCharmHad(aod::FemtoFullCollision const& col,
                        aod::FemtoHFTracks const& tracks,
                        soa::Filtered<CandidateLc> const& candidates)
  {
    fillCharmHadronTable<false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processDataCharmHad,
                 "Provide experimental data for charm hadron femto", false);

  void
    processDataCharmHadWithML(aod::FemtoFullCollision const& col,
                              aod::FemtoHFTracks const& tracks,
                              soa::Filtered<soa::Join<CandidateLc, aod::HfMlLcToPKPi>> const& candidates)
  {

    fillCharmHadronTable<false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processDataCharmHadWithML,
                 "Provide experimental data for charm hadron femto with ml", true);

  void processMCCharmHad(aod::FemtoFullCollisionMC const& col,
                         soa::Join<aod::FemtoHFTracks, aod::McTrackLabels> const& tracks,
                         soa::Filtered<CandidateLcMC> const& candidates,
                         GenratedMC const& particles)
  {
    fillCharmHadronTable<true, false>(col, tracks, candidates);
    if (fillLcMCGenTable)
      fillCharmHadMCGen(col, tracks, particles);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMCCharmHad, "Provide MC for charm hadron", false);

  void processMCCharmHadWithML(aod::FemtoFullCollisionMC const& col,
                               soa::Join<aod::FemtoHFTracks, aod::McTrackLabels> const& tracks,
                               soa::Filtered<soa::Join<CandidateLcMC, aod::HfMlLcToPKPi>> const& candidates,
                               GenratedMC const& particles)
  {
    fillCharmHadronTable<true, true>(col, tracks, candidates);
    if (fillLcMCGenTable)
      fillCharmHadMCGen(col, tracks, particles);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMCCharmHadWithML, "Provide MC for charm hadron with ml", false);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::McCollisions const& mcCollisions,
                 aod::McParticles const& mcParticles,
                 soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<true, true>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMC, "Provide MC data", false);

  void processMCNoCentrality(aod::FemtoFullCollision_noCent_MC const& col,
                             aod::BCsWithTimestamps const&,
                             soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                             aod::McCollisions const& mcCollisions,
                             aod::McParticles const& mcParticles,
                             soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<true, false>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMCNoCentrality, "Provide MC data without requiring a centrality calibration", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamProducerTask>(cfgc)};
  return workflow;
}
