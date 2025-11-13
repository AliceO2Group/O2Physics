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

/// \file femtoWorldProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#include "PWGCF/FemtoWorld/Core/FemtoWorldCollisionSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPairCleaner.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPhiSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldTrackSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldV0Selection.h"
#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <TDatabasePDG.h> // FIXME

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollisionRun2 = soa::Join<aod::Collisions,
                                         aod::EvSels,
                                         aod::Mults,
                                         aod::CentRun2V0Ms>::iterator;
using FemtoFullCollisionRun3 = soa::Join<aod::Collisions,
                                         aod::EvSels,
                                         aod::Mults,
                                         aod::CentFT0Ms>::iterator;
using FemtoFullTracks = soa::Join<aod::FullTracks,
                                  aod::TracksDCA, aod::TOFSignal,
                                  aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                  aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                  aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                  aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe, aod::pidTOFbeta>;
using FemtoPhiTracks = soa::Join<aod::FullTracks,
                                 aod::TracksDCA, aod::TOFSignal,
                                 aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                 aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                 aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                 aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe, aod::pidTOFbeta>;
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

struct femtoWorldProducerTask {

  Produces<aod::FemtoWorldCollisions> outputCollision;
  Produces<aod::FemtoWorldParticles> outputParts;

  /// Event cuts
  FemtoWorldCollisionSelection colCuts;
  FemtoWorldTrackSelection trackCuts;
  FemtoWorldV0Selection v0Cuts;
  TrackSelection* o2PhysicsTrackSelection;

  struct : ConfigurableGroup {
    // Configurables for D0/D0bar mesons
    Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
    Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
    Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};

    // Choose if filtering or skimming version is run
    Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};

    // Choose if running on converted data or pilot beam
    Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"};
    Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

    Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
    Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
    Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

    Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0 table"};
    Configurable<bool> ConfStorePhi{"ConfStorePhi", true, "True: store Phi table"};

    // just sanity check to make sure in case there are problems in conversion or MC production it does not affect results
    Configurable<bool> ConfRejectNotPropagatedTracks{"ConfRejectNotPropagatedTracks", false, "True: reject not propagated tracks"};
    Configurable<bool> ConfRejectITSHitandTOFMissing{"ConfRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};

    Configurable<std::vector<float>> ConfTrkCharge{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kSign, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkPtmin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.4f, 0.6f, 0.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kpTMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkEta{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kEtaMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCfClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCcRowsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kITSnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kITSnClsIbMin, "Track selection: ")};
    Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kDCAxyMax, "Track selection: ")}; /// here we need an open cut to do the DCA fits later on!
    Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 3.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kDCAzMax, "Track selection: ")};
    /// \todo Reintegrate PID to the general selection container
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kPIDnSigmaMax, "Track selection: ")};
    Configurable<std::vector<int>> ConfTrkTPIDspecies{"ConfTrkTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

    /// \todo Labeled array (see Track-Track task)
    /// V0
    Configurable<std::vector<float>> ConfV0Sign{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0PtMin{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kpTV0Min, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kpTV0Min, "V0 selection: ")};
    Configurable<std::vector<float>> ConfDCAV0DaughMax{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kDCAV0DaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kDCAV0DaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> ConfCPAV0Min{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kCPAV0Min, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kCPAV0Min, "V0 selection: ")};

    Configurable<std::vector<float>> V0TranRadV0Min{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kTranRadV0Min, "ConfV0"), std::vector<float>{0.2f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kTranRadV0Min, "V0 selection: ")};
    Configurable<std::vector<float>> V0TranRadV0Max{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kTranRadV0Max, "ConfV0"), std::vector<float>{100.f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kTranRadV0Max, "V0 selection: ")};
    Configurable<std::vector<float>> V0DecVtxMax{FemtoWorldV0Selection::getSelectionName(femtoWorldV0Selection::kDecVtxMax, "ConfV0"), std::vector<float>{100.f}, FemtoWorldV0Selection::getSelectionHelper(femtoWorldV0Selection::kDecVtxMax, "V0 selection: ")};

    Configurable<std::vector<float>> ConfV0DaughCharge{"ConfV0DaughCharge", std::vector<float>{-1, 1}, "V0 Daugh sel: Charge"};
    Configurable<std::vector<float>> ConfDaughEta{"ConfDaughEta", std::vector<float>{0.8f}, "V0 Daugh sel: max eta"};
    Configurable<std::vector<float>> ConfV0DaughTPCnclsMin{"ConfV0DaughTPCnclsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Daugh sel: Min. nCls TPC"};
    Configurable<std::vector<float>> ConfV0DaughDCAMin{"ConfV0DaughDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> ConfV0DaughPIDnSigmaMax{"ConfV0DaughPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Daugh sel: Max. PID nSigma TPC"};

    Configurable<std::vector<int>> ConfV0DaughTPIDspecies{"ConfV0DaughTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Daugh sel: Particles species for PID"};

    Configurable<float> ConfInvMassLowLimit{"ConfInvMassLowLimit", 1.005, "Lower limit of the V0 invariant mass"};
    Configurable<float> ConfInvMassUpLimit{"ConfInvMassUpLimit", 1.035, "Upper limit of the V0 invariant mass"};

    Configurable<bool> ConfRejectKaons{"ConfRejectKaons", false, "Switch to reject kaons"};
    Configurable<float> ConfInvKaonMassLowLimit{"ConfInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> ConfInvKaonMassUpLimit{"ConfInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    // PHI Daughters (Kaons)
    Configurable<float> ConfInvMassLowLimitPhi{"ConfInvMassLowLimitPhi", 1.011, "Lower limit of the Phi invariant mass"}; // change that to do invariant mass cut
    Configurable<float> ConfInvMassUpLimitPhi{"ConfInvMassUpLimitPhi", 1.027, "Upper limit of the Phi invariant mass"};
    Configurable<bool> ConfRejectKaonsPhi{"ConfRejectKaonsPhi", false, "Switch to reject kaons"};
    Configurable<float> ConfInvKaonMassLowLimitPhi{"ConfInvKaonMassLowLimitPhi", 0.48, "Lower limit of the Phi invariant mass for Kaon rejection"};
    Configurable<float> ConfInvKaonMassUpLimitPhi{"ConfInvKaonMassUpLimitPhi", 0.515, "Upper limit of the Phi invariant mass for Kaon rejection"};
    Configurable<bool> ConfNsigmaTPCTOFKaon{"ConfNsigmaTPCTOFKaon", true, "Use TPC and TOF for PID of Kaons"};
    Configurable<float> ConfNsigmaCombinedKaon{"ConfNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined) for momentum > 0.4"};
    Configurable<float> ConfNsigmaTPCKaon{"ConfNsigmaTPCKaon", 3.0, "TPC Kaon Sigma for momentum < 0.4"};
    // For Phi (daughter 1)
    Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
    Configurable<float> cfgPtLowPart1{"cfgPtLowPart1", 0.14, "Lower limit for Pt for the first particle"};
    Configurable<float> cfgPtHighPart1{"cfgPtHighPart1", 4.0, "Higher limit for Pt for the first particle"};
    Configurable<float> cfgPLowPart1{"cfgPLowPart1", 0.14, "Lower limit for P for the first particle"};
    Configurable<float> cfgPHighPart1{"cfgPHighPart1", 1.5, "Higher limit for P for the first particle"};
    Configurable<float> cfgEtaLowPart1{"cfgEtaLowPart1", -0.8, "Lower limit for Eta for the first particle"};
    Configurable<float> cfgEtaHighPart1{"cfgEtaHighPart1", 0.8, "Higher limit for Eta for the first particle"};
    Configurable<float> cfgDcaXYPart1{"cfgDcaXYPart1", 2.4, "Value for DCA_XY for the first particle"};
    Configurable<float> cfgDcaZPart1{"cfgDcaZPart1", 3.2, "Value for DCA_Z for the first particle"};
    Configurable<int> cfgTpcClPart1{"cfgTpcClPart1", 88, "Number of tpc clasters for the first particle"};             // min number of found TPC clusters
    Configurable<int> cfgTpcCrosRoPart1{"cfgTpcCrosRoPart1", 70, "Number of tpc crossed rows for the first particle"}; // min number of crossed rows
    Configurable<float> cfgChi2TpcPart1{"cfgChi2TpcPart1", 4.0, "Chi2 / cluster for the TPC track segment for the first particle"};
    Configurable<float> cfgChi2ItsPart1{"cfgChi2ItsPart1", 36.0, "Chi2 / cluster for the ITS track segment for the first particle"};

    // For Phi (daughter 2)
    Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};
    Configurable<float> cfgPtLowPart2{"cfgPtLowPart2", 0.14, "Lower limit for Pt for the second particle"};
    Configurable<float> cfgPtHighPart2{"cfgPtHighPart2", 4.0, "Higher limit for Pt for the second particle"};
    Configurable<float> cfgPLowPart2{"cfgPLowPart2", 0.14, "Lower limit for P for the second particle"};
    Configurable<float> cfgPHighPart2{"cfgPHighPart2", 1.5, "Higher limit for P for the second particle"};
    Configurable<float> cfgEtaLowPart2{"cfgEtaLowPart2", -0.8, "Lower limit for Eta for the second particle"};
    Configurable<float> cfgEtaHighPart2{"cfgEtaHighPart2", 0.8, "Higher limit for Eta for the second particle"};
    Configurable<float> cfgDcaXYPart2{"cfgDcaXYPart2", 2.4, "Value for DCA_XY for the second particle"};
    Configurable<float> cfgDcaZPart2{"cfgDcaZPart2", 3.2, "Value for DCA_Z for the second particle"};
    Configurable<int> cfgTpcClPart2{"cfgTpcClPart2", 88, "Number of tpc clasters for the second particle"};             // min number of found TPC clusters
    Configurable<int> cfgTpcCrosRoPart2{"cfgTpcCrosRoPart2", 70, "Number of tpc crossed rows for the second particle"}; // min number of crossed rows
    Configurable<float> cfgChi2TpcPart2{"cfgChi2TpcPart2", 4.0, "Chi2 / cluster for the TPC track segment for the second particle"};
    Configurable<float> cfgChi2ItsPart2{"cfgChi2ItsPart2", 36.0, "Chi2 / cluster for the ITS track segment for the second particle"};
  } config;

  // PHI Candidates
  FemtoWorldPhiSelection PhiCuts;
  // Configurable<std::vector<float>> ConfPhiSign{FemtoWorldPhiSelection::getSelectionName(femtoWorldPhiSelection::kPhiSign, "ConfPhi"), std::vector<float>{-1, 1}, FemtoWorldPhiSelection::getSelectionHelper(femtoWorldPhiSelection::kPhiSign, "Phi selection: ")};
  // Configurable<std::vector<float>> ConfPhiPtMin{FemtoWorldPhiSelection::getSelectionName(femtoWorldPhiSelection::kpTPhiMin, "ConfPhi"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoWorldPhiSelection::getSelectionHelper(femtoWorldPhiSelection::kpTPhiMin, "Phi selection: ")};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  HfHelper hfHelper;

  void init(InitContext&)
  {
    colCuts.setCuts(config.ConfEvtZvtx, config.ConfEvtTriggerCheck, config.ConfEvtTriggerSel, config.ConfEvtOfflineCheck, config.ConfIsRun3);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(config.ConfTrkCharge, femtoWorldTrackSelection::kSign, femtoWorldSelection::kEqual);
    trackCuts.setSelection(config.ConfTrkPtmin, femtoWorldTrackSelection::kpTMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkEta, femtoWorldTrackSelection::kEtaMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(config.ConfTrkTPCnclsMin, femtoWorldTrackSelection::kTPCnClsMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkTPCfCls, femtoWorldTrackSelection::kTPCfClsMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkTPCcRowsMin, femtoWorldTrackSelection::kTPCcRowsMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkTPCsCls, femtoWorldTrackSelection::kTPCsClsMax, femtoWorldSelection::kUpperLimit);
    trackCuts.setSelection(config.ConfTrkITSnclsMin, femtoWorldTrackSelection::kITSnClsMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkITSnclsIbMin, femtoWorldTrackSelection::kITSnClsIbMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(config.ConfTrkDCAxyMax, femtoWorldTrackSelection::kDCAxyMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(config.ConfTrkDCAzMax, femtoWorldTrackSelection::kDCAzMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(config.ConfTrkPIDnSigmaMax, femtoWorldTrackSelection::kPIDnSigmaMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(config.ConfTrkTPIDspecies);
    trackCuts.init<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild, aod::femtoworldparticle::cutContainerType>(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0), femtoWorldV0Selection::kDecVtxMax, femtoWorldSelection::kAbsUpperLimit);
    if (config.ConfStoreV0) {
      v0Cuts.setSelection(config.ConfV0Sign, femtoWorldV0Selection::kV0Sign, femtoWorldSelection::kEqual);
      v0Cuts.setSelection(config.ConfV0PtMin, femtoWorldV0Selection::kpTV0Min, femtoWorldSelection::kLowerLimit);
      v0Cuts.setSelection(config.ConfDCAV0DaughMax, femtoWorldV0Selection::kDCAV0DaughMax, femtoWorldSelection::kUpperLimit);
      v0Cuts.setSelection(config.ConfCPAV0Min, femtoWorldV0Selection::kCPAV0Min, femtoWorldSelection::kLowerLimit);

      v0Cuts.setChildCuts(femtoWorldV0Selection::kPosTrack, config.ConfV0DaughCharge, femtoWorldTrackSelection::kSign, femtoWorldSelection::kEqual);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kPosTrack, config.ConfDaughEta, femtoWorldTrackSelection::kEtaMax, femtoWorldSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kPosTrack, config.ConfV0DaughTPCnclsMin, femtoWorldTrackSelection::kTPCnClsMin, femtoWorldSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kPosTrack, config.ConfV0DaughDCAMin, femtoWorldTrackSelection::kDCAMin, femtoWorldSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kPosTrack, config.ConfV0DaughPIDnSigmaMax, femtoWorldTrackSelection::kPIDnSigmaMax, femtoWorldSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kNegTrack, config.ConfV0DaughCharge, femtoWorldTrackSelection::kSign, femtoWorldSelection::kEqual);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kNegTrack, config.ConfDaughEta, femtoWorldTrackSelection::kEtaMax, femtoWorldSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kNegTrack, config.ConfV0DaughTPCnclsMin, femtoWorldTrackSelection::kTPCnClsMin, femtoWorldSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kNegTrack, config.ConfV0DaughDCAMin, femtoWorldTrackSelection::kDCAMin, femtoWorldSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoWorldV0Selection::kNegTrack, config.ConfV0DaughPIDnSigmaMax, femtoWorldTrackSelection::kPIDnSigmaMax, femtoWorldSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoWorldV0Selection::kPosTrack, config.ConfV0DaughTPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoWorldV0Selection::kNegTrack, config.ConfV0DaughTPIDspecies);
      v0Cuts.init<aod::femtoworldparticle::ParticleType::kV0, aod::femtoworldparticle::ParticleType::kV0Child, aod::femtoworldparticle::cutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(config.ConfInvMassLowLimit, config.ConfInvMassUpLimit);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoWorldV0Selection::kPosTrack, config.ConfRejectNotPropagatedTracks);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoWorldV0Selection::kNegTrack, config.ConfRejectNotPropagatedTracks);

      if (config.ConfRejectKaons) {
        v0Cuts.setKaonInvMassLimits(config.ConfInvKaonMassLowLimit, config.ConfInvKaonMassUpLimit);
      }
      if (config.ConfRejectITSHitandTOFMissing) {
        o2PhysicsTrackSelection = new TrackSelection(getGlobalTrackSelection());
        o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      }
    }

    if (config.ConfStorePhi) {
      PhiCuts.init<aod::femtoworldparticle::ParticleType::kPhi, aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::cutContainerType>(&qaRegistry);
      if (config.ConfRejectKaonsPhi) {
        //! PhiCuts.setKaonInvMassLimits(ConfInvKaonMassLowLimitPhi, ConfInvKaonMassUpLimitPhi);
      }
      if (config.ConfRejectITSHitandTOFMissing) {
        o2PhysicsTrackSelection = new TrackSelection(getGlobalTrackSelection());
        o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      }
    }

    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // changed long to float because of the MegaLinter
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  // PID
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    //|nsigma_TPC| < 5 for p < 0.4 GeV/c
    //|nsigma_combined| < 5 for p > 0.4

    // using configurables:
    // ConfNsigmaTPCTOFKaon -> are we doing TPC TOF PID for Kaons? (boolean)
    // ConfNsigmaTPCKaon -> TPC Kaon Sigma for momentum < 0.4
    // ConfNsigmaCombinedKaon -> TPC and TOF Kaon Sigma (combined) for momentum > 0.4
    if (config.ConfNsigmaTPCTOFKaon) {
      if (mom < 0.4) {
        if (TMath::Abs(nsigmaTPCK) < config.ConfNsigmaTPCKaon) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.4) {
        if (TMath::Hypot(nsigmaTOFK, nsigmaTPCK) < config.ConfNsigmaCombinedKaon) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // get magnetic field for run
    if (mRunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();
    float output = -999;

    if (config.ConfIsRun3 && !config.ConfIsForceGRP) {
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
  void processProdTrackRun2(aod::FemtoFullCollisionRun2 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    /*auto bc = col.bc_as<aod::BCsWithTimestamps>();
    if (mRunNumber != bc.runNumber()) {
      mMagField = getMagneticFieldTesla(bc.timestamp());
      mRunNumber = bc.runNumber();
    }*/
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centRun2V0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centRun2V0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, mult, spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackRun2, "Produce Femto tables", true);

  void processProdTrackV0Run2(aod::FemtoFullCollisionRun2 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks,
                              o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centRun2V0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centRun2V0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, mult, spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    for (auto& v0 : fullV0s) {
      auto postrack = v0.posTrack_as<aod::FemtoFullTracks>();
      auto negtrack = v0.negTrack_as<aod::FemtoFullTracks>(); ///\tocheck funnily enough if we apply the filter the sign of Pos and Neg track is always negative
      // const auto dcaXYpos = postrack.dcaXY();
      // const auto dcaZpos = postrack.dcaZ();
      // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
      v0Cuts.fillLambdaQA(col, v0, postrack, negtrack);

      if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
        continue;
      }

      if (config.ConfRejectITSHitandTOFMissing) {
        // Uncomment only when TOF timing is solved
        // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack, TrackSelection::TrackCuts::kITSHits);
        // bool itsHit = o2PhysicsTrackSelection->IsSelected(negtrack, TrackSelection::TrackCuts::kITSHits);
      }

      v0Cuts.fillQA<aod::femtoworldparticle::ParticleType::kV0, aod::femtoworldparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
      auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtoworldparticle::cutContainerType>(col, v0, postrack, negtrack);

      if ((cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0) > 0) && (cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts) > 0) && (cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts) > 0)) {
        int postrackID = v0.posTrackId();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    v0.positivept(),
                    v0.positiveeta(),
                    v0.positivephi(),
                    0,     // v0.p(),
                    0,     // mass
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kV0Child,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    postrack.sign(),
                    postrack.beta(),
                    postrack.itsChi2NCl(),
                    postrack.tpcChi2NCl(),
                    postrack.tpcNSigmaKa(),
                    postrack.tofNSigmaKa(),
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
                    -999.,
                    -999.,
                    -999.);
        const int rowOfPosTrack = outputParts.lastIndex();
        int negtrackID = v0.negTrackId();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;
        outputParts(outputCollision.lastIndex(),
                    v0.negativept(),
                    v0.negativeeta(),
                    v0.negativephi(),
                    0,     // momentum
                    0,     // mass
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kV0Child,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    negtrack.sign(),
                    negtrack.beta(),
                    negtrack.itsChi2NCl(),
                    negtrack.tpcChi2NCl(),
                    negtrack.tpcNSigmaKa(),
                    negtrack.tofNSigmaKa(),
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
                    -999.,
                    -999.,
                    -999.);
        const int rowOfNegTrack = outputParts.lastIndex();
        int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    0,     // momentum
                    0,     // mass
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kV0,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    v0.mLambda(),
                    v0.mAntiLambda(),
                    postrack.sign(),
                    postrack.beta(),
                    postrack.itsChi2NCl(),
                    postrack.tpcChi2NCl(),
                    postrack.tpcNSigmaKa(),
                    postrack.tofNSigmaKa(),
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
                    -999.,
                    -999.,
                    -999.);
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackV0Run2, "Produce Femto tables", true);

  void processProdTrackPhiRun2(aod::FemtoFullCollisionRun2 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, o2::aod::V0Datas const&)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centRun2V0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centRun2V0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, mult, spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      if ((p1.trackType() == o2::aod::track::TrackTypeEnum::Run2Tracklet) || (p2.trackType() == o2::aod::track::TrackTypeEnum::Run2Tracklet)) {
        continue;
      } else if (p1.globalIndex() == p2.globalIndex()) { // checking not to correlate same particles
        continue;
      } else if ((p1.pt() < config.cfgPtLowPart1) || (p1.pt() > config.cfgPtHighPart1)) { // pT cuts for part1
        continue;
      } else if ((p1.p() < config.cfgPLowPart1) || (p1.p() > config.cfgPHighPart1)) { // p cuts for part1
        continue;
      } else if ((p1.eta() < config.cfgEtaLowPart1) || (p1.eta() > config.cfgEtaHighPart1)) { // eta cuts for part1
        continue;
      } else if ((p2.pt() < config.cfgPtLowPart2) || (p2.pt() > config.cfgPtHighPart2)) { // pT cuts for part2
        continue;
      } else if ((p2.p() < config.cfgPLowPart2) || (p2.p() > config.cfgPHighPart2)) { // p cuts for part2
        continue;
      } else if ((p2.eta() < config.cfgEtaLowPart2) || (p2.eta() > config.cfgEtaHighPart2)) { // eta for part2
        continue;
      } else if (!(IsKaonNSigma(p1.p(), p1.tpcNSigmaKa(), p1.tofNSigmaKa()))) { // PID for Kaons
        continue;
      } else if (!(IsKaonNSigma(p2.p(), p2.tpcNSigmaKa(), p2.tofNSigmaKa()))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;
      float mMassOne = TDatabasePDG::Instance()->GetParticle(config.ConfPDGCodePartOne)->Mass(); // FIXME: Get from the PDG service of the common header
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(config.ConfPDGCodePartTwo)->Mass(); // FIXME: Get from the PDG service of the common header

      part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      float phiPt = sumVec.Pt();
      float phiP = sumVec.P();
      float phiM = sumVec.M();

      if (((phiM < config.ConfInvMassLowLimitPhi) || (phiM > config.ConfInvMassUpLimitPhi))) {
        continue;
      }

      PhiCuts.fillQA<aod::femtoworldparticle::ParticleType::kPhi, aod::femtoworldparticle::ParticleType::kPhiChild>(col, p1, p1, p2); ///\todo fill QA also for daughters
      auto cutContainerV0 = PhiCuts.getCutContainer<aod::femtoworldparticle::cutContainerType>(col, p1, p2);
      if (true) { // temporary true value, we are doing simpler version first
        int postrackID = p1.globalIndex();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    p1.pt(),
                    p1.eta(),
                    p1.phi(),
                    p1.p(),
                    mMassOne,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhiChild,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    p1.sign(),
                    p1.beta(),
                    p1.itsChi2NCl(),
                    p1.tpcChi2NCl(),
                    p1.tpcNSigmaKa(),
                    p1.tofNSigmaKa(),
                    (uint8_t)p1.tpcNClsFound(),
                    p1.tpcNClsFindable(),
                    (uint8_t)p1.tpcNClsCrossedRows(),
                    p1.tpcNClsShared(),
                    p1.tpcInnerParam(),
                    p1.itsNCls(),
                    p1.itsNClsInnerBarrel(),
                    p1.dcaXY(),
                    p1.dcaZ(),
                    p1.tpcSignal(),
                    p1.tpcNSigmaStoreEl(),
                    p1.tpcNSigmaStorePi(),
                    p1.tpcNSigmaStoreKa(),
                    p1.tpcNSigmaStorePr(),
                    p1.tpcNSigmaStoreDe(),
                    p1.tofNSigmaStoreEl(),
                    p1.tofNSigmaStorePi(),
                    p1.tofNSigmaStoreKa(),
                    p1.tofNSigmaStorePr(),
                    p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
        const int rowOfPosTrack = outputParts.lastIndex();
        int negtrackID = p2.globalIndex();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;
        outputParts(outputCollision.lastIndex(),
                    p2.pt(),
                    p2.eta(),
                    p2.phi(),
                    p2.p(),
                    mMassTwo,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhiChild,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    p2.sign(),
                    p2.beta(),
                    p2.itsChi2NCl(),
                    p2.tpcChi2NCl(),
                    p2.tpcNSigmaKa(),
                    p2.tofNSigmaKa(),
                    (uint8_t)p2.tpcNClsFound(),
                    p2.tpcNClsFindable(),
                    (uint8_t)p2.tpcNClsCrossedRows(),
                    p2.tpcNClsShared(),
                    p2.tpcInnerParam(),
                    p2.itsNCls(),
                    p2.itsNClsInnerBarrel(),
                    p2.dcaXY(),
                    p2.dcaZ(),
                    p2.tpcSignal(),
                    p2.tpcNSigmaStoreEl(),
                    p2.tpcNSigmaStorePi(),
                    p2.tpcNSigmaStoreKa(),
                    p2.tpcNSigmaStorePr(),
                    p2.tpcNSigmaStoreDe(),
                    p2.tofNSigmaStoreEl(),
                    p2.tofNSigmaStorePi(),
                    p2.tofNSigmaStoreKa(),
                    p2.tofNSigmaStorePr(),
                    p2.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);

        const int rowOfNegTrack = outputParts.lastIndex();
        int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
        float phiPhi = sumVec.Phi();
        if (sumVec.Phi() < 0) {
          phiPhi = sumVec.Phi() + 2 * o2::constants::math::PI;
        } else if (sumVec.Phi() >= 0) {
          phiPhi = sumVec.Phi();
        }
        outputParts(outputCollision.lastIndex(),
                    phiPt,
                    phiEta,
                    phiPhi,
                    phiP,
                    phiM,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhi,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                    0,
                    0, // p1.v0cosPA(),
                    indexChildID,
                    0, // v0.mLambda(),
                    0, // v0.mAntiLambda(),
                    p1.sign(),
                    p1.beta(),
                    p1.itsChi2NCl(),
                    p1.tpcChi2NCl(),
                    p1.tpcNSigmaKa(),
                    p1.tofNSigmaKa(),
                    (uint8_t)p1.tpcNClsFound(),
                    0, // p1.tpcNClsFindable(),
                    0, //(uint8_t)p1.tpcNClsCrossedRows(),
                    p1.tpcNClsShared(),
                    p1.tpcInnerParam(),
                    p1.itsNCls(),
                    p1.itsNClsInnerBarrel(),
                    0, // p1.dcaXY(),
                    0, // p1.dcaZ(),
                    p1.tpcSignal(),
                    p1.tpcNSigmaStoreEl(),
                    p1.tpcNSigmaStorePi(),
                    p1.tpcNSigmaStoreKa(),
                    p1.tpcNSigmaStorePr(),
                    p1.tpcNSigmaStoreDe(),
                    p1.tofNSigmaStoreEl(),
                    p1.tofNSigmaStorePi(),
                    p1.tofNSigmaStoreKa(),
                    p1.tofNSigmaStorePr(),
                    p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackPhiRun2, "Process data", false);

  void processProdTrackD0Run2(aod::FemtoFullCollisionRun2 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  // void processCharmMesons(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centRun2V0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centRun2V0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, mult, spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    /*int nTracks = 0;
    if (col.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) {
          continue;
        }
        nTracks++;
      }
    }*/

    for (const auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    // loop over 2-prong candidates
    for (auto& candidate : candidates) { // selectedD0Candidates

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (config.yCandMax >= 0. && std::abs(hfHelper.yD0(candidate)) > config.yCandMax) {
        continue;
      }

      // we retrieve also the event index from one of the daughters
      auto posDauTrack = candidate.prong0_as<aod::FemtoFullTracks>(); // positive daughter
      auto negDauTrack = candidate.prong1_as<aod::FemtoFullTracks>(); // negative daughter

      if ((candidate.isSelD0() >= config.selectionFlagD0) || (candidate.isSelD0bar() >= config.selectionFlagD0bar)) { // temporary true value, we are doing simpler version first
        int postrackID = candidate.prong0Id();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    candidate.ptProng0(),
                    -999.,                                                                          // eta
                    -999.,                                                                          // phi
                    RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()), // momentum
                    0.,                                                                             // posDauTrack.m(), //general mass
                    -999.,                                                                          // D0mass
                    -999.,                                                                          // D0bar mass
                    candidate.isSelD0(),                                                            // D0 flag
                    candidate.isSelD0bar(),                                                         // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0barChild,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    posDauTrack.sign(),
                    posDauTrack.beta(),
                    -999., // p1.itsChi2NCl(),
                    -999., // p1.tpcChi2NCl(),
                    -999., // p1.tpcNSigmaKa(),
                    -999., // p1.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    posDauTrack.dcaXY(),
                    posDauTrack.dcaZ(),
                    posDauTrack.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tpcNSigmaPi(), //p1.tpcNSigmaStorePi(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tpcNSigmaKa(), //p1.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tofNSigmaPi(), //p1.tofNSigmaStorePi(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tofNSigmaKa(), //p1.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
        const int rowOfPosTrack = outputParts.lastIndex();
        int negtrackID = candidate.prong1Id();
        // int negtrackID = candidate.prong1().globalIndex();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;
        outputParts(outputCollision.lastIndex(),
                    candidate.ptProng1(),
                    0.,                                                                             // negDauTrack.eta(),
                    0.,                                                                             // negDauTrack.phi(),
                    RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()), // negDauTrack.p(),
                    0.,                                                                             // negDauTrack.m(),
                    -999.,                                                                          // D0mass
                    -999.,                                                                          // D0bar mass
                    candidate.isSelD0(),                                                            // D0 flag
                    candidate.isSelD0bar(),                                                         // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0barChild,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    negDauTrack.sign(),
                    negDauTrack.beta(),
                    -999., // negDauTrack.itsChi2NCl(),
                    -999., // negDauTrack.tpcChi2NCl(),
                    -999., // negDauTrack.tpcNSigmaKa(),
                    -999., // negDauTrack.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    negDauTrack.dcaXY(),
                    negDauTrack.dcaZ(),
                    negDauTrack.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // negDauTrack.tpcNSigmaStorePi(),
                    -999., // negDauTrack.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // negDauTrack.tofNSigmaStorePi(),
                    -999., // negDauTrack.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);

        const int rowOfNegTrack = outputParts.lastIndex();
        int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
        outputParts(outputCollision.lastIndex(),
                    candidate.pt(),
                    candidate.eta(),
                    candidate.phi(),
                    candidate.p(),
                    0.,                                    // general mass
                    hfHelper.invMassD0ToPiK(candidate),    // D0mass
                    hfHelper.invMassD0barToKPi(candidate), // D0bar mass
                    candidate.isSelD0(),                   // D0 flag
                    candidate.isSelD0bar(),                // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0bar,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                    0,
                    0, // p1.v0cosPA(),
                    indexChildID,
                    -999,  // v0.mLambda(),
                    -999,  // v0.mAntiLambda(),
                    -999., // p1.sign(),
                    -999., // p1.beta(),
                    -999., // p1.itsChi2NCl(),
                    -999., // p1.tpcChi2NCl(),
                    -999., // p1.tpcNSigmaKa(),
                    -999., // p1.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    -999., // candidate.dcaXY(),
                    -999., // candidate.dcaZ(),
                    -999., // p1.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // p1.tpcNSigmaStorePi(),
                    -999., // p1.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // p1.tofNSigmaStorePi(),
                    -999., // p1.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999., // daughDCA
                    candidate.decayLength(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackD0Run2, "Process data", false);

  void processProdTrackRun3(aod::FemtoFullCollisionRun3 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centFT0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    // const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centFT0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, col.multFT0M(), spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackRun3, "Produce Femto tables", true);

  void processProdTrackV0Run3(aod::FemtoFullCollisionRun3 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks,
                              o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centFT0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    // const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centFT0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, col.multFT0M(), spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    if (config.ConfStoreV0) {
      for (auto& v0 : fullV0s) {
        auto postrack = v0.posTrack_as<aod::FemtoFullTracks>();
        auto negtrack = v0.negTrack_as<aod::FemtoFullTracks>(); ///\tocheck funnily enough if we apply the filter the sign of Pos and Neg track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        v0Cuts.fillLambdaQA(col, v0, postrack, negtrack);

        if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }

        if (config.ConfRejectITSHitandTOFMissing) {
          // Uncomment only when TOF timing is solved
          // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack, TrackSelection::TrackCuts::kITSHits);
          // bool itsHit = o2PhysicsTrackSelection->IsSelected(negtrack, TrackSelection::TrackCuts::kITSHits);
        }

        v0Cuts.fillQA<aod::femtoworldparticle::ParticleType::kV0, aod::femtoworldparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtoworldparticle::cutContainerType>(col, v0, postrack, negtrack);

        if ((cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0) > 0) && (cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts) > 0) && (cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts) > 0)) {
          int postrackID = v0.posTrackId();
          int rowInPrimaryTrackTablePos = -1;
          rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
          childIDs[0] = rowInPrimaryTrackTablePos;
          childIDs[1] = 0;
          outputParts(outputCollision.lastIndex(),
                      v0.positivept(),
                      v0.positiveeta(),
                      v0.positivephi(),
                      0,     // v0.p(),
                      0,     // mass
                      -999., // D0mass
                      -999., // D0bar mass
                      -999., // D0 flag
                      -999., // D0bar flag
                      aod::femtoworldparticle::ParticleType::kV0Child,
                      cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                      cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                      0.,
                      childIDs,
                      0,
                      0,
                      postrack.sign(),
                      postrack.beta(),
                      postrack.itsChi2NCl(),
                      postrack.tpcChi2NCl(),
                      postrack.tpcNSigmaKa(),
                      postrack.tofNSigmaKa(),
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
                      -999.,
                      -999.,
                      -999.);
          const int rowOfPosTrack = outputParts.lastIndex();
          int negtrackID = v0.negTrackId();
          int rowInPrimaryTrackTableNeg = -1;
          rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNeg;
          outputParts(outputCollision.lastIndex(),
                      v0.negativept(),
                      v0.negativeeta(),
                      v0.negativephi(),
                      0,     // momentum
                      0,     // mass
                      -999., // D0mass
                      -999., // D0bar mass
                      -999., // D0 flag
                      -999., // D0bar flag
                      aod::femtoworldparticle::ParticleType::kV0Child,
                      cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts),
                      cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegPID),
                      0.,
                      childIDs,
                      0,
                      0,
                      negtrack.sign(),
                      negtrack.beta(),
                      negtrack.itsChi2NCl(),
                      negtrack.tpcChi2NCl(),
                      negtrack.tpcNSigmaKa(),
                      negtrack.tofNSigmaKa(),
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
                      -999.,
                      -999.,
                      -999.);
          const int rowOfNegTrack = outputParts.lastIndex();
          int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
          outputParts(outputCollision.lastIndex(),
                      v0.pt(),
                      v0.eta(),
                      v0.phi(),
                      0,     // momentum
                      0,     // mass
                      -999., // D0mass
                      -999., // D0bar mass
                      -999., // D0 flag
                      -999., // D0bar flag
                      aod::femtoworldparticle::ParticleType::kV0,
                      cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                      0,
                      v0.v0cosPA(),
                      indexChildID,
                      v0.mLambda(),
                      v0.mAntiLambda(),
                      postrack.sign(),
                      postrack.beta(),
                      postrack.itsChi2NCl(),
                      postrack.tpcChi2NCl(),
                      postrack.tpcNSigmaKa(),
                      postrack.tofNSigmaKa(),
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
                      -999.,
                      -999.,
                      -999.);
        }
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackV0Run3, "Produce Femto tables", true);

  void processProdTrackPhiRun3(aod::FemtoFullCollisionRun3 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, o2::aod::V0Datas const&)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centFT0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    // const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centFT0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, col.multFT0M(), spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      if ((p1.trackType() == o2::aod::track::TrackTypeEnum::Run2Tracklet) || (p2.trackType() == o2::aod::track::TrackTypeEnum::Run2Tracklet)) {
        continue;
      } else if (p1.globalIndex() == p2.globalIndex()) { // checking not to correlate same particles
        continue;
      } else if ((p1.pt() < config.cfgPtLowPart1) || (p1.pt() > config.cfgPtHighPart1)) { // pT cuts for part1
        continue;
      } else if ((p1.p() < config.cfgPLowPart1) || (p1.p() > config.cfgPHighPart1)) { // p cuts for part1
        continue;
      } else if ((p1.eta() < config.cfgEtaLowPart1) || (p1.eta() > config.cfgEtaHighPart1)) { // eta cuts for part1
        continue;
      } else if ((p2.pt() < config.cfgPtLowPart2) || (p2.pt() > config.cfgPtHighPart2)) { // pT cuts for part2
        continue;
      } else if ((p2.p() < config.cfgPLowPart2) || (p2.p() > config.cfgPHighPart2)) { // p cuts for part2
        continue;
      } else if ((p2.eta() < config.cfgEtaLowPart2) || (p2.eta() > config.cfgEtaHighPart2)) { // eta for part2
        continue;
      } else if (!(IsKaonNSigma(p1.p(), p1.tpcNSigmaKa(), p1.tofNSigmaKa()))) { // PID for Kaons
        continue;
      } else if (!(IsKaonNSigma(p2.p(), p2.tpcNSigmaKa(), p2.tofNSigmaKa()))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;
      float mMassOne = TDatabasePDG::Instance()->GetParticle(config.ConfPDGCodePartOne)->Mass(); // FIXME: Get from the PDG service of the common header
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(config.ConfPDGCodePartTwo)->Mass(); // FIXME: Get from the PDG service of the common header

      part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      float phiPt = sumVec.Pt();
      float phiP = sumVec.P();
      float phiM = sumVec.M();

      if (((phiM < config.ConfInvMassLowLimitPhi) || (phiM > config.ConfInvMassUpLimitPhi))) {
        continue;
      }

      PhiCuts.fillQA<aod::femtoworldparticle::ParticleType::kPhi, aod::femtoworldparticle::ParticleType::kPhiChild>(col, p1, p1, p2); ///\todo fill QA also for daughters
      auto cutContainerV0 = PhiCuts.getCutContainer<aod::femtoworldparticle::cutContainerType>(col, p1, p2);
      if (true) { // temporary true value, we are doing simpler version first
        int postrackID = p1.globalIndex();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    p1.pt(),
                    p1.eta(),
                    p1.phi(),
                    p1.p(),
                    mMassOne,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhiChild,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    p1.sign(),
                    p1.beta(),
                    p1.itsChi2NCl(),
                    p1.tpcChi2NCl(),
                    p1.tpcNSigmaKa(),
                    p1.tofNSigmaKa(),
                    (uint8_t)p1.tpcNClsFound(),
                    p1.tpcNClsFindable(),
                    (uint8_t)p1.tpcNClsCrossedRows(),
                    p1.tpcNClsShared(),
                    p1.tpcInnerParam(),
                    p1.itsNCls(),
                    p1.itsNClsInnerBarrel(),
                    p1.dcaXY(),
                    p1.dcaZ(),
                    p1.tpcSignal(),
                    p1.tpcNSigmaStoreEl(),
                    p1.tpcNSigmaStorePi(),
                    p1.tpcNSigmaStoreKa(),
                    p1.tpcNSigmaStorePr(),
                    p1.tpcNSigmaStoreDe(),
                    p1.tofNSigmaStoreEl(),
                    p1.tofNSigmaStorePi(),
                    p1.tofNSigmaStoreKa(),
                    p1.tofNSigmaStorePr(),
                    p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
        const int rowOfPosTrack = outputParts.lastIndex();
        int negtrackID = p2.globalIndex();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;
        outputParts(outputCollision.lastIndex(),
                    p2.pt(),
                    p2.eta(),
                    p2.phi(),
                    p2.p(),
                    mMassTwo,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhiChild,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegCuts),
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kNegPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    p2.sign(),
                    p2.beta(),
                    p2.itsChi2NCl(),
                    p2.tpcChi2NCl(),
                    p2.tpcNSigmaKa(),
                    p2.tofNSigmaKa(),
                    (uint8_t)p2.tpcNClsFound(),
                    p2.tpcNClsFindable(),
                    (uint8_t)p2.tpcNClsCrossedRows(),
                    p2.tpcNClsShared(),
                    p2.tpcInnerParam(),
                    p2.itsNCls(),
                    p2.itsNClsInnerBarrel(),
                    p2.dcaXY(),
                    p2.dcaZ(),
                    p2.tpcSignal(),
                    p2.tpcNSigmaStoreEl(),
                    p2.tpcNSigmaStorePi(),
                    p2.tpcNSigmaStoreKa(),
                    p2.tpcNSigmaStorePr(),
                    p2.tpcNSigmaStoreDe(),
                    p2.tofNSigmaStoreEl(),
                    p2.tofNSigmaStorePi(),
                    p2.tofNSigmaStoreKa(),
                    p2.tofNSigmaStorePr(),
                    p2.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);

        const int rowOfNegTrack = outputParts.lastIndex();
        int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
        float phiPhi = sumVec.Phi();
        if (sumVec.Phi() < 0) {
          phiPhi = sumVec.Phi() + 2 * o2::constants::math::PI;
        } else if (sumVec.Phi() >= 0) {
          phiPhi = sumVec.Phi();
        }
        outputParts(outputCollision.lastIndex(),
                    phiPt,
                    phiEta,
                    phiPhi,
                    phiP,
                    phiM,
                    -999., // D0mass
                    -999., // D0bar mass
                    -999., // D0 flag
                    -999., // D0bar flag
                    aod::femtoworldparticle::ParticleType::kPhi,
                    cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                    0,
                    0, // p1.v0cosPA(),
                    indexChildID,
                    0, // v0.mLambda(),
                    0, // v0.mAntiLambda(),
                    p1.sign(),
                    p1.beta(),
                    p1.itsChi2NCl(),
                    p1.tpcChi2NCl(),
                    p1.tpcNSigmaKa(),
                    p1.tofNSigmaKa(),
                    (uint8_t)p1.tpcNClsFound(),
                    0, // p1.tpcNClsFindable(),
                    0, //(uint8_t)p1.tpcNClsCrossedRows(),
                    p1.tpcNClsShared(),
                    p1.tpcInnerParam(),
                    p1.itsNCls(),
                    p1.itsNClsInnerBarrel(),
                    0, // p1.dcaXY(),
                    0, // p1.dcaZ(),
                    p1.tpcSignal(),
                    p1.tpcNSigmaStoreEl(),
                    p1.tpcNSigmaStorePi(),
                    p1.tpcNSigmaStoreKa(),
                    p1.tpcNSigmaStorePr(),
                    p1.tpcNSigmaStoreDe(),
                    p1.tofNSigmaStoreEl(),
                    p1.tofNSigmaStorePi(),
                    p1.tofNSigmaStoreKa(),
                    p1.tofNSigmaStorePr(),
                    p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackPhiRun3, "Process data", false);

  void processProdTrackD0Run3(aod::FemtoFullCollisionRun3 const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  // void processCharmMesons(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (config.ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFT0M(), colCuts.computeSphericity(col, tracks), mMagField, col.centFT0M());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    // const auto mult = collision.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    const auto centrality = col.centFT0M();
    colCuts.fillQA(col);

    outputCollision(vtxZ, col.multFT0M(), spher, mMagField, centrality);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (const auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  track.p(),
                  1,
                  -999., // D0mass
                  -999., // D0bar mass
                  -999., // D0 flag
                  -999., // D0bar flag
                  aod::femtoworldparticle::ParticleType::kTrack,
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                  0, // cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0, // początek nowej części
                  track.sign(),
                  track.beta(),
                  track.itsChi2NCl(),
                  track.tpcChi2NCl(),
                  track.tpcNSigmaKa(),
                  track.tofNSigmaKa(),
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
                  -999.,
                  -999.,
                  -999.);
      tmpIDtrack.push_back(track.globalIndex());
    }

    // loop over 2-prong candidates
    for (const auto& candidate : candidates) { // selectedD0Candidates
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (config.yCandMax >= 0. && std::abs(hfHelper.yD0(candidate)) > config.yCandMax) {
        continue;
      }

      // we retrieve also the event index from one of the daughters
      auto posDauTrack = candidate.prong0_as<aod::FemtoFullTracks>(); // positive daughter
      auto negDauTrack = candidate.prong1_as<aod::FemtoFullTracks>(); // negative daughter
      // candidate.prong1_as<aod::BigTracksPID>().tofNSigmaPi()
      // candidate.prong1_as<aod::BigTracksPID>().tofNSigmaKa()

      if ((candidate.isSelD0() == 1) || (candidate.isSelD0bar() == 1)) { // temporary true value, we are doing simpler version first
        int postrackID = candidate.prong0Id();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    candidate.ptProng0(),
                    -999.,                                                                          // eta
                    -999.,                                                                          // phi
                    RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()), // momentum
                    0.,                                                                             // posDauTrack.m(), //general mass
                    -999.,                                                                          // D0mass
                    -999.,                                                                          // D0bar mass
                    candidate.isSelD0(),                                                            // D0 flag
                    candidate.isSelD0bar(),                                                         // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0barChild,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    posDauTrack.sign(),
                    posDauTrack.beta(),
                    -999., // p1.itsChi2NCl(),
                    -999., // p1.tpcChi2NCl(),
                    -999., // p1.tpcNSigmaKa(),
                    -999., // p1.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    posDauTrack.dcaXY(),
                    posDauTrack.dcaZ(),
                    posDauTrack.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tpcNSigmaPi(), //p1.tpcNSigmaStorePi(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tpcNSigmaKa(), //p1.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tofNSigmaPi(), //p1.tofNSigmaStorePi(),
                    -999., // candidate.prong0_as<aod::BigTracksPID>().tofNSigmaKa(), //p1.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
        const int rowOfPosTrack = outputParts.lastIndex();
        int negtrackID = candidate.prong1Id();
        // int negtrackID = candidate.prong1().globalIndex();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;
        outputParts(outputCollision.lastIndex(),
                    candidate.ptProng1(),
                    0.,                                                                             // negDauTrack.eta(),
                    0.,                                                                             // negDauTrack.phi(),
                    RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()), // negDauTrack.p(),
                    0.,                                                                             // negDauTrack.m(),
                    -999.,                                                                          // D0mass
                    -999.,                                                                          // D0bar mass
                    candidate.isSelD0(),                                                            // D0 flag
                    candidate.isSelD0bar(),                                                         // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0barChild,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosCuts),
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kPosPID),
                    0.,
                    childIDs,
                    0,
                    0,
                    negDauTrack.sign(),
                    negDauTrack.beta(),
                    -999., // negDauTrack.itsChi2NCl(),
                    -999., // negDauTrack.tpcChi2NCl(),
                    -999., // negDauTrack.tpcNSigmaKa(),
                    -999., // negDauTrack.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    negDauTrack.dcaXY(),
                    negDauTrack.dcaZ(),
                    negDauTrack.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // negDauTrack.tpcNSigmaStorePi(),
                    -999., // negDauTrack.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // negDauTrack.tofNSigmaStorePi(),
                    -999., // negDauTrack.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);

        const int rowOfNegTrack = outputParts.lastIndex();
        int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
        outputParts(outputCollision.lastIndex(),
                    candidate.pt(),
                    candidate.eta(),
                    candidate.phi(),
                    candidate.p(),
                    0.,                                    // general mass
                    hfHelper.invMassD0ToPiK(candidate),    // D0mass
                    hfHelper.invMassD0barToKPi(candidate), // D0bar mass
                    candidate.isSelD0(),                   // D0 flag
                    candidate.isSelD0bar(),                // D0bar flag
                    aod::femtoworldparticle::ParticleType::kD0D0bar,
                    0., // cutContainerV0.at(femtoWorldV0Selection::V0ContainerPosition::kV0),
                    0,
                    0, // p1.v0cosPA(),
                    indexChildID,
                    -999,  // v0.mLambda(),
                    -999,  // v0.mAntiLambda(),
                    -999., // p1.sign(),
                    -999., // p1.beta(),
                    -999., // p1.itsChi2NCl(),
                    -999., // p1.tpcChi2NCl(),
                    -999., // p1.tpcNSigmaKa(),
                    -999., // p1.tofNSigmaKa(),
                    -999., //(uint8_t)p1.tpcNClsFound(),
                    -999., // p1.tpcNClsFindable(),
                    -999., //(uint8_t)p1.tpcNClsCrossedRows(),
                    -999., // p1.tpcNClsShared(),
                    -999., // p1.tpcInnerParam(),
                    -999., // p1.itsNCls(),
                    -999., // p1.itsNClsInnerBarrel(),
                    -999., // candidate.dcaXY(),
                    -999., // candidate.dcaZ(),
                    -999., // p1.tpcSignal(),
                    -999., // p1.tpcNSigmaStoreEl(),
                    -999., // p1.tpcNSigmaStorePi(),
                    -999., // p1.tpcNSigmaStoreKa(),
                    -999., // p1.tpcNSigmaStorePr(),
                    -999., // p1.tpcNSigmaStoreDe(),
                    -999., // p1.tofNSigmaStoreEl(),
                    -999., // p1.tofNSigmaStorePi(),
                    -999., // p1.tofNSigmaStoreKa(),
                    -999., // p1.tofNSigmaStorePr(),
                    -999., // p1.tofNSigmaStoreDe(),
                    -999., // daughDCA
                    candidate.decayLength(),
                    -999.,
                    -999.,
                    -999.,
                    -999.,
                    -999.);
      }
    }
  }
  PROCESS_SWITCH(femtoWorldProducerTask, processProdTrackD0Run3, "Process data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoWorldProducerTask>(cfgc)};
  return workflow;
}
