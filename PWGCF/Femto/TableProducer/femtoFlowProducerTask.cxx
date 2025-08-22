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
/// \author Wenya Wu, TUM, wenya.wu@cern.ch
/// \NOTE:  The femtoflow framework borrows and copies the framework of femtouniverse and femtodream

#include "PWGCF/Femto/Core/FemtoFlowCollisionSelection.h"
// #include "PWGCF/Femto/Core/FemtoFlowPhiSelection.h"
#include "PWGCF/Femto/Core/FemtoFlowTrackSelection.h"
// #include "PWGCF/Femto/Core/FemtoFlowV0Selection.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/DataModel/FemtoDerived.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
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
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_flow;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollision_CentPbPb =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFV0As, aod::QvectorFT0CVecs>::iterator;
using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA, aod::TOFSignal, aod::pidTOFbeta, aod::pidTPCEl, aod::TrackSelection,
            aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
            aod::pidTPCDe, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
            aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;
} // namespace o2::aod

namespace softwareTriggers
{
static const int nTriggers = 6;
static const std::vector<std::string> triggerNames{"fPPP", "fPPL", "fPLL", "fLLL", "fPD", "fLD"};
static const float triggerSwitches[1][nTriggers]{
  {0, 0, 0, 0, 0, 0}};
} // namespace softwareTriggers

/// \todo fix how to pass array to setSelection, getRow() passing a different
/// type!
// static constexpr float arrayV0Sel[3][3] = {{100.f, 100.f, 100.f}, {0.2f,
// 0.2f, 0.2f}, {100.f, 100.f, 100.f}}; unsigned int rows = sizeof(arrayV0Sel) /
// sizeof(arrayV0Sel[0]); unsigned int columns = sizeof(arrayV0Sel[0]) /
// sizeof(arrayV0Sel[0][0]);

struct FemtoFlowProducerTask {
  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDExtParticles> outputDebugParts;

  Configurable<bool> confIsDebug{"confIsDebug", true, "Enable Debug tables"};
  Configurable<bool> confIsRun3{"confIsRun3", true, "Running on Run3 or pilot"};
  Configurable<bool> confDoSpher{"confDoSpher", true, "Calculate sphericity. If false sphericity will take value of 2."};
  Configurable<bool> confDoqnVec{"confDoqnVec", true, "Calculate qn vector. If false sphericity will take value of 10."};
  Configurable<bool> confIsForceGRP{"confIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  /// Event cuts
  FemtoFlowCollisionSelection colCuts;
  // Event cuts - Triggers
  Configurable<bool> ConfEnableTriggerSelection{"ConfEnableTriggerSelection", false, "Should the trigger selection be enabled for collisions?"};
  Configurable<LabeledArray<float>> ConfTriggerSwitches{
    "ConfTriggerSwitches",
    {softwareTriggers::triggerSwitches[0], 1, softwareTriggers::nTriggers, std::vector<std::string>{"Switch"}, softwareTriggers::triggerNames},
    "Turn on which trigger should be checked for recorded events to pass selection"};
  Configurable<std::string> ConfBaseCCDBPathForTriggers{"ConfBaseCCDBPathForTriggers", "Users/m/mpuccio/EventFiltering/OTS/Chunked/", "Provide ccdb path for trigger table; default - trigger coordination"};

  // Event cuts - usual selection criteria
  Configurable<bool> confEvtUseTPCmult{"confEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<float> confCentFT0Min{"confCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> confCentFT0Max{"confCentFT0Max", 100.f, "Max CentFT0 value for centrality selection"};
  Configurable<bool> confIsUsePileUp{"confIsUsePileUp", true, "Required for choosing whether to run the pile-up cuts"};
  Configurable<bool> confEvNoSameBunchPileup{"confEvNoSameBunchPileup", true, "Require kNoSameBunchPileup selection on Events."};
  Configurable<bool> confEvIsGoodZvtxFT0vsPV{"confEvIsGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV selection on Events."};
  Configurable<bool> confEvIsGoodITSLayersAll{"confEvIsGoodITSLayersAll", true, "Require kIsGoodITSLayersAll selection on Events."};
  Configurable<bool> confEvNoCollInRofStandard{"confEvNoCollInRofStandard", true, "Require kNoCollInRofStandard selection on Events."};
  Configurable<bool> confEvNoHighMultCollInPrevRof{"confEvNoHighMultCollInPrevRof", true, "Require kNoHighMultCollInPrevRof selection on Events."};
  Configurable<bool> confEvNoCollInTimeRangeStandard{"confEvNoCollInTimeRangeStandard", true, "Require kNoCollInTimeRangeStandard selection on Events."};
  Configurable<bool> confEvIsVertexITSTPC{"confEvIsVertexITSTPC", true, "Require kIsVertexITSTPC selection on Events"};
  Configurable<int> confTPCOccupancyMin{"confTPCOccupancyMin", 0, "Minimum value for TPC Occupancy selection"};
  Configurable<int> confTPCOccupancyMax{"confTPCOccupancyMax", 1000, "Maximum value for TPC Occupancy selection"};

  // Track cuts
  FemtoFlowTrackSelection trackCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confTrkCharge{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kSign, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmin{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kpTMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmax{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kpTMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkEta{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kEtaMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCnclsMin{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kTPCnClsMin, "ConfTrk"), std::vector<float>{70.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kTPCnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCfCls{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.83f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kTPCfClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCcRowsMin{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kTPCcRowsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCsCls{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kTPCsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCfracsCls{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kTPCfracsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kTPCfracsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsMin{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kITSnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsIbMin{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kITSnClsIbMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAxyMax{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kDCAxyMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAzMax{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kDCAzMax, "Track selection: ")}; /// \todo Reintegrate PID to the general selection container
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{FemtoFlowTrackSelection::getSelectionName(femto_flow_track_selection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoFlowTrackSelection::getSelectionHelper(femto_flow_track_selection::kPIDnSigmaMax, "Track selection: ")};
    // Configurable<std::vector<int>> confTrkPIDspecies{"confTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID (Pion=2, Kaon=3, Proton=4, Deuteron=5)"};
    Configurable<std::vector<int>> confTrkPIDspecies{"confTrkPIDspecies", std::vector<int>{o2::track::PID::Proton}, "Trk sel: Only select track Proton"};
    // Numbers from ~/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h //static constexpr ID Pion = 2; static constexpr ID Kaon = 3; static constexpr ID Proton = 4; static constexpr ID Deuteron = 5;
  } ConfTrkSelection;

  Configurable<float> confTrkPIDnSigmaOffsetTPC{"confTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> confTrkPIDnSigmaOffsetTOF{"confTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<float> confTOFpTmin{"confTOFpTmin", 500, "TOF pT min"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfDcaXYCustom0Cut{"ConfDcaXYCustom0Cut", false, "Enable Custom Dcaxy < [0] cut."};
    Configurable<float> ConfDcaXYFilterCut{"ConfDcaXYFilterCut", 2.4, "Value for DCA_XY for the global track"}; // max dca to vertex XY
    Configurable<bool> ConfDcaXYCustom1Cut{"ConfDcaXYCustom1Cut", true, "Enable Custom |DCAxy| < [1] + [2]/pt cut."};
    Configurable<float> ConfDcaXYCustom11FilterCut{"ConfDcaXYCustom11FilterCut", 0.004, "Value for [1] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> ConfDcaXYCustom12FilterCut{"ConfDcaXYCustom12FilterCut", 0.013, "Value for [2] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> ConfDcaZFilterCut{"ConfDcaZFilterCut", 3.2, "Value for DCA_Z for the global track"}; // max dca to vertex Z
    Configurable<float> ConfTrkMinChi2PerClusterTPC{"ConfTrkMinChi2PerClusterTPC", 0.f, "Lower limit for chi2 of TPC; currently for testing only"};
    Configurable<float> ConfTrkMaxChi2PerClusterTPC{"ConfTrkMaxChi2PerClusterTPC", 4.f, "Upper limit for chi2 of TPC; currently for testing only"};
    Configurable<float> ConfTrkMaxChi2PerClusterITS{"ConfTrkMaxChi2PerClusterITS", 10.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default
    Configurable<bool> ConfTrkTPCRefit{"ConfTrkTPCRefit", false, "True: require TPC refit"};
    Configurable<bool> ConfTrkITSRefit{"ConfTrkITSRefit", false, "True: require ITS refit"};
  } ConfTrackSpecialFilters;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  HistogramRegistry TrackRegistry{"Tracks", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber = 0;
  float mMagField;
  std::string zorroTriggerNames = "";
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    colCuts.setCuts(confEvtZvtx, confEvtTriggerCheck, confEvtTriggerSel, confEvtOfflineCheck, confIsRun3, confCentFT0Min, confCentFT0Max);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkSelection.confTrkCharge, femto_flow_track_selection::kSign, femto_flow_selection::kEqual);
    trackCuts.setSelection(ConfTrkSelection.confTrkPtmin, femto_flow_track_selection::kpTMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkPtmax, femto_flow_track_selection::kpTMax, femto_flow_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkEta, femto_flow_track_selection::kEtaMax, femto_flow_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCnclsMin, femto_flow_track_selection::kTPCnClsMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCfCls, femto_flow_track_selection::kTPCfClsMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCcRowsMin, femto_flow_track_selection::kTPCcRowsMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCsCls, femto_flow_track_selection::kTPCsClsMax, femto_flow_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCfracsCls, femto_flow_track_selection::kTPCfracsClsMax, femto_flow_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkITSnclsMin, femto_flow_track_selection::kITSnClsMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkITSnclsIbMin, femto_flow_track_selection::kITSnClsIbMin, femto_flow_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkDCAxyMax, femto_flow_track_selection::kDCAxyMax, femto_flow_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkDCAzMax, femto_flow_track_selection::kDCAzMax, femto_flow_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkPIDnSigmaMax, femto_flow_track_selection::kPIDnSigmaMax, femto_flow_selection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkSelection.confTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(confTrkPIDnSigmaOffsetTPC, confTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtoflowparticle::ParticleType::kTrack, aod::femtoflowparticle::TrackType::kNoChild, aod::femtoflowparticle::CutContainerType>(&TrackRegistry);

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

  template <typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    outputDebugParts(particle.sign(), (uint8_t)particle.tpcNClsFound(),
                     particle.tpcNClsFindable(),
                     (uint8_t)particle.tpcNClsCrossedRows(),
                     particle.tpcNClsShared(), particle.tpcFractionSharedCls(), particle.tpcInnerParam(),
                     particle.itsNCls(), particle.itsNClsInnerBarrel(),
                     particle.dcaXY(), particle.dcaZ(), particle.tpcSignal(), particle.beta(),
                     particle.tpcNSigmaStoreEl(), particle.tpcNSigmaStorePi(),
                     particle.tpcNSigmaStoreKa(), particle.tpcNSigmaStorePr(),
                     particle.tpcNSigmaStoreDe(), particle.tofNSigmaStoreEl(),
                     particle.tofNSigmaStorePi(), particle.tofNSigmaStoreKa(),
                     particle.tofNSigmaStorePr(), particle.tofNSigmaStoreDe(),
                     -999., -999., -999., -999., -999., -999.);
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
                                     /// FemtoFlowRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (confEvtUseTPCmult) {
      multNtr = col.multTPC();
    }
    const auto occupancy = col.trackOccupancyInTimeRange();

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col) || (occupancy < confTPCOccupancyMin || occupancy > confTPCOccupancyMax)) {
      return false;
    }

    float sphericity = colCuts.computeSphericity(col, tracks);
    int qnbin = colCuts.myqnBin(col);

    if (!confIsUsePileUp) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? sphericity : 2, (confDoqnVec && qnbin >= 0 && qnbin < 10) ? qnbin : -999, mMagField);
      colCuts.fillQA(col);
      return true;
    } else if ((!confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!confEvIsGoodITSLayersAll || col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) && (!confEvNoCollInRofStandard || col.selection_bit(aod::evsel::kNoCollInRofStandard)) && (!confEvNoHighMultCollInPrevRof || col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) && (!confEvNoCollInTimeRangeStandard || col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
               // && (!confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))
    ) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? sphericity : 2, (confDoqnVec && qnbin >= 0 && qnbin < 10) ? qnbin : -999, mMagField);
      colCuts.fillQA(col);
      return true;
    } else {
      return false;
    }
  }

  // template <bool isMC, typename CollisionType>
  // bool fillCollisionsCentRun3(CollisionType const& col)
  // {
  //   const auto vtxZ = col.posZ();
  //   const auto multNtr = col.multNTracksPV();
  //   const auto cent = col.centFT0C();
  //   const auto occupancy = col.trackOccupancyInTimeRange();

  //   // check whether the basic event selection criteria are fulfilled
  //   // if the basic selection is NOT fulfilled:
  //   // in case of skimming run - don't store such collisions
  //   // in case of trigger run - store such collisions but don't store any
  //   // particle candidates for such collisions

  //   if (!colCuts.isSelectedRun3(col) || (occupancy < confTPCOccupancyMin || occupancy > confTPCOccupancyMax)) {
  //     return false;
  //   } else {
  //     if (col.selection_bit(aod::evsel::kNoSameBunchPileup) &&
  //         col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
  //         col.selection_bit(aod::evsel::kIsGoodITSLayersAll) &&
  //         col.selection_bit(aod::evsel::kNoCollInRofStandard) &&
  //         col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof) &&
  //         col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) {
  //       outputCollision(vtxZ, cent, multNtr, 2, mMagField);
  //       return true;
  //     } else {
  //       return false;
  //     }
  //   }
  // }

  template <bool isMC, typename TrackType>
  void fillTracks(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (const auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (track.pt() > confTOFpTmin) {
        if (!track.hasTOF()) {
          continue;
        }
      }
      if (ConfTrackSpecialFilters.ConfDcaXYCustom0Cut && fabs(track.dcaXY()) > ConfTrackSpecialFilters.ConfDcaXYFilterCut) {
        continue;
      }
      if (ConfTrackSpecialFilters.ConfDcaXYCustom1Cut && fabs(track.dcaXY()) > ConfTrackSpecialFilters.ConfDcaXYCustom11FilterCut + ConfTrackSpecialFilters.ConfDcaXYCustom12FilterCut / track.pt()) {
        continue;
      }
      if (fabs(track.dcaZ()) > ConfTrackSpecialFilters.ConfDcaZFilterCut) {
        continue;
      }
      if (track.tpcChi2NCl() < ConfTrackSpecialFilters.ConfTrkMinChi2PerClusterTPC || track.tpcChi2NCl() > ConfTrackSpecialFilters.ConfTrkMaxChi2PerClusterTPC) {
        continue;
      }
      if (track.itsChi2NCl() > ConfTrackSpecialFilters.ConfTrkMaxChi2PerClusterITS) {
        continue;
      }
      if ((ConfTrackSpecialFilters.ConfTrkTPCRefit && !track.hasTPC()) || (ConfTrackSpecialFilters.ConfTrkITSRefit && !track.hasITS())) {
        continue;
      }
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      trackCuts.fillQA<aod::femtoflowparticle::ParticleType::kTrack,
                       aod::femtoflowparticle::TrackType::kNoChild>(track);

      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtoflowparticle::CutContainerType>(track);

      // now the table is filled
      outputParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                  track.phi(), aod::femtoflowparticle::ParticleType::kTrack,
                  cutContainer.at(
                    femto_flow_track_selection::TrackContainerPosition::kCuts),
                  cutContainer.at(
                    femto_flow_track_selection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0,
                  track.sign()); // sign getter is mAntiLambda()

      if (confIsDebug) {
        fillDebugParticle(track);
      }

      tmpIDtrack.push_back(track.globalIndex());
    }
  }

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillCollisionsAndTracks(CollisionType const& col, TrackType const& tracks)
  {
    const auto colcheck = fillCollisions<isMC>(col, tracks);
    if (colcheck)
      fillTracks<isMC>(tracks);
  }

  void processFullData(aod::FemtoFullCollision_CentPbPb const& col,
                       aod::BCsWithTimestamps const&,
                       aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracks<false>(col, tracks);
  }
  PROCESS_SWITCH(FemtoFlowProducerTask, processFullData, "Provd experimental data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoFlowProducerTask>(cfgc)};
  return workflow;
}
