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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/VarManager.h"

namespace o2::aod
{
namespace dqcuts
{
AnalysisCompositeCut* GetCompositeCut(const char* cutName);
AnalysisCut* GetAnalysisCut(const char* cutName);
} // namespace dqcuts
} // namespace o2::aod

AnalysisCompositeCut* o2::aod::dqcuts::GetCompositeCut(const char* cutName)
{
  //
  // define composie cuts, typically combinations of all the ingredients needed for a full cut
  //
  // TODO: Agree on some conventions for the naming
  //       Think of possible customization of the predefined cuts via names

  AnalysisCompositeCut* cut = new AnalysisCompositeCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("jpsiO2MCdebugCuts")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1"));
    return cut;
  }

  if (!nameStr.compare("jpsiBenchmarkCuts")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityBenchmark"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut1")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine"));
    cut->AddCut(GetAnalysisCut("SPDany"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine"));
    cut->AddCut(GetAnalysisCut("SPDany"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut1SPDfirst")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2SPDfirst")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts2")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts3")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts4")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaLoose"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts5")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryLoose"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts6")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryVeryLoose"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts7")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts8")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts9")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftDown"));
    return cut;
  }

  if (!nameStr.compare("jpsiKineAndQuality")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("jpsiPID1")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1"));
    return cut;
  }

  if (!nameStr.compare("jpsiPID2")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID2"));
    return cut;
  }

  if (!nameStr.compare("jpsiPIDnsigma")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (!nameStr.compare("PIDCalibElectron")) {
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (!nameStr.compare("PIDCalibPion")) {
    cut->AddCut(GetAnalysisCut("pidcalib_pion"));
    return cut;
  }

  if (!nameStr.compare("PIDCalibProton")) {
    cut->AddCut(GetAnalysisCut("pidcalib_proton"));
    return cut;
  }

  if (!nameStr.compare("PIDCalib_basic")) {
    cut->AddCut(GetAnalysisCut("pidbasic"));
    return cut;
  }

  if (!nameStr.compare("highPtHadron")) {
    cut->AddCut(GetAnalysisCut("highPtHadron"));
    return cut;
  }

  if (!nameStr.compare("PIDCalib")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (!nameStr.compare("NoPID")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("kaonPID")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  if (!nameStr.compare("kaonPID")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  //---------------------------------------------------------------------------------------
  // NOTE: Below there are several TPC pid cuts used for studies of the dE/dx degradation
  //    and its impact on the high lumi pp quarkonia triggers
  //  To be removed when not needed anymore
  if (!nameStr.compare("jpsiPID1Randomized")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1randomized"));
    return cut;
  }

  if (!nameStr.compare("jpsiPID2Randomized")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID2randomized"));
    return cut;
  }

  if (!nameStr.compare("jpsiPIDnsigmaRandomized")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaRandomized"));
    return cut;
  }

  if (!nameStr.compare("jpsiPIDworseRes")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDworseRes"));
    return cut;
  }

  if (!nameStr.compare("jpsiPIDshift")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDshift"));
    return cut;
  }

  if (!nameStr.compare("jpsiPID1shiftUp")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (!nameStr.compare("jpsiPID1shiftDown")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftDown"));
    return cut;
  }
  // -------------------------------------------------------------------------------------------------

  if (!nameStr.compare("lmeePID_TPChadrejTOFrec")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));

    AnalysisCompositeCut* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_hadrej);
    cut_pid_OR->AddCut(cut_tof_rec);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmeePID_TPChadrejTOFrecRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_muon_band_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection_highp"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));

    AnalysisCompositeCut* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron_loose"));
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_pion_rejection_highp"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_hadrej);
    cut_pid_OR->AddCut(cut_tof_rec);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmeePID_TPChadrej")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));
    cut->AddCut(cut_tpc_hadrej);
    return cut;
  }

  if (!nameStr.compare("lmee_eNSigma")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmeePID_TOFrec")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    AnalysisCompositeCut* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron"));

    cut->AddCut(cut_tof_rec);
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrack")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_TPCTrackRun3_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3_TPC_ePID_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("tpc_electron"));
    return cut;
  }

  if (!nameStr.compare("muonQualityCuts")) {
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt2")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonHighPt")) {
    cut->AddCut(GetAnalysisCut("muonHighPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonTightQualityCutsForTests")) {
    cut->AddCut(GetAnalysisCut("muonTightQualityCutsForTests"));
    return cut;
  }

  if (!nameStr.compare("mchTrack")) {
    cut->AddCut(GetAnalysisCut("mchTrack"));
    return cut;
  }

  if (!nameStr.compare("matchedMchMid")) {
    cut->AddCut(GetAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (!nameStr.compare("matchedFwd")) {
    cut->AddCut(GetAnalysisCut("matchedFwd"));
    return cut;
  }

  if (!nameStr.compare("matchedGlobal")) {
    cut->AddCut(GetAnalysisCut("matchedGlobal"));
    return cut;
  }

  if (!nameStr.compare("pairNoCut")) {
    cut->AddCut(GetAnalysisCut("pairNoCut"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow1")) {
    cut->AddCut(GetAnalysisCut("pairMassLow1"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow2")) {
    cut->AddCut(GetAnalysisCut("pairMassLow2"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow3")) {
    cut->AddCut(GetAnalysisCut("pairMassLow3"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz1")) {
    cut->AddCut(GetAnalysisCut("pairDalitz1"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz2")) {
    cut->AddCut(GetAnalysisCut("pairDalitz2"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz3")) {
    cut->AddCut(GetAnalysisCut("pairDalitz3"));
    return cut;
  }

  if (!nameStr.compare("pairJpsi")) {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    return cut;
  }

  if (!nameStr.compare("pairPsi2S")) {
    cut->AddCut(GetAnalysisCut("pairPsi2S"));
    return cut;
  }

  if (!nameStr.compare("pairUpsilon")) {
    cut->AddCut(GetAnalysisCut("pairUpsilon"));
    return cut;
  }

  if (!nameStr.compare("pairRapidityForward")) {
    cut->AddCut(GetAnalysisCut("pairRapidityForward"));
    return cut;
  }

  if (!nameStr.compare("pairJpsiLowPt1")) {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    cut->AddCut(GetAnalysisCut("pairPtLow1"));
    return cut;
  }

  if (!nameStr.compare("pairJpsiLowPt2")) {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    cut->AddCut(GetAnalysisCut("pairPtLow2"));
    return cut;
  }

  delete cut;
  return nullptr;
}

AnalysisCut* o2::aod::dqcuts::GetAnalysisCut(const char* cutName)
{
  //
  // define here cuts which are likely to be used often
  //
  AnalysisCut* cut = new AnalysisCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("eventStandard")) {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("eventStandardNoINT7")) {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    return cut;
  }

  if (!nameStr.compare("eventDimuonStandard")) {
    cut->AddCut(VarManager::kIsMuonUnlikeLowPt7, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("eventMuonStandard")) {
    cut->AddCut(VarManager::kIsMuonSingleLowPt7, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("int7vtxZ5")) {
    cut->AddCut(VarManager::kVtxZ, -5.0, 5.0);
    cut->AddCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("jpsiStandardKine")) {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("electronDalitzKine")) {
    cut->AddCut(VarManager::kPt, 0.2, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("highPtHadron")) {
    cut->AddCut(VarManager::kPt, 4.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  if (!nameStr.compare("lmeeStandardKine")) {
    cut->AddCut(VarManager::kPt, 0.2, 10.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("lmeeLowBKine")) {
    cut->AddCut(VarManager::kPt, 0.075, 10.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("PIDStandardKine")) {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("TightGlobalTrack")) {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (!nameStr.compare("TightGlobalTrackRun3")) {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (!nameStr.compare("TightTPCTrackRun3")) {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    return cut;
  }

  if (!nameStr.compare("electronStandardQuality")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (!nameStr.compare("electronStandardQualityBenchmark")) {
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  if (!nameStr.compare("SPDfirst")) {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("SPDany")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    return cut;
  }

  if (!nameStr.compare("electronStandardQualityForO2MCdebug")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (!nameStr.compare("electronStandardQualityForO2MCdebug2")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (!nameStr.compare("standardPrimaryTrack")) {
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  TF1* cutLow1 = new TF1("cutLow1", "pol1", 0., 10.);
  if (!nameStr.compare("electronPID1")) {
    cutLow1->SetParameters(130., -40.0);
    cut->AddCut(VarManager::kTPCsignal, 70., 100.);
    cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID1shiftUp")) {
    cut->AddCut(VarManager::kTPCsignal, 70. - 0.85, 100. - 0.85);
    cutLow1->SetParameters(130. - 0.85, -40.0);
    cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0 - 0.85, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID1shiftDown")) {
    cut->AddCut(VarManager::kTPCsignal, 70.0 + 0.85, 100.0 + 0.85);
    cutLow1->SetParameters(130. + 0.85, -40.0);
    cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0 + 0.85, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID1randomized")) {
    cutLow1->SetParameters(130., -40.0);
    cut->AddCut(VarManager::kTPCsignalRandomized, 70., 100.);
    cut->AddCut(VarManager::kTPCsignalRandomized, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID2")) {
    cutLow1->SetParameters(130., -40.0);
    cut->AddCut(VarManager::kTPCsignal, 73., 100.);
    cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID3")) {
    cutLow1->SetParameters(130., -40.0);
    cut->AddCut(VarManager::kTPCsignal, 60., 110.);
    cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPID2randomized")) {
    cutLow1->SetParameters(130., -40.0);
    cut->AddCut(VarManager::kTPCsignalRandomized, 73., 100.);
    cut->AddCut(VarManager::kTPCsignalRandomized, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPID_TPCnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaOpen")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.0, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaVeryVeryLoose")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaVeryLoose")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaLoose")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaMedium")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.7, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.7, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDPrKaPiRej")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (!nameStr.compare("electronPIDOnly")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (!nameStr.compare("kaonPIDnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (!nameStr.compare("kaonPIDnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDnsigmaRandomized")) {
    cut->AddCut(VarManager::kTPCnSigmaElRandomized, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPrRandomized, 3.0, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPiRandomized, 3.0, 3000.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDworseRes")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0 * 0.8, 3000.0); // emulates a 20% degradation in PID resolution
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0 * 0.8, 3000.0); // proton and pion rejections are effectively relaxed by 20%
    return cut;
  }

  if (!nameStr.compare("electronPIDshift")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0 - 0.2, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0 - 0.2, 3000.0);
    return cut;
  }

  if (!nameStr.compare("tpc_pion_rejection")) {
    TF1* f1maxPi = new TF1("f1maxPi", "[0]+[1]*x", 0, 10);
    f1maxPi->SetParameters(85, -50);
    cut->AddCut(VarManager::kTPCsignal, 70, f1maxPi, true, VarManager::kPin, 0.0, 0.4, false);
    return cut;
  }

  if (!nameStr.compare("tpc_pion_band_rejection")) {
    TF1* f1minPi = new TF1("f1minPi", "[0]+[1]*log(x)", 0, 10);
    f1minPi->SetParameters(-115, -90);
    TF1* f1maxPi = new TF1("f1maxPi", "[0]+[1]*log(x)", 0, 10);
    f1maxPi->SetParameters(-70, -90);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.05, 0.3, false);
    return cut;
  }

  if (!nameStr.compare("tpc_pion_muon_band_rejection")) {
    TF1* f1minPi = new TF1("f1minPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPi->SetParameters(37, -18, 5.5);
    TF1* f1maxPi = new TF1("f1maxPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPi->SetParameters(67, -17, 5.9);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (!nameStr.compare("tpc_pion_rejection_highp")) {
    TF1* f1minPi = new TF1("f1minPi", "[0]+[1]*x", 0, 10);
    f1minPi->SetParameters(65, 4.);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, 110., false, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (!nameStr.compare("tpc_kaon_rejection")) {
    TF1* f1minKa = new TF1("f1minKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1minKa->SetParameters(37, -4, 5.6);
    TF1* f1maxKa = new TF1("f1maxKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxKa->SetParameters(60, -4.1, 6.);
    cut->AddCut(VarManager::kTPCsignal, f1minKa, f1maxKa, true, VarManager::kPin, 0.0, 10.0, false);
    return cut;
  }

  if (!nameStr.compare("tpc_proton_rejection")) {
    TF1* f1minPr = new TF1("f1minPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPr->SetParameters(37, -2.6, 6.1);
    TF1* f1maxPr = new TF1("f1maxPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPr->SetParameters(60, -2.4, 6.2);
    cut->AddCut(VarManager::kTPCsignal, f1minPr, f1maxPr, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (!nameStr.compare("tpc_electron")) {
    cut->AddCut(VarManager::kTPCsignal, 60, 110, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("tof_electron")) {
    cut->AddCut(VarManager::kTOFbeta, 0.99, 1.01, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("tof_electron_loose")) {
    cut->AddCut(VarManager::kTOFbeta, 0.95, 1.05, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("pidcalib_ele")) {
    cut->AddCut(VarManager::kIsLegFromGamma, 0.5, 1.5, false);
    return cut;
  }

  if (!nameStr.compare("electronPID_TOFnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.4, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("pidcalib_pion")) {
    cut->AddCut(VarManager::kIsLegFromK0S, 0.5, 1.5, false);
    return cut;
  }

  if (!nameStr.compare("pidcalib_proton")) {
    cut->AddCut(VarManager::kIsProtonFromLambdaAndAntiLambda, 0.5, 1.5, false);
    return cut;
  }

  if (!nameStr.compare("pidbasic")) {
    cut->AddCut(VarManager::kTPCnclsCR, 70, 161);
    cut->AddCut(VarManager::kTPCchi2, 0, 4);
    return cut;
  }

  if (!nameStr.compare("muonQualityCuts")) {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (!nameStr.compare("muonQualityCutsMatchingOnly")) {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (!nameStr.compare("muonLowPt")) {
    cut->AddCut(VarManager::kPt, 0.5, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonLowPt2")) {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonHighPt")) {
    cut->AddCut(VarManager::kPt, 4.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonTightQualityCutsForTests")) {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 20.0, 60.0);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    return cut;
  }

  if (!nameStr.compare("mchTrack")) {
    cut->AddCut(VarManager::kMuonTrackType, 3.5, 4.5);
    return cut;
  }

  if (!nameStr.compare("matchedMchMid")) {
    cut->AddCut(VarManager::kMuonTrackType, 2.5, 3.5);
    return cut;
  }

  if (!nameStr.compare("matchedFwd")) {
    cut->AddCut(VarManager::kMuonTrackType, 1.5, 2.5);
    return cut;
  }

  if (!nameStr.compare("matchedGlobal")) {
    cut->AddCut(VarManager::kMuonTrackType, -0.5, 0.5);
    return cut;
  }

  if (!nameStr.compare("pairDalitz1")) {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    TF1* fcutHigh = new TF1("f1", "[0] - [0]/[1]*x", -1.5, 1.5);
    fcutHigh->SetParameters(0.6, 0.12);
    TF1* fcutLow = new TF1("f2", "-[0] + [0]/[1]*x", -1.5, 1.5);
    fcutLow->SetParameters(0.6, 0.12);
    cut->AddCut(VarManager::kPsiPair, fcutLow, fcutHigh, true, VarManager::kDeltaPhiPair, 0, 0.12);
    return cut;
  }

  if (!nameStr.compare("pairDalitz2")) {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    return cut;
  }

  if (!nameStr.compare("pairDalitz3")) {
    cut->AddCut(VarManager::kMass, 0.0, 0.15);
    return cut;
  }

  if (!nameStr.compare("pairNoCut")) {
    cut->AddCut(VarManager::kMass, 0.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow1")) {
    cut->AddCut(VarManager::kMass, 2.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow2")) {
    cut->AddCut(VarManager::kMass, 2.2, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow3")) {
    cut->AddCut(VarManager::kMass, 2.5, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairJpsi")) {
    cut->AddCut(VarManager::kMass, 2.8, 3.3);
    return cut;
  }

  if (!nameStr.compare("pairPsi2S")) {
    cut->AddCut(VarManager::kMass, 3.4, 3.9);
    return cut;
  }

  if (!nameStr.compare("pairUpsilon")) {
    cut->AddCut(VarManager::kMass, 8.0, 11.0);
    return cut;
  }

  if (!nameStr.compare("pairPtLow1")) {
    cut->AddCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairPtLow2")) {
    cut->AddCut(VarManager::kPt, 5.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairRapidityForward")) {
    cut->AddCut(VarManager::kRap, 2.5, 4.0);
    return cut;
  }

  if (!nameStr.compare("pairDCA")) {
    cut->AddCut(VarManager::kQuadDCAabsXY, .0, .50);
    return cut;
  }

  if (!nameStr.compare("singleDCA")) {
    cut->AddCut(VarManager::kTrackDCAsigXY, 0.0, 5.);
    return cut;
  }

  delete cut;
  return nullptr;
}
