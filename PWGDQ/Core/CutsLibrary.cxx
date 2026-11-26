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
/// \author Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
/// \brief Library with cuts to use
/// \file CutsLibrary.cxx

#include "PWGDQ/Core/CutsLibrary.h"

#include "AnalysisCompositeCut.h"
#include "AnalysisCut.h"
#include "VarManager.h"

#include <Framework/Array2D.h>
#include <Framework/Logger.h>

#include <TF1.h>
#include <TString.h>

#include <rapidjson/document.h>
#include <rapidjson/error/error.h>

#include <cstddef>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

AnalysisCompositeCut* o2::aod::dqcuts::getCompositeCut(const char* cutName)
{
  //
  // define composie cuts, typically combinations of all the ingredients needed for a full cut
  //
  // TODO: Agree on some conventions for the naming
  //       Think of possible customization of the predefined cuts via names

  auto* cut = new AnalysisCompositeCut(cutName, cutName);
  std::string nameStr = cutName;

  // ///////////////////////////////////////////////
  //   These are the Cuts used in the CEFP Task   //
  //   to select tracks in the event selection    //
  //                                              //
  //    see CutsLubrary.h for the description     //
  // ///////////////////////////////////////////////
  if (nameStr == "Electron2022") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5"));
    return cut;
  }
  if (nameStr == "Electron2023") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5_noCorr"));
    return cut;
  }
  if (nameStr == "Electron2025_1") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -3.0, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_2") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -3.0, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_3") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_4") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_5") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.25, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_4_ldong") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kP, 1.0, 1000.0);
    kineCut->addCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kITSchi2, 0.0, 5.0);
    qualityCuts->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->addCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->addCut(VarManager::kTrackDCAz, -0.5, 0.5);
    qualityCuts->addCut(VarManager::kTrackDCAxy, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.7, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPi, 2.7, 999, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->addCut(VarManager::kTPCnSigmaPr, 2.7, 999, false, VarManager::kPin, 5.0, 1000.0);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(pidCuts);
    return cut;
  }

  if (nameStr == "LowMassElectron2023") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("lmee_pp_502TeV_TOFloose_pionrej"));
    return cut;
  }
  if (nameStr == "MuonLow2022") {
    cut->addCut(getAnalysisCut("muonLowPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }
  if (nameStr == "MuonHigh2022") {
    cut->addCut(getAnalysisCut("muonHighPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }
  if (nameStr == "MuonLow2023") {
    cut->addCut(getAnalysisCut("muonLowPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }
  if (nameStr == "MuonHigh2023") {
    cut->addCut(getAnalysisCut("muonHighPt6"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }
  if (nameStr == "ElectronForEMu") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaLoose"));
    return cut;
  }
  if (nameStr == "MuonForEMu") {
    cut->addCut(getAnalysisCut("muonLowPt5"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }
  // ///////////////////////////////////////////////
  //           End of Cuts for CEFP               //
  // ///////////////////////////////////////////////

  if (nameStr == "jpsiO2MCdebugCuts") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPID1"));
    return cut;
  }

  if (nameStr == "jpsiBenchmarkCuts") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityBenchmark"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2_prefiltered1") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigma"));
    cut->addCut(getAnalysisCut("notDalitzLeg1"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts3") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts4") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaLoose"));
    return cut;
  }

  if (nameStr == "electronSelection1_ionut") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }
  if (nameStr == "jpsiKineDcaQualitynoPID") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    return cut;
  }
  if (nameStr == "jpsiKineDcaQualitynoPID2") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    return cut;
  }
  if (nameStr == "electronSelection1_ionut_withTOFPID") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium_withLargeTOFPID"));
    return cut;
  }
  if (nameStr == "electronSelection1_idstoreh") { // same as electronSelection1_ionut, but with kIsSPDAny -> kIsITSibAny
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "electronSelection1pos_ionut") {
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }
  if (nameStr == "electronSelection1neg_ionut") {
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "electronSelection2_ionut") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    cut->addCut(getAnalysisCut("insideTPCsector"));
    return cut;
  }
  if (nameStr == "electronSelection2pos_ionut") {
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    cut->addCut(getAnalysisCut("insideTPCsector"));
    return cut;
  }
  if (nameStr == "electronSelection2neg_ionut") {
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaMedium"));
    cut->addCut(getAnalysisCut("insideTPCsector"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts4_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts5") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaVeryLoose"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts6") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaVeryVeryLoose"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaOpen"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_noCorr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5_noCorr"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr_2") {
    cut->addCut(getAnalysisCut("jpsiStandardKine2"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr_3") {
    cut->addCut(getAnalysisCut("jpsiStandardKine3"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts8") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts9") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPID1shiftDown"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts10_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts10_Corr_Amb") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    cut->addCut(getAnalysisCut("ambiguousTrack")); // IsAmbiguous
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts11_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug3")); // cut on 1 ITS cluster
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts12") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->addCut(getAnalysisCut("electronPIDnsigmaVeryLoose"));     // with 3 sigma El TOF
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine4"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("pidCut_lowP_Corr"));

    auto* pidCutHighP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", false);
    pidCutHighP->addCut(getAnalysisCut("EleInclusion_highP_Corr"));
    pidCutHighP->addCut(getAnalysisCut("PionExclusion_highP_Corr"));
    cut->addCut(pidCutHighP);
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent") {
    cut->addCut(getAnalysisCut("jpsiStandardKine4"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("pidCut_lowP"));

    auto* pidCutHighP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", false);
    pidCutHighP->addCut(getAnalysisCut("EleInclusion_highP"));
    pidCutHighP->addCut(getAnalysisCut("PionExclusion_highP"));
    cut->addCut(pidCutHighP);
    return cut;
  }
  if (nameStr == "jpsiO2MCdebugCuts_Pdependent2_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine4"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("pidCut_lowP_Corr"));

    auto* pidCutHighP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", false);
    pidCutHighP->addCut(getAnalysisCut("EleInclusion_highP2_Corr"));
    pidCutHighP->addCut(getAnalysisCut("PionExclusion_highP_Corr"));
    cut->addCut(pidCutHighP);
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent2") {
    cut->addCut(getAnalysisCut("jpsiStandardKine4"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("pidCut_lowP"));

    auto* pidCutHighP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", false);
    pidCutHighP->addCut(getAnalysisCut("EleInclusion_highP2"));
    pidCutHighP->addCut(getAnalysisCut("PionExclusion_highP"));
    cut->addCut(pidCutHighP);
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts1") { // please do not remove or modify, this is used for the common Skimmed tree production, (Xiaozhi Bai)
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("electronTrackQualitySkimmed"));
    cut->addCut(getAnalysisCut("electronPIDLooseSkimmed"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts2") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("electronTrackQualitySkimmed"));
    cut->addCut(getAnalysisCut("electronPIDLooseSkimmed2"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts3") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("electronTrackQualitySkimmed2"));
    cut->addCut(getAnalysisCut("electronPIDLooseSkimmed2"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts4") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("electronTrackQualitySkimmed2"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug9")); // loose cut
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts5") {
    cut->addCut(getAnalysisCut("electronTrackQualitySkimmed3"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug8"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut") {
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    cut->addCut(getAnalysisCut("jpsiStandardKine3"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut_posEta") {
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    cut->addCut(getAnalysisCut("jpsiPIDcalibKine_posEta"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut_negEta") {
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    cut->addCut(getAnalysisCut("jpsiPIDcalibKine_negEta"));
    return cut;
  }

  if (nameStr == "pidPion_ionut") {
    cut->addCut(getAnalysisCut("pidcalib_pion"));
    cut->addCut(getAnalysisCut("jpsiStandardKine3"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts13_Corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCA")); // with DCA cut
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts14") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts14andDCA") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed"));
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "emu_electronCuts") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tof") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed"));
    cut->addCut(getAnalysisCut("tof_electron_sigma_2"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tightTPC") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed_2"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tof_tightTPC") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaSkewed_2"));
    cut->addCut(getAnalysisCut("tof_electron_sigma_2"));
    return cut;
  }

  if (nameStr == "jpsiKineAndQuality") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "jpsiPID1") {
    cut->addCut(getAnalysisCut("jpsiStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPID1"));
    return cut;
  }

  if (nameStr == "jpsiPID2") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPID2"));
    return cut;
  }

  if (nameStr == "jpsiPIDnsigma") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDCut1") {
    cut->addCut(getAnalysisCut("pionQualityCut1"));
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDCut2") {
    cut->addCut(getAnalysisCut("pionQualityCut2"));
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "PIDCalibElectron") {
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "PIDCalibPion") {
    cut->addCut(getAnalysisCut("pidcalib_pion"));
    return cut;
  }

  if (nameStr == "PIDCalibKaon") {
    cut->addCut(getAnalysisCut("pidcalib_kaon"));
    return cut;
  }

  if (nameStr == "PIDCalibProton") {
    cut->addCut(getAnalysisCut("pidcalib_proton"));
    return cut;
  }

  if (nameStr == "PIDCalib_basic") {
    cut->addCut(getAnalysisCut("pidbasic"));
    return cut;
  }

  if (nameStr == "PIDefficiency_wPID") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "PIDefficiency_woPID") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "highPtHadron") {
    cut->addCut(getAnalysisCut("highPtHadron"));
    return cut;
  }

  if (nameStr == "rho0Cuts") {
    cut->addCut(getAnalysisCut("rho0Kine"));
    cut->addCut(getAnalysisCut("pionQuality"));
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "rho0Kine") {
    cut->addCut(getAnalysisCut("rho0Kine"));
    return cut;
  }

  if (nameStr == "openEtaSel") {
    cut->addCut(getAnalysisCut("openEtaSel"));
    return cut;
  }

  if (nameStr == "hasTOF") {
    cut->addCut(getAnalysisCut("hasTOF"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts1") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("SPDany"));
    cut->addCut(getAnalysisCut("openEtaSel"));
    cut->addCut(getAnalysisCut("pionQuality"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts2") {
    cut->addCut(getAnalysisCut("muonLowPt3"));
    cut->addCut(getAnalysisCut("ITSiball"));
    cut->addCut(getAnalysisCut("openEtaSel"));
    cut->addCut(getAnalysisCut("pionQuality"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts3") {
    cut->addCut(getAnalysisCut("PIDStandardKine2"));
    cut->addCut(getAnalysisCut("SPDany"));

    auto* cutNoTPC = new AnalysisCompositeCut("NoTPC", "NoTPC", true);
    cutNoTPC->addCut(getAnalysisCut("noTPC"));

    auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
    cutTPCPID->addCut(getAnalysisCut("pionQuality"));

    auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
    cutOR->addCut(cutNoTPC);
    cutOR->addCut(cutTPCPID);
    cut->addCut(cutOR);
    return cut;
  }

  if (nameStr == "singleGapTrackCuts4") {
    cut->addCut(getAnalysisCut("PIDStandardKine2"));
    cut->addCut(getAnalysisCut("ITSibany"));

    auto* cutNoTPC = new AnalysisCompositeCut("NoTPC", "NoTPC", true);
    cutNoTPC->addCut(getAnalysisCut("noTPC"));

    auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
    cutTPCPID->addCut(getAnalysisCut("pionQuality"));

    auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
    cutOR->addCut(cutNoTPC);
    cutOR->addCut(cutTPCPID);
    cut->addCut(cutOR);
    return cut;
  }

  if (nameStr == "PIDCalib") {
    cut->addCut(getAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "NoPID") {
    cut->addCut(getAnalysisCut("PIDStandardKine2")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly2"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    return cut;
  }

  if (nameStr == "KineCutOnly") {
    cut->addCut(getAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "KineCutOnly2") {
    cut->addCut(getAnalysisCut("PIDStandardKine2")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "KineCutOnly3") {
    cut->addCut(getAnalysisCut("PIDStandardKine3")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "kaonPID") {
    cut->addCut(getAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  if (nameStr == "kaonPID2") {
    cut->addCut(getAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("kaonPIDnsigma2"));
    return cut;
  }

  if (nameStr == "kaonPID3") {
    cut->addCut(getAnalysisCut("AssocKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID3_withDCA") {        // same as kaonPID3 but with cut on DCA and SPDAny->ITSAny
    cut->addCut(getAnalysisCut("AssocKine")); // standard kine cuts usually are applied via Filter in the task
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID4") {
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID5") {
    cut->addCut(getAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  if (nameStr == "kaonPID6") {
    cut->addCut(getAnalysisCut("kaonPIDnsigma700"));
    return cut;
  }

  if (nameStr == "kaonPIDTPCTOForTPC") {
    auto* cutTPCTOFnSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", true);
    cutTPCTOFnSigma->addCut(getAnalysisCut("hasTOF"));
    cutTPCTOFnSigma->addCut(getAnalysisCut("kaonPID_TPCnTOF"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("noTOF"));
    cutTPCnSigma->addCut(getAnalysisCut("kaonPIDnsigma"));

    auto* cutPIDOR = new AnalysisCompositeCut("kaon_nsigma", "kaon_nsigma", false);
    cutPIDOR->addCut(cutTPCTOFnSigma);
    cutPIDOR->addCut(cutTPCnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "kaonPIDTPCTOForTPC700") {
    auto* cutTPCTOFnSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", true);
    cutTPCTOFnSigma->addCut(getAnalysisCut("hasTOF"));
    cutTPCTOFnSigma->addCut(getAnalysisCut("kaonPID_TPCnTOF"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("noTOF"));
    cutTPCnSigma->addCut(getAnalysisCut("kaonPIDnsigma700"));

    auto* cutPIDOR = new AnalysisCompositeCut("kaon_nsigma", "kaon_nsigma", false);
    cutPIDOR->addCut(cutTPCTOFnSigma);
    cutPIDOR->addCut(cutTPCnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "kaonPosPID4") {
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "kaonPosPID4Pt05") {
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "kaonNegPID4") {
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "kaonNegPID4Pt05") {
    cut->addCut(getAnalysisCut("kaonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "pionPID") {
    cut->addCut(getAnalysisCut("pionPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "pionPID2") {
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDTPCTOForTPC") {
    auto* cutTPCTOFnSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", true);
    cutTPCTOFnSigma->addCut(getAnalysisCut("pionPID_TPCnTOF"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("pionPIDnsigma"));

    auto* cutPIDOR = new AnalysisCompositeCut("pion_nsigma", "pion_nsigma", false);
    cutPIDOR->addCut(cutTPCTOFnSigma);
    cutPIDOR->addCut(cutTPCnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "pionPosPID") {
    cut->addCut(getAnalysisCut("pionPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "pionPosPID2") {
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "pionNegPID") {
    cut->addCut(getAnalysisCut("pionPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "pionNegPID2") {
    cut->addCut(getAnalysisCut("pionPIDnsigma"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "protonPosPID") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "protonPosPIDPt05") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "protonNegPID") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "protonNegPIDPt05") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF"));
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "protonPIDPV") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF2"));
    cut->addCut(getAnalysisCut("protonPVcut"));
    return cut;
  }

  if (nameStr == "protonPIDPV2") {
    cut->addCut(getAnalysisCut("protonPID_TPCnTOF2"));
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCAz") {
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "posPrimaryTrack_DCAz") {
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negPrimaryTrack_DCAz") {
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posStandardPrimaryTrackDCA") {
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCA"));
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negStandardPrimaryTrackDCA") {
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCA"));
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posTrack") {
    cut->addCut(getAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negTrack") {
    cut->addCut(getAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posTrackKaonRej") {
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("kaonRejNsigma"));
    return cut;
  }

  if (nameStr == "negTrackKaonRej") {
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("kaonRejNsigma"));
    return cut;
  }

  if (nameStr == "pTLow05") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "pTLow04") {
    cut->addCut(getAnalysisCut("pTLow04"));
    return cut;
  }

  if (nameStr == "pTLow03") {
    cut->addCut(getAnalysisCut("pTLow03"));
    return cut;
  }

  if (nameStr == "pTLow02") {
    cut->addCut(getAnalysisCut("pTLow02"));
    return cut;
  }

  if (nameStr == "pTLow05DCAzHigh03") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "pTLow04DCAzHigh03") {
    cut->addCut(getAnalysisCut("pTLow04"));
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "pTLow03DCAzHigh03") {
    cut->addCut(getAnalysisCut("pTLow03"));
    cut->addCut(getAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  // NOTE Below there are several TPC pid cuts used for studies of the Run3 TPC post PID calib.
  if (nameStr == "Jpsi_TPCPost_calib_debug1") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug2") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_noITSCuts_debug2") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug3") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug3"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug4") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_noITSCuts_debug4") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug6") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug2"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug6"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug7") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug2"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug7"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug8") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug5"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug8"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug9") {
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug4"));
    cut->addCut(getAnalysisCut("electronPIDLooseSkimmed3"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug10") {
    cut->addCut(getAnalysisCut("jpsiKineSkimmed"));
    cut->addCut(getAnalysisCut("jpsi_trackCut_debug6"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug10"));
    return cut;
  }

  if (nameStr == "LMee_TPCPost_calib_debug1") {
    cut->addCut(getAnalysisCut("lmee_trackCut_debug"));
    cut->addCut(getAnalysisCut("lmee_TPCPID_debug1"));
    return cut;
  }

  if (nameStr == "ITSalone_prefilter") {
    cut->addCut(getAnalysisCut("lmeePrefilterKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "ITSalonebAny_prefilter") {
    cut->addCut(getAnalysisCut("lmeePrefilterKine"));
    cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "TPCalone_prefilter") {
    cut->addCut(getAnalysisCut("lmeePrefilterKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "ITSTPC_prefilter") {
    cut->addCut(getAnalysisCut("lmeePrefilterKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  for (int iCut = 0; iCut < 10; iCut++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("jpsiEleSel%d_ionut", iCut)) {
      cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
      cut->addCut(getAnalysisCut("dcaCut1_ionut"));
      cut->addCut(getAnalysisCut("trackQuality_ionut"));
      cut->addCut(getAnalysisCut(Form("pidJpsiEle%d_ionut", iCut)));
      return cut;
    }

    if (nameStr == Form("jpsiEleSelTight%d_ionut", iCut)) {
      cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
      cut->addCut(getAnalysisCut("dcaCut1_ionut"));
      cut->addCut(getAnalysisCut("trackQualityTight_ionut"));
      cut->addCut(getAnalysisCut(Form("pidJpsiEle%d_ionut", iCut)));
      return cut;
    }
  }

  // Magnus composite cuts -----------------------------------------------------------------------------------------------------------------

  auto* magnusPID111 = new AnalysisCompositeCut("magnus_PID111", "");
  magnusPID111->addCut(getAnalysisCut("pidJpsi_magnus_ele1"));
  magnusPID111->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID111->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization111") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID111);
    return cut;
  }

  auto* magnusPID211 = new AnalysisCompositeCut("magnus_PID211", "");
  magnusPID211->addCut(getAnalysisCut("pidJpsi_magnus_ele2"));
  magnusPID211->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID211->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization211") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID211);
    return cut;
  }

  auto* magnusPID311 = new AnalysisCompositeCut("magnus_PID311", "");
  magnusPID311->addCut(getAnalysisCut("pidJpsi_magnus_ele3"));
  magnusPID311->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID311->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization311") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID311);
    return cut;
  }

  auto* magnusPID121 = new AnalysisCompositeCut("magnus_PID121", "");
  magnusPID121->addCut(getAnalysisCut("pidJpsi_magnus_ele1"));
  magnusPID121->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID121->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization121") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID121);
    return cut;
  }

  auto* magnusPID112 = new AnalysisCompositeCut("magnus_PID112", "");
  magnusPID112->addCut(getAnalysisCut("pidJpsi_magnus_ele1"));
  magnusPID112->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID112->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization112") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID112);
    return cut;
  }

  auto* magnusPID122 = new AnalysisCompositeCut("magnus_PID122", "");
  magnusPID122->addCut(getAnalysisCut("pidJpsi_magnus_ele1"));
  magnusPID122->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID122->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization122") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID122);
    return cut;
  }

  auto* magnusPID222 = new AnalysisCompositeCut("magnus_PID222", "");
  magnusPID222->addCut(getAnalysisCut("pidJpsi_magnus_ele2"));
  magnusPID222->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID222->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization222") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID222);
    return cut;
  }

  auto* magnusPID212 = new AnalysisCompositeCut("magnus_PID212", "");
  magnusPID212->addCut(getAnalysisCut("pidJpsi_magnus_ele2"));
  magnusPID212->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID212->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization212") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID212);
    return cut;
  }

  auto* magnusPID221 = new AnalysisCompositeCut("magnus_PID221", "");
  magnusPID221->addCut(getAnalysisCut("pidJpsi_magnus_ele2"));
  magnusPID221->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID221->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization221") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID221);
    return cut;
  }

  auto* magnusPID321 = new AnalysisCompositeCut("magnus_PID321", "");
  magnusPID321->addCut(getAnalysisCut("pidJpsi_magnus_ele3"));
  magnusPID321->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID321->addCut(getAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization321") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID321);
    return cut;
  }

  auto* magnusPID312 = new AnalysisCompositeCut("magnus_PID312", "");
  magnusPID312->addCut(getAnalysisCut("pidJpsi_magnus_ele3"));
  magnusPID312->addCut(getAnalysisCut("pidJpsi_magnus_pion1"));
  magnusPID312->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization312") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID312);
    return cut;
  }

  auto* magnusPID322 = new AnalysisCompositeCut("magnus_PID322", "");
  magnusPID322->addCut(getAnalysisCut("pidJpsi_magnus_ele1"));
  magnusPID322->addCut(getAnalysisCut("pidJpsi_magnus_pion2"));
  magnusPID322->addCut(getAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization322") {
    cut->addCut(getAnalysisCut("kineJpsiEle_ionut"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("trackQuality_ionut"));
    cut->addCut(magnusPID322);
    return cut;
  }
  //-------------------------------------------------------------------------------------------------------

  //---------------------------------------------------------------
  // Cuts for the selection of legs from dalitz decay
  //
  if (nameStr == "DalitzCut1") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (nameStr == "DalitzCut2") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (nameStr == "DalitzCut2_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut3") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (nameStr == "DalitzCut3_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut1SPDfirst") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (nameStr == "DalitzCut1SPDfirst_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDOnly_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut2SPDfirst") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (nameStr == "DalitzCut2SPDfirst_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut3SPDfirst") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (nameStr == "DalitzCut3SPDfirst_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  if (nameStr == "Dalitz_WithTOF_SPDfirst") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose"));

    auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
    cutTOFnSigma->addCut(getAnalysisCut("electronPID_TOFnsigma"));

    auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
    cutPIDOR->addCut(cutTPCnSigma);
    cutPIDOR->addCut(cutTOFnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "Dalitz_WithTOF_SPDfirst_Corr") {
    cut->addCut(getAnalysisCut("dalitzStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("SPDfirst"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));

    auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
    cutTOFnSigma->addCut(getAnalysisCut("electronPID_TOFnsigma_Corr"));

    auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
    cutPIDOR->addCut(cutTPCnSigma);
    cutPIDOR->addCut(cutTOFnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("dalitzSelected%d", i)) {
      cut->addCut(getAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }
  }

  if (nameStr == "electronPrimaryTag0") {
    // with tight 3 sigma DCA cut for selecting primary electrons
    cut->addCut(getAnalysisCut("electronPID_TPCnsigma_loose")); // 3 sigma inclusion, 3sigma rejection
    cut->addCut(getAnalysisCut("electronPrimary_dca3sigma"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryTag1") {
    // with 7 sigma DCA cut for selecting primary electrons
    cut->addCut(getAnalysisCut("electronPID_TPCnsigma_loose")); // 3 sigma inclusion, 3sigma rejection
    cut->addCut(getAnalysisCut("electronPrimary_dca7sigma"));
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_TPC") {
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_ITS") {
    cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCA"));
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_ITSTPC") {
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
    cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCA"));
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "jpsiPIDworseRes") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPIDworseRes"));
    return cut;
  }

  if (nameStr == "jpsiPIDshift") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPIDshift"));
    return cut;
  }

  if (nameStr == "jpsiPID1shiftUp") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (nameStr == "jpsiPID1shiftDown") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQuality"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    cut->addCut(getAnalysisCut("electronPID1shiftDown"));
    return cut;
  }
  // -------------------------------------------------------------------------------------------------
  //
  // Q vector contributor cut
  //
  if (nameStr == "selTPCCentral") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->addCut(VarManager::kEta, -0.8, 0.8);
    kineCut->addCut(VarManager::kPt, 0.15, 5);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->addCut(VarManager::kTPCchi2, 0., 4.);
    qualityCuts->addCut(VarManager::kTPCnCRoverFindCls, 0.8, 1.);
    qualityCuts->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->addCut(VarManager::kITSchi2, 0., 36.);

    auto* dcaCuts = new AnalysisCut("dcaCuts", "DCA cuts");
    std::shared_ptr<TF1> f1dcaxyHigh = std::make_shared<TF1>("f1dcaxy", "[0]+[1]/pow(x,[2])", 0., 10.);
    f1dcaxyHigh->SetParameters(0.0105, 0.035, 1.1);
    std::shared_ptr<TF1> f1dcaxyLow = std::make_shared<TF1>("f1dcaxy_low", "[0]+[1]/pow(x,[2])", 0., 10.);
    f1dcaxyLow->SetParameters(-0.0105, -0.035, 1.1);
    dcaCuts->addCut(VarManager::kTrackDCAxy, f1dcaxyLow, f1dcaxyHigh);
    dcaCuts->addCut(VarManager::kTrackDCAz, -2., 2.);

    cut->addCut(kineCut);
    cut->addCut(qualityCuts);
    cut->addCut(dcaCuts);
  }

  // -------------------------------------------------------------------------------------------------
  //
  // LMee cuts
  // List of cuts used for low mass dielectron analyses
  //
  // Skimming cuts:
  if (nameStr == "lmee_skimming") {
    cut->addCut(getAnalysisCut("lmee_skimming_cuts"));
    return cut;
  }

  // LMee Run2 PID cuts

  if (nameStr == "lmeePID_TPChadrejTOFrec") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));

    auto* cutTPCHadRej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", true);
    cutTPCHadRej->addCut(getAnalysisCut("tpc_electron"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_pion_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_kaon_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_proton_rejection"));

    auto* cutTOFRec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", true);
    cutTOFRec->addCut(getAnalysisCut("tpc_electron"));
    cutTOFRec->addCut(getAnalysisCut("tof_electron"));

    auto* cutPIDOR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", false);
    cutPIDOR->addCut(cutTPCHadRej);
    cutPIDOR->addCut(cutTOFRec);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "lmeePID_TPChadrej") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));

    auto* cutTPCHadRej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", true);
    cutTPCHadRej->addCut(getAnalysisCut("tpc_electron"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_pion_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_kaon_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_proton_rejection"));
    cut->addCut(cutTPCHadRej);
    return cut;
  }

  if (nameStr == "lmee_eNSigmaRun2") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));

    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("electronPID_TPCnsigma"));

    auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
    cutTOFnSigma->addCut(getAnalysisCut("electronPID_TOFnsigma"));

    auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
    cutPIDOR->addCut(cutTPCnSigma);
    cutPIDOR->addCut(cutTOFnSigma);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "lmeePID_TOFrec") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));

    auto* cutTOFRec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", true);
    cutTOFRec->addCut(getAnalysisCut("tpc_electron"));
    cutTOFRec->addCut(getAnalysisCut("tof_electron"));

    cut->addCut(cutTOFRec);
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun2") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun2_lowPt") {
    cut->addCut(getAnalysisCut("lmeeLowBKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun2_lowPt") {
    cut->addCut(getAnalysisCut("lmeeLowBKine"));
    cut->addCut(getAnalysisCut("TightTPCTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun2") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightTPCTrack"));
    cut->addCut(getAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  // LMee Run3 PID cuts

  if (nameStr == "lmeePID_TPChadrejTOFrecRun3") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

    auto* cutTPCHadRej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", true);
    cutTPCHadRej->addCut(getAnalysisCut("tpc_electron"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_pion_muon_band_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_pion_rejection_highp"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_kaon_rejection"));
    cutTPCHadRej->addCut(getAnalysisCut("tpc_proton_rejection"));

    auto* cutTOFRec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", true);
    cutTOFRec->addCut(getAnalysisCut("tpc_electron"));
    cutTOFRec->addCut(getAnalysisCut("tof_electron_loose"));
    cutTOFRec->addCut(getAnalysisCut("tpc_pion_rejection_highp"));

    auto* cutPIDOR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", false);
    cutPIDOR->addCut(cutTPCHadRej);
    cutPIDOR->addCut(cutTOFRec);
    cut->addCut(cutPIDOR);
    return cut;
  }

  if (nameStr == "lmee_Run3_TPCelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_TOFelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("tof_electron_sigma"));
    return cut;
  }

  if (nameStr == "lmee_Run3_posTrack_TPCelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_posTrack_TOFelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("posTrack"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("tof_electron_sigma"));
    return cut;
  }

  if (nameStr == "lmee_Run3_negTrack_TPCelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_negTrack_TOFelectron") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("negTrack"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
    cut->addCut(getAnalysisCut("tof_electron_sigma"));
    return cut;
  }

  std::vector<TString> vecTypetrack;
  vecTypetrack.emplace_back("");                  // default TightGlobalTrackRun3
  vecTypetrack.emplace_back("_7ITSncls");         // default TightGlobalTrackRun3 but with 7 ITS clusters
  vecTypetrack.emplace_back("_ITS");              // Ask only for ITS requirements
  vecTypetrack.emplace_back("_ITSalone");         // Ask only for ITS requirements + ITSalone (no TPC matching)
  vecTypetrack.emplace_back("_TPC");              // Ask only for TPC requirements
  vecTypetrack.emplace_back("_TPCalone");         // Ask only for TPC requirements + TPCalone (no ITS matching)
  vecTypetrack.emplace_back("_TPCnoTRD");         // Ask only for TPC requirements no TRD matching
  vecTypetrack.emplace_back("_TPCstrongncls");    // default TightGlobalTrackRun3 but with 130 TPC clusters
  vecTypetrack.emplace_back("_ITSanyfirsttwo");   // default TightGlobalTrackRun3 but with a cluster in any of the first two layers
  vecTypetrack.emplace_back("_ITSanyfirstthree"); // default TightGlobalTrackRun3 but with a cluster in any of the first three layers

  // loop to define PID cuts with and without post calibration
  for (size_t icase = 0; icase < vecTypetrack.size(); icase++) {
    // Tracking cuts of Pb--Pb analysis
    if (nameStr == Form("lmee%s_PbPb_selection", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee%s_PbPb_selection_pt04", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee%s_TrackCuts_Resol", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("openEtaSel")); // No pt cut and wider eta cut to produce resolution maps
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 4 cuts to separate pos & neg tracks in pos & neg eta range
    if (nameStr == Form("lmee_posTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("posEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("posEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_posTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("negEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("negEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 2 cuts to separate pos & neg tracks
    if (nameStr == Form("lmee_posTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("etaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("etaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 4 cuts to separate pos & neg tracks in pos & neg eta range low B field
    if (nameStr == Form("lmee_lowB_posTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("posEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("posEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_posTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("negEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("negEtaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    // 2 cuts to separate pos & neg tracks in low B field
    if (nameStr == Form("lmee_lowB_posTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("posTrack"));
      cut->addCut(getAnalysisCut("etaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->addCut(getAnalysisCut("negTrack"));
      cut->addCut(getAnalysisCut("etaSel"));
      cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }
  }

  std::vector<TString> vecPIDcase;
  vecPIDcase.emplace_back("");              // without post calibration
  vecPIDcase.emplace_back("_Corr");         // case of using post calibrated PID spectra
  vecPIDcase.emplace_back("_CorrWithKaon"); // case of using post calibrated PID spectra with also the kaons

  std::vector<TString> vecTypetrackWithPID;
  vecTypetrackWithPID.emplace_back("");                  // default TightGlobalTrackRun3
  vecTypetrackWithPID.emplace_back("_7ITSncls");         // default TightGlobalTrackRun3 but with 7 ITS clusters
  vecTypetrackWithPID.emplace_back("_TPCstrongncls");    // default TightGlobalTrackRun3 but with 130 TPC clusters
  vecTypetrackWithPID.emplace_back("_ITSanyfirsttwo");   // default TightGlobalTrackRun3 but with a cluster in any of the first two layers
  vecTypetrackWithPID.emplace_back("_ITSanyfirstthree"); // default TightGlobalTrackRun3 but with a cluster in any of the first three layers

  // loop to define PID cuts with and without post calibration
  for (size_t icase = 0; icase < vecPIDcase.size(); icase++) {
    if (nameStr == Form("lmee_onlyTPCPID%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITSTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
      cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      cut->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITS_ifTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutNoTPC = new AnalysisCompositeCut("NoTPC", "NoTPC", true);
      cutNoTPC->addCut(getAnalysisCut("noTPC"));

      auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
      cutTPCPID->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cutTPCPID->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
      cutOR->addCut(cutNoTPC);
      cutOR->addCut(cutTPCPID);
      cut->addCut(cutOR);
      return cut;
    }

    if (nameStr == Form("ITS_ifTPCStandard_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualityITSOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutNoTPCStandard = new AnalysisCompositeCut("NoTPCstandard", "NoTPCstandard", true);
      cutNoTPCStandard->addCut(getAnalysisCut("NoelectronStandardQualityTPCOnly"));

      auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
      cutTPCPID->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cutTPCPID->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
      cutOR->addCut(cutNoTPCStandard);
      cutOR->addCut(cutTPCPID);
      cut->addCut(cutOR);
      return cut;
    }

    if (nameStr == Form("ITSTPCbAny_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      cut->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITSbAny_ifTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutNoTPC = new AnalysisCompositeCut("NoTPC", "NoTPC", true);
      cutNoTPC->addCut(getAnalysisCut("noTPC"));

      auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
      cutTPCPID->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cutTPCPID->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
      cutOR->addCut(cutNoTPC);
      cutOR->addCut(cutTPCPID);
      cut->addCut(cutOR);
      return cut;
    }

    if (nameStr == Form("ITSbAny_ifTPCStandard_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeePrefilterKine"));
      cut->addCut(getAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutNoTPCStandard = new AnalysisCompositeCut("NoTPCstandard", "NoTPCstandard", true);
      cutNoTPCStandard->addCut(getAnalysisCut("NoelectronStandardQualityTPCOnly"));

      auto* cutTPCPID = new AnalysisCompositeCut("pid_TPC", "pid_TPC", true);
      cutTPCPID->addCut(getAnalysisCut("electronStandardQualityTPCOnly"));
      cutTPCPID->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cutOR = new AnalysisCompositeCut("OR", "OR", false);
      cutOR->addCut(cutNoTPCStandard);
      cutOR->addCut(cutTPCPID);
      cut->addCut(cutOR);
      return cut;
    }

    for (unsigned int i = 0; i < 30; i++) { // o2-linter: disable=magic-number (number of cuts)
      if (nameStr == Form("ElSelCutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeCutVarTrackCuts%i", i)));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }
    }

    for (size_t jcase = 0; jcase < vecTypetrackWithPID.size(); jcase++) {
      // All previous cut with TightGlobalTrackRun3
      if (nameStr == Form("ITSTPC%s_TPCPIDalone%s_PbPb", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_loose", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_loose", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_loose", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongHadRej", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TOFreq", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
        cut->addCut(getAnalysisCut(Form("electronPID_TPC_TOFnsigma%s", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_Resol", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("openEtaSel")); // No pt cut and wider eta cut to produce resolution maps
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_tightNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_tightNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_eNSigmaRun3%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeLowBKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz")); // to reject looper using DCAz

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_eNSigmaRun3%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeLowBKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz")); // to reject looper using DCAz

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3%s_tightNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      // 8 cuts for QC
      if (nameStr == Form("lmee%s_NSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt02Sel"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt02Sel"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt04Sel"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt04Sel"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      // 6 cuts for QC
      if (nameStr == Form("lmee%s_posTOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_posTOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt04Sel"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("pt04Sel"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_negTOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_negTOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));
        cut->addCut(getAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TPC_PID", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        cut->addCut(cutTPCnSigma);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TOF_PID", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        cut->addCut(cutTOFnSigma);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE_DCA05", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_DCA05"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID
      if (nameStr == Form("lmee%s_posNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID for low B field
      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID for low B field with bad TOF rejection
      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_posEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_posEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("posEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_negEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("posTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_negEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("negTrack"));
        cut->addCut(getAnalysisCut("negEtaSel"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      // some older cuts
      if (nameStr == Form("lmee%s_pp502TeV_PID%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
        if (nameStr == Form("lmee%s_pp502TeV_PID%s_UsePrefilter%d", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data(), i)) {
          cut->addCut(getAnalysisCut(Form("notDalitzLeg%d", i)));
          cut->addCut(getAnalysisCut("lmeeStandardKine"));
          cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
          cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

          auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
          cutTPCnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())));

          auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
          cutTOFnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())));

          auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
          cutPIDOR->addCut(cutTPCnSigma);
          cutPIDOR->addCut(cutTOFnSigma);
          cut->addCut(cutPIDOR);
          return cut;
        }
      }

      if (nameStr == Form("lmee%s_pp502TeV_lowB_PID%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeLowBKine"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("standardPrimaryTrackDCAz")); // DCAz to reject loopers

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_lowB_TPC%s", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_lowB_TOF%s", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->addCut(getAnalysisCut("lmeeStandardKine_pt04"));
        cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
        cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
        cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
        cutPIDOR->addCut(cutTPCnSigma);
        cutPIDOR->addCut(cutTOFnSigma);
        cut->addCut(cutPIDOR);
        return cut;
      }

      for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
        if (nameStr == Form("lmee%s_eNSigmaRun3%s_UsePrefilter%d", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data(), i)) {
          cut->addCut(getAnalysisCut(Form("notDalitzLeg%d", i)));
          cut->addCut(getAnalysisCut("lmeeStandardKine"));
          cut->addCut(getAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
          cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

          auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
          cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

          auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
          cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

          auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
          cutPIDOR->addCut(cutTPCnSigma);
          cutPIDOR->addCut(cutTOFnSigma);
          cut->addCut(cutPIDOR);
          return cut;
        }
      }
    }

    if (nameStr == Form("lmee_eNSigmaRun3%s_strongTPC", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("TightGlobalTrackRun3_strongTPC"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
      cutTPCnSigma->addCut(getAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

      auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
      cutTOFnSigma->addCut(getAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

      auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
      cutPIDOR->addCut(cutTPCnSigma);
      cutPIDOR->addCut(cutTOFnSigma);
      cut->addCut(cutPIDOR);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_PID%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
      cutTPCnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TPCloose%s", vecPIDcase.at(icase).Data())));

      auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
      cutTOFnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())));

      auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
      cutPIDOR->addCut(cutTPCnSigma);
      cutPIDOR->addCut(cutTOFnSigma);
      cut->addCut(cutPIDOR);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtestb_PID%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
      cutTPCnSigma->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TPCloosenopkrej%s", vecPIDcase.at(icase).Data())));

      cut->addCut(cutTPCnSigma);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      cut->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF_pionrej%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
      cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
      cut->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF_pionrej_noDCA%s", vecPIDcase.at(icase).Data())) {
      cut->addCut(getAnalysisCut("lmeeStandardKine"));
      cut->addCut(getAnalysisCut("LooseGlobalTrackRun3"));
      cut->addCut(getAnalysisCut(Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())));
      return cut;
    }
  }

  if (nameStr == "testCut_chic") {
    cut->addCut(getAnalysisCut("jpsiStandardKine5"));
    cut->addCut(getAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->addCut(getAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun3") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun3_lowPt") {
    cut->addCut(getAnalysisCut("lmeeLowBKine"));
    cut->addCut(getAnalysisCut("TightGlobalTrackRun3"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun3_lowPt") {
    cut->addCut(getAnalysisCut("lmeeLowBKine"));
    cut->addCut(getAnalysisCut("TightTPCTrackRun3"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun3") {
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("TightTPCTrackRun3"));
    cut->addCut(getAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "trackCut_compareDQEMframework") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->addCut(getAnalysisCut("lmeeStandardKine"));
    cut->addCut(getAnalysisCut("trackQuality_compareDQEMframework"));
    cut->addCut(getAnalysisCut("trackDCA1cm"));
    auto* cutTPCnSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", true);
    cutTPCnSigma->addCut(getAnalysisCut("lmee_commonDQEM_PID_TPC"));

    auto* cutTOFnSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", true);
    cutTOFnSigma->addCut(getAnalysisCut("lmee_commonDQEM_PID_TOF"));

    auto* cutPIDOR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", false);
    cutPIDOR->addCut(cutTPCnSigma);
    cutPIDOR->addCut(cutTOFnSigma);
    return cut;
  }

  if (nameStr == "jpsi_debug_TPCTOF3_rejBadTOF") {
    cut->addCut(getAnalysisCut("jpsiStandardKine5"));
    cut->addCut(getAnalysisCut("electronStandardQualityTPCOnly3"));
    cut->addCut(getAnalysisCut("SPDfirst"));
    cut->addCut(getAnalysisCut("dcaCut1_ionut"));
    cut->addCut(getAnalysisCut("pidJpsi_TPCpion0"));
    cut->addCut(getAnalysisCut("pidJpsi_beta"));
    cut->addCut(getAnalysisCut("pidJpsi_noTOF_prot"));
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // lmee pair cuts

  if (nameStr == "pairPhiV") {
    auto* cutPairPhiV = new AnalysisCompositeCut("cut_pairPhiV", "cut_pairPhiV", true);
    cutPairPhiV->addCut(getAnalysisCut("pairLowMass"));
    cutPairPhiV->addCut(getAnalysisCut("pairPhiV"));
    cut->addCut(cutPairPhiV);
    return cut;
  }

  if (nameStr == "excludePairPhiV") {
    auto* cutPairLowPhiV = new AnalysisCompositeCut("cut_pairlowPhiV", "cut_pairlowPhiV", false);
    cutPairLowPhiV->addCut(getAnalysisCut("excludePairLowMass"));
    cutPairLowPhiV->addCut(getAnalysisCut("excludePairPhiV"));
    cut->addCut(cutPairLowPhiV);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts

  if (nameStr == "muonQualityCutsMatchingOnly") {
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonEtaCut") {
    cut->addCut(getAnalysisCut("muonEtaCut"));
    return cut;
  }

  if (nameStr == "muonMinimalCuts") {
    cut->addCut(getAnalysisCut("muonMinimalCuts"));
    return cut;
  }

  if (nameStr == "muonMinimalCuts10SigmaPDCA") {
    cut->addCut(getAnalysisCut("muonMinimalCuts10SigmaPDCA"));
    return cut;
  }

  if (nameStr == "muonQualityCuts") {
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "matchedQualityCuts") {
    cut->addCut(getAnalysisCut("matchedQualityCuts"));
    return cut;
  }

  if (nameStr == "matchedQualityCutsMFTeta") {
    cut->addCut(getAnalysisCut("matchedQualityCutsMFTeta"));
    return cut;
  }

  if (nameStr == "muonQualityCuts5SigmaPDCA_Run3") {
    cut->addCut(getAnalysisCut("muonQualityCuts5SigmaPDCA_Run3"));
    return cut;
  }

  if (nameStr == "muonLowPt5SigmaPDCA_Run3") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonQualityCuts5SigmaPDCA_Run3"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonQualityCuts10SigmaPDCA_MCHMID") {
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt1p5_MchMidmatching") {
    cut->addCut(getAnalysisCut("muonLowPt1p5"));
    cut->addCut(getAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (nameStr == "muonLowPt1p7_MchMidmatching") {
    cut->addCut(getAnalysisCut("muonLowPt1p7"));
    cut->addCut(getAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (nameStr == "muonLowPt6_MchMidmatching") {
    cut->addCut(getAnalysisCut("muonLowPt6"));
    cut->addCut(getAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (nameStr == "muonLowPt10SigmaPDCA") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt210SigmaPDCA") {
    cut->addCut(getAnalysisCut("muonLowPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt510SigmaPDCA") {
    cut->addCut(getAnalysisCut("muonLowPt5"));
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt610SigmaPDCA") {
    cut->addCut(getAnalysisCut("muonLowPt6"));
    cut->addCut(getAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->addCut(getAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt2") {
    cut->addCut(getAnalysisCut("muonLowPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt3") {
    cut->addCut(getAnalysisCut("muonLowPt3"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt4") {
    cut->addCut(getAnalysisCut("muonLowPt4"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt5") {
    cut->addCut(getAnalysisCut("muonLowPt5"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt6") {
    cut->addCut(getAnalysisCut("muonLowPt6"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt1p5") {
    cut->addCut(getAnalysisCut("muonLowPt1p5"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt1p7") {
    cut->addCut(getAnalysisCut("muonLowPt1p7"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly2") {
    cut->addCut(getAnalysisCut("muonLowPt2"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly3") {
    cut->addCut(getAnalysisCut("muonLowPt3"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly4") {
    cut->addCut(getAnalysisCut("muonLowPt4"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly5") {
    cut->addCut(getAnalysisCut("muonLowPt5"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly6") {
    cut->addCut(getAnalysisCut("muonLowPt6"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPt") {
    cut->addCut(getAnalysisCut("muonHighPt"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt2") {
    cut->addCut(getAnalysisCut("muonHighPt2"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt3") {
    cut->addCut(getAnalysisCut("muonHighPt3"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt4") {
    cut->addCut(getAnalysisCut("muonHighPt4"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt5") {
    cut->addCut(getAnalysisCut("muonHighPt5"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt6") {
    cut->addCut(getAnalysisCut("muonHighPt6"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly2") {
    cut->addCut(getAnalysisCut("muonHighPt2"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly3") {
    cut->addCut(getAnalysisCut("muonHighPt3"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonTightQualityCutsForTests") {
    cut->addCut(getAnalysisCut("muonTightQualityCutsForTests"));
    return cut;
  }

  if (nameStr == "mchTrack") {
    cut->addCut(getAnalysisCut("mchTrack"));
    return cut;
  }

  if (nameStr == "matchedMchMid") {
    cut->addCut(getAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (nameStr == "muonMchAndMchMid") {
    cut->addCut(getAnalysisCut("muonMchAndMchMid"));
    return cut;
  }

  if (nameStr == "muonStreamCuts") {
    cut->addCut(getAnalysisCut("muonMchAndMchMid"));
    cut->addCut(getAnalysisCut("muonMinimalCuts10SigmaPDCA"));
    return cut;
  }

  if (nameStr == "matchedFwd") {
    cut->addCut(getAnalysisCut("matchedFwd"));
    return cut;
  }

  if (nameStr == "matchedGlobal") {
    cut->addCut(getAnalysisCut("matchedGlobal"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut1") {
    cut->addCut(getAnalysisCut("Chi2MCHMFTCut1"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut2") {
    cut->addCut(getAnalysisCut("Chi2MCHMFTCut2"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut3") {
    cut->addCut(getAnalysisCut("Chi2MCHMFTCut3"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut4") {
    cut->addCut(getAnalysisCut("Chi2MCHMFTCut4"));
    return cut;
  }

  if (nameStr == "muonQualityCutsMUONStandalone") {
    cut->addCut(getAnalysisCut("matchedMchMid"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonQualityCutsGlobal") {
    cut->addCut(getAnalysisCut("matchedGlobal"));
    cut->addCut(getAnalysisCut("muonQualityCuts"));
    return cut;
  }
  // -----------------------------------------------------------
  // Pair cuts
  if (nameStr == "pairNoCut") {
    cut->addCut(getAnalysisCut("pairNoCut"));
    return cut;
  }

  if (nameStr == "pairMassLow1") {
    cut->addCut(getAnalysisCut("pairMassLow1"));
    return cut;
  }

  if (nameStr == "pairMassLow2") {
    cut->addCut(getAnalysisCut("pairMassLow2"));
    return cut;
  }

  if (nameStr == "pairPtLow3") {
    cut->addCut(getAnalysisCut("pairPtLow3"));
    return cut;
  }

  if (nameStr == "pairPtLow4") {
    cut->addCut(getAnalysisCut("pairPtLow4"));
    return cut;
  }

  if (nameStr == "pairPtLow5") {
    cut->addCut(getAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairMassLow3") {
    cut->addCut(getAnalysisCut("pairMassLow3"));
    return cut;
  }

  if (nameStr == "pairMassLow4") {
    cut->addCut(getAnalysisCut("pairMassLow4"));
    return cut;
  }

  if (nameStr == "pairMassLow5") {
    cut->addCut(getAnalysisCut("pairMassLow5"));
    return cut;
  }

  if (nameStr == "pairMassLow6") {
    cut->addCut(getAnalysisCut("pairMassLow6"));
    return cut;
  }

  if (nameStr == "pairMassLow7") {
    cut->addCut(getAnalysisCut("pairMassLow7"));
    return cut;
  }

  if (nameStr == "pairMassLow8") {
    cut->addCut(getAnalysisCut("pairMassLow8"));
    return cut;
  }

  if (nameStr == "pairMassLow9") {
    cut->addCut(getAnalysisCut("pairMassLow9"));
    return cut;
  }

  if (nameStr == "pairMassLow10") {
    cut->addCut(getAnalysisCut("pairMassLow10"));
    return cut;
  }

  if (nameStr == "pairMassLow11") {
    cut->addCut(getAnalysisCut("pairMassLow11"));
    return cut;
  }

  if (nameStr == "pairMassLow12") {
    cut->addCut(getAnalysisCut("pairMassLow12"));
    return cut;
  }

  if (nameStr == "pairMass1to2") {
    cut->addCut(getAnalysisCut("pairMass1to2"));
    return cut;
  }

  if (nameStr == "pairMassIMR") {
    cut->addCut(getAnalysisCut("pairMassIMR"));
    return cut;
  }

  if (nameStr == "pairMass1_5to2_7") {
    cut->addCut(getAnalysisCut("pairMass1_5to2_7"));
    return cut;
  }

  if (nameStr == "pairMass1_3to3_5") {
    cut->addCut(getAnalysisCut("pairMass1_3to3_5"));
    return cut;
  }

  if (nameStr == "pairMass1_3") {
    cut->addCut(getAnalysisCut("pairMass1_3"));
    return cut;
  }

  if (nameStr == "pairMass1_5to3_5") {
    cut->addCut(getAnalysisCut("pairMass1_5to3_5"));
    return cut;
  }

  if (nameStr == "pairDalitz1") {
    cut->addCut(getAnalysisCut("pairDalitz1"));
    return cut;
  }

  if (nameStr == "pairDalitz1Strong") {
    cut->addCut(getAnalysisCut("pairDalitz1Strong"));
    return cut;
  }

  if (nameStr == "pairDalitz2") {
    cut->addCut(getAnalysisCut("pairDalitz2"));
    return cut;
  }

  if (nameStr == "pairDalitz3") {
    cut->addCut(getAnalysisCut("pairDalitz3"));
    return cut;
  }

  if (nameStr == "paira_prefilter1") {
    cut->addCut(getAnalysisCut("paira_prefilter1"));
    return cut;
  }

  if (nameStr == "paira_prefilter2") {
    cut->addCut(getAnalysisCut("paira_prefilter2"));
    return cut;
  }

  if (nameStr == "paira_prefilter3") {
    cut->addCut(getAnalysisCut("paira_prefilter3"));
    return cut;
  }

  if (nameStr == "paira_prefilter4") {
    cut->addCut(getAnalysisCut("paira_prefilter4"));
    return cut;
  }

  if (nameStr == "paira_prefilter5") {
    cut->addCut(getAnalysisCut("paira_prefilter5"));
    return cut;
  }

  if (nameStr == "paira_prefilter6") {
    cut->addCut(getAnalysisCut("paira_prefilter6"));
    return cut;
  }

  if (nameStr == "paira_prefilter7") {
    cut->addCut(getAnalysisCut("paira_prefilter7"));
    return cut;
  }

  if (nameStr == "pairb_prefilter1") {
    cut->addCut(getAnalysisCut("pairb_prefilter1"));
    return cut;
  }

  if (nameStr == "pairb_prefilter2") {
    cut->addCut(getAnalysisCut("pairb_prefilter2"));
    return cut;
  }

  if (nameStr == "pairb_prefilter3") {
    cut->addCut(getAnalysisCut("pairb_prefilter3"));
    return cut;
  }

  if (nameStr == "pairb_prefilter4") {
    cut->addCut(getAnalysisCut("pairb_prefilter4"));
    return cut;
  }

  if (nameStr == "pairb_prefilter5") {
    cut->addCut(getAnalysisCut("pairb_prefilter5"));
    return cut;
  }

  if (nameStr == "pairb_prefilter6") {
    cut->addCut(getAnalysisCut("pairb_prefilter6"));
    return cut;
  }

  if (nameStr == "pairb_prefilter7") {
    cut->addCut(getAnalysisCut("pairb_prefilter7"));
    return cut;
  }

  if (nameStr == "pairc_prefilter1") {
    cut->addCut(getAnalysisCut("pairc_prefilter1"));
    return cut;
  }

  if (nameStr == "pairc_prefilter2") {
    cut->addCut(getAnalysisCut("pairc_prefilter2"));
    return cut;
  }

  if (nameStr == "pairc_prefilter3") {
    cut->addCut(getAnalysisCut("pairc_prefilter3"));
    return cut;
  }

  if (nameStr == "pairc_prefilter4") {
    cut->addCut(getAnalysisCut("pairc_prefilter4"));
    return cut;
  }

  if (nameStr == "pairc_prefilter5") {
    cut->addCut(getAnalysisCut("pairc_prefilter5"));
    return cut;
  }

  if (nameStr == "pairc_prefilter6") {
    cut->addCut(getAnalysisCut("pairc_prefilter6"));
    return cut;
  }

  if (nameStr == "pairc_prefilter7") {
    cut->addCut(getAnalysisCut("pairc_prefilter7"));
    return cut;
  }

  if (nameStr == "paird_prefilter1") {
    cut->addCut(getAnalysisCut("paird_prefilter1"));
    return cut;
  }

  if (nameStr == "paire_prefilter1") {
    cut->addCut(getAnalysisCut("paire_prefilter1"));
    return cut;
  }

  if (nameStr == "pairf_prefilter1") {
    cut->addCut(getAnalysisCut("pairf_prefilter1"));
    return cut;
  }

  if (nameStr == "pairg_prefilter1") {
    cut->addCut(getAnalysisCut("pairg_prefilter1"));
    return cut;
  }

  if (nameStr == "pairh_prefilter1") {
    cut->addCut(getAnalysisCut("pairh_prefilter1"));
    return cut;
  }

  if (nameStr == "pairi_prefilter1") {
    cut->addCut(getAnalysisCut("pairi_prefilter1"));
    return cut;
  }

  if (nameStr == "pairJpsi") {
    cut->addCut(getAnalysisCut("pairJpsi"));
    return cut;
  }

  if (nameStr == "pairJpsi2") {
    cut->addCut(getAnalysisCut("pairJpsi2"));
    return cut;
  }

  if (nameStr == "pairJpsi3") {
    cut->addCut(getAnalysisCut("pairJpsi3"));
    return cut;
  }

  if (nameStr == "pairPsi2S") {
    cut->addCut(getAnalysisCut("pairPsi2S"));
    return cut;
  }

  if (nameStr == "pairUpsilon") {
    cut->addCut(getAnalysisCut("pairUpsilon"));
    return cut;
  }

  if (nameStr == "pairX3872Cut1") {
    cut->addCut(getAnalysisCut("pairX3872"));
    return cut;
  }

  if (nameStr == "pairX3872Cut2") {
    cut->addCut(getAnalysisCut("pairX3872_2"));
    return cut;
  }

  if (nameStr == "pairX3872Cut3") {
    cut->addCut(getAnalysisCut("pairX3872_3"));
    return cut;
  }

  if (nameStr == "pairX3872Minitree") {
    cut->addCut(getAnalysisCut("pairX3872_minitree"));
    return cut;
  }

  if (nameStr == "DipionPairCut1") {
    cut->addCut(getAnalysisCut("DipionMassCut1"));
    return cut;
  }

  if (nameStr == "DipionPairCut2") {
    cut->addCut(getAnalysisCut("DipionMassCut2"));
    return cut;
  }

  if (nameStr == "pairRapidityForward") {
    cut->addCut(getAnalysisCut("pairRapidityForward"));
    return cut;
  }

  if (nameStr == "pairJpsiLowPt1") {
    cut->addCut(getAnalysisCut("pairJpsi"));
    cut->addCut(getAnalysisCut("pairPtLow1"));
    return cut;
  }

  if (nameStr == "pairJpsiLowPt2") {
    cut->addCut(getAnalysisCut("pairJpsi"));
    cut->addCut(getAnalysisCut("pairPtLow2"));
    return cut;
  }

  if (nameStr == "pairCoherentRho0") {
    cut->addCut(getAnalysisCut("pairPtLow3"));
    return cut;
  }

  if (nameStr == "pairD0") {
    cut->addCut(getAnalysisCut("pairD0"));
    return cut;
  }

  if (nameStr == "pairD0Topology") {
    cut->addCut(getAnalysisCut("pairD0Topology"));
    return cut;
  }

  if (nameStr == "pairD0TopologyTight") {
    cut->addCut(getAnalysisCut("pairD0TopologyTight"));
    return cut;
  }

  if (nameStr == "pairLambdaCTopology") {
    cut->addCut(getAnalysisCut("pairLambdaCTopology"));
    return cut;
  }

  if (nameStr == "pairD0HighPt1") {
    cut->addCut(getAnalysisCut("pairLxyzProjected3sigma"));
    cut->addCut(getAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairD0HighPt2") {
    cut->addCut(getAnalysisCut("pairTauxyzProjected1"));
    cut->addCut(getAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairD0HighPt3") {
    cut->addCut(getAnalysisCut("pairTauxyzProjected1sigma"));
    cut->addCut(getAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1") {
    cut->addCut(getAnalysisCut("pairTauxyzProjected1"));
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaLambdacCand") {
    cut->addCut(getAnalysisCut("pairLxyProjected3sigmaLambdacCand"));
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaDplusCand") {
    cut->addCut(getAnalysisCut("pairLxyProjected3sigmaDplusCand"));
    return cut;
  }

  if (nameStr == "pairCosPointingPos") {
    cut->addCut(getAnalysisCut("pairCosPointingPos"));
    return cut;
  }

  if (nameStr == "pairCosPointingNeg90") {
    cut->addCut(getAnalysisCut("pairCosPointingNeg90"));
    return cut;
  }

  if (nameStr == "pairCosPointingNeg85") {
    cut->addCut(getAnalysisCut("pairCosPointingNeg85"));
    return cut;
  }

  if (nameStr == "pairTauxyzProjectedCosPointing1") {
    cut->addCut(getAnalysisCut("pairCosPointingNeg"));
    cut->addCut(getAnalysisCut("pairTauxyzProjected1"));
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and in order or optimize the trigger
  // trigger selection cuts

  if (nameStr == "jpsiO2TriggerTestCuts_LooseNsigma") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_LooseNsigma_corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug5"));
    return cut;
  }
  if (nameStr == "jpsiO2TriggerTestCuts_MediumNsigma") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_MediumNsigma_corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (nameStr == "jpsiO2TriggerTestCuts_TightNsigma") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TightNsigma_corr") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID1") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_TriggerTest1"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID2") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_TriggerTest2"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID3") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_TriggerTest3"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID4") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronStandardQualityTriggerTest"));
    cut->addCut(getAnalysisCut("jpsi_TPCPID_TriggerTest4"));
    return cut;
  }

  if (nameStr == "emu_electron_test") {
    cut->addCut(getAnalysisCut("jpsiStandardKine"));
    cut->addCut(getAnalysisCut("electronTrackQuality_Maolin"));
    cut->addCut(getAnalysisCut("electronPIDnsigmaEMu"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts") {
    cut->addCut(getAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_LowPt") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_HighPt2") {
    cut->addCut(getAnalysisCut("muonHighPt2"));
    cut->addCut(getAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_HighPt3") {
    cut->addCut(getAnalysisCut("muonHighPt3"));
    cut->addCut(getAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly2") {
    cut->addCut(getAnalysisCut("muonHighPt2"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly3") {
    cut->addCut(getAnalysisCut("muonHighPt3"));
    cut->addCut(getAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts") {
    cut->addCut(getAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts_LowPt") {
    cut->addCut(getAnalysisCut("muonLowPt"));
    cut->addCut(getAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  //---------------------------------------------------------------
  // ALICE 3 studies composite cuts

  if (nameStr == "alice3StandardKine") {
    cut->addCut(getAnalysisCut("alice3StandardKine"));
    return cut;
  }

  if (nameStr == "alice3KineSkim") {
    cut->addCut(getAnalysisCut("alice3KineSkim"));
    return cut;
  }

  if (nameStr == "alice3TrackQuality") {
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TrackQualityTightDCA") {
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrack") {
    cut->addCut(getAnalysisCut("alice3StandardKine"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackTOFAcceptance") {
    cut->addCut(getAnalysisCut("alice3KineTOFAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackRICHAcceptance") {
    cut->addCut(getAnalysisCut("alice3KineRICHAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackD0") {
    cut->addCut(getAnalysisCut("alice3StandardKineD0"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TOFEta") {
    cut->addCut(getAnalysisCut("alice3TOFEta"));
    cut->addCut(getAnalysisCut("alice3TrackQualityHF"));
    return cut;
  }

  if (nameStr == "alice3RICHEta") {
    cut->addCut(getAnalysisCut("alice3RICHEta"));
    cut->addCut(getAnalysisCut("alice3TrackQualityHF"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackTOFAcceptanceD0") {
    cut->addCut(getAnalysisCut("alice3KineTOFAcceptanceD0"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackRICHAcceptanceD0") {
    cut->addCut(getAnalysisCut("alice3KineRICHAcceptanceD0"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrack") {
    cut->addCut(getAnalysisCut("alice3StandardKine"));
    cut->addCut(getAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrackTOFAcceptance") {
    cut->addCut(getAnalysisCut("alice3KineTOFAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrackRICHAcceptance") {
    cut->addCut(getAnalysisCut("alice3KineRICHAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDEl") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPi") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDKa") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPr") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDEl") {
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPi") {
    cut->addCut(getAnalysisCut("alice3oTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDKa") {
    cut->addCut(getAnalysisCut("alice3oTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPr") {
    cut->addCut(getAnalysisCut("alice3oTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDEl") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDPi") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPi"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDKa") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDKa"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDPr") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPr"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3EitherTOFPIDPi") {
    auto* pidCut = new AnalysisCompositeCut(
      "alice3EitherTOFPIDPi",
      "iTOF OR oTOF pion PID",
      false);

    pidCut->addCut(getAnalysisCut("alice3iTOFPIDPi"));
    pidCut->addCut(getAnalysisCut("alice3oTOFPIDPi"));

    delete cut;
    return pidCut;
  }

  if (nameStr == "alice3EitherTOFPIDKa") {
    auto* pidCut = new AnalysisCompositeCut(
      "alice3EitherTOFPIDKa",
      "iTOF OR oTOF kaon PID",
      false);

    pidCut->addCut(getAnalysisCut("alice3iTOFPIDKa"));
    pidCut->addCut(getAnalysisCut("alice3oTOFPIDKa"));

    delete cut;
    return pidCut;
  }

  if (nameStr == "alice3RICHPIDEl") {
    cut->addCut(getAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDPi") {
    cut->addCut(getAnalysisCut("alice3RICHPIDPi"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDKa") {
    cut->addCut(getAnalysisCut("alice3RICHPIDKa"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDPr") {
    cut->addCut(getAnalysisCut("alice3RICHPIDPr"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDEl") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDPi") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPi"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDPi"));
    cut->addCut(getAnalysisCut("alice3RICHPIDPi"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDKa") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDKa"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDKa"));
    cut->addCut(getAnalysisCut("alice3RICHPIDKa"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDPr") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPr"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDPr"));
    cut->addCut(getAnalysisCut("alice3RICHPIDPr"));
    return cut;
  }

  if (nameStr == "alice3KaonTOFRICH") {
    cut->addCut(getAnalysisCut("alice3KaonTOFRICH"));
    return cut;
  }

  if (nameStr == "alice3PionTOFRICH") {
    cut->addCut(getAnalysisCut("alice3PionTOFRICH"));
    return cut;
  }

  if (nameStr == "alice3ProtonTOFRICH") {
    cut->addCut(getAnalysisCut("alice3ProtonTOFRICH"));
    return cut;
  }

  if (nameStr == "alice3ElectronTOFRICH") {
    cut->addCut(getAnalysisCut("alice3ElectronTOFRICH"));
    return cut;
  }

  if (nameStr == "alice3DielectronPID") {
    cut->addCut(getAnalysisCut("alice3JpsiKine"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3CharmoniumPID"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFOnly") {
    cut->addCut(getAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDRICHOnly") {
    cut->addCut(getAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFOnly") {
    cut->addCut(getAnalysisCut("alice3JpsiKine"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFAcceptance") {
    cut->addCut(getAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDRICHAcceptance") {
    cut->addCut(getAnalysisCut("alice3JpsiKineRICHAcceptance"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    cut->addCut(getAnalysisCut("alice3iTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDEl"));
    cut->addCut(getAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3LambdaCQualityCuts") {
    cut->addCut(getAnalysisCut("alice3LambdaCKine"));
    cut->addCut(getAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TOFPIDPi") {
    cut->addCut(getAnalysisCut("alice3iTOFPIDPi"));
    cut->addCut(getAnalysisCut("alice3oTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3TOFPIDKaRej") {
    cut->addCut(getAnalysisCut("alice3D0KaonTOFPID"));
    return cut;
  }

  if (nameStr == "alice3TOFPIDPiRej") {
    cut->addCut(getAnalysisCut("alice3D0PionTOFPID"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDPiRej") {
    cut->addCut(getAnalysisCut("alice3D0PionRICHTOFPID"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDKaRej") {
    cut->addCut(getAnalysisCut("alice3D0KaonRICHTOFPID"));
    return cut;
  }

  delete cut;
  LOGF(fatal, Form("Did not find cut %s. Returning nullptr", cutName));
  return nullptr;
}

AnalysisCut* o2::aod::dqcuts::getAnalysisCut(const char* cutName)
{
  //
  // define here cuts which are likely to be used often
  //
  auto* cut = new AnalysisCut(cutName, cutName);
  std::string nameStr = cutName;
  // ---------------------------------------------------------------
  // Event cuts
  if (nameStr == "noEventCut") {
    return cut;
  }

  if (nameStr == "eventNoTFBorder") {
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventIsTVXTriggered") {
    cut->addCut(VarManager::kIsTVXTriggered, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandard") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    return cut;
  }

  if (nameStr == "eventStandardtest") {
    cut->addCut(VarManager::kVtxZ, -30.0, 30.0);
    return cut;
  }

  if (nameStr == "eventSel8") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    return cut;
  }
  if (nameStr == "eventStandardSel8") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    return cut;
  }
  if (nameStr == "eventStandardSel8WithITSROFRecomputedCut") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorderRecomputed, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBorder") { // Redundant w.r.t. eventStandardSel8, to be removed
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBNoITSROFB") { // Redundant w.r.t. eventStandardSel8, to be removed
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBNoITSROFBrecomp") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorderRecomputed, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8NoSameBunch") {
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8TriggerZNAZNC") {
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsTriggerZNAZNC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8TriggerZNAZNCNoPileUp") {
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsTriggerZNAZNC, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8NoSameBunchGoodZvtx") {
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQuality") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoPileup") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityCent90") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityGoodITSLayersAll") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodITSLayersAll, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancy") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0., 1000);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityFirmTrackOccupancy") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0., 2000);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityLooseTrackOccupancy") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0., 5000);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancyCollInTime") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0., 500);
    cut->addCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancyCollInTime2") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0., 1000);
    cut->addCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);

    return cut;
  }

  std::vector<double> vecOccupancies = {0.,
                                        250.,
                                        500.,
                                        750.,
                                        1000.,
                                        1500.,
                                        2000.,
                                        3000.,
                                        4500.,
                                        6000.,
                                        8000.,
                                        10000.,
                                        50000.};

  for (size_t icase = 0; icase < vecOccupancies.size() - 1; icase++) {
    if (nameStr == Form("eventStandardSel8PbPbQualityTrackOccupancySlice%lu", icase)) {
      cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
      cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
      cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
      cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
      cut->addCut(VarManager::kTrackOccupancyInTimeRange, vecOccupancies[icase], vecOccupancies[icase + 1]);

      return cut;
    }
  }

  for (size_t icase = 0; icase < vecOccupancies.size() - 1; icase++) {
    if (nameStr == Form("eventStandardSel8PbPbQualityTrackOccupancySlice_0_%lu", icase)) {
      cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
      cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
      cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
      cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
      cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
      cut->addCut(VarManager::kTrackOccupancyInTimeRange, 0, vecOccupancies[icase]);

      return cut;
    }
  }

  if (nameStr == "eventStandardSel8ppQuality") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8ppQualityNoVtxZ") {
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8multAnalysis") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8VtxQuality1") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->addCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbMultCorr") {
    std::shared_ptr<TF1> fMultPVCutLow = std::make_shared<TF1>("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
    std::shared_ptr<TF1> fMultPVCutHigh = std::make_shared<TF1>("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

    std::shared_ptr<TF1> fMultCutLow = std::make_shared<TF1>("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    std::shared_ptr<TF1> fMultCutHigh = std::make_shared<TF1>("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    cut->addCut(VarManager::kVtxNcontribReal, fMultPVCutLow, fMultPVCutHigh, false, VarManager::kCentFT0C, 0.0, 100.0, false);
    cut->addCut(VarManager::kMultA, fMultCutLow, fMultCutHigh, false, VarManager::kCentFT0C, 0.0, 100.0, false);
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->addCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->addCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventDimuonStandard") {
    cut->addCut(VarManager::kIsMuonUnlikeLowPt7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventMuonStandard") {
    cut->addCut(VarManager::kIsMuonSingleLowPt7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventTPCMultLow") {
    cut->addCut(VarManager::kMultTPC, 0, 50);
    return cut;
  }

  if (nameStr == "eventExclusivePair") {
    cut->addCut(VarManager::kVtxNcontrib, 2, 2);
    return cut;
  }

  if (nameStr == "eventVtxNContrib") {
    cut->addCut(VarManager::kVtxNcontrib, 0, 10);
    return cut;
  }

  if (nameStr == "eventTPCMult3") {
    cut->addCut(VarManager::kMultTPC, 3, 3);
    return cut;
  }

  if (nameStr == "int7vtxZ5") {
    cut->addCut(VarManager::kVtxZ, -5.0, 5.0);
    cut->addCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventDoubleGap") {
    cut->addCut(VarManager::kIsDoubleGap, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGap") {
    cut->addCut(VarManager::kIsSingleGap, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapA") {
    cut->addCut(VarManager::kIsSingleGapA, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapAZDC") {
    cut->addCut(VarManager::kIsSingleGapA, 0.5, 1.5);
    cut->addCut(VarManager::kEnergyCommonZNA, -1000., 1.);
    cut->addCut(VarManager::kEnergyCommonZNC, 1., 1000.);
    return cut;
  }

  if (nameStr == "eventSingleGapC") {
    cut->addCut(VarManager::kIsSingleGapC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapCZDC") {
    cut->addCut(VarManager::kIsSingleGapC, 0.5, 1.5);
    cut->addCut(VarManager::kEnergyCommonZNC, -1000., 1.);
    cut->addCut(VarManager::kEnergyCommonZNA, 1., 1000.);
    return cut;
  }

  if (nameStr == "eventSingleGapACZDC") {
    auto* cutA = new AnalysisCompositeCut("singleGapAZDC", "singleGapAZDC", true);
    cutA->addCut(getAnalysisCut("eventSingleGapAZDC"));

    auto* cutC = new AnalysisCompositeCut("singleGapCZDC", "singleGapCZDC", true);
    cutC->addCut(getAnalysisCut("eventSingleGapCZDC"));

    auto* cutAorC = new AnalysisCompositeCut("singleGapACZDC", "singleGapACZDC", false);
    cutAorC->addCut(cutA);
    cutAorC->addCut(cutC);
    return cutAorC;
  }

  if (nameStr == "eventUPCMode") {
    cut->addCut(VarManager::kIsITSUPCMode, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapACZDC_UPCMode") {
    auto* cutA = new AnalysisCompositeCut("singleGapAZDC", "singleGapAZDC", true);
    cutA->addCut(getAnalysisCut("eventSingleGapAZDC"));
    cutA->addCut(getAnalysisCut("eventUPCMode"));

    auto* cutC = new AnalysisCompositeCut("singleGapCZDC", "singleGapCZDC", true);
    cutC->addCut(getAnalysisCut("eventSingleGapCZDC"));
    cutC->addCut(getAnalysisCut("eventUPCMode"));

    auto* cutAorC = new AnalysisCompositeCut("singleGapACZDC", "singleGapACZDC", false);
    cutAorC->addCut(cutA);
    cutAorC->addCut(cutC);
    return cutAorC;
  }

  if (nameStr == "eventXn0nTime") {
    cut->addCut(VarManager::kTimeZNA, -2.0, 2.0);
    cut->addCut(VarManager::kTimeZNC, -2.0, 2.0, true);
    return cut;
  }

  if (nameStr == "event0nXnTime") {
    cut->addCut(VarManager::kTimeZNA, -2.0, 2.0, true);
    cut->addCut(VarManager::kTimeZNC, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "event0n0nTime") {
    cut->addCut(VarManager::kTimeZNA, -2.0, 2.0, true);
    cut->addCut(VarManager::kTimeZNC, -2.0, 2.0, true);
    return cut;
  }

  if (nameStr == "eventXnXnTime") {
    cut->addCut(VarManager::kTimeZNA, -2.0, 2.0);
    cut->addCut(VarManager::kTimeZNC, -2.0, 2.0);
    return cut;
  }

  // Event cuts based on centrality
  if (nameStr == "eventStandardNoINT7Cent090") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 0.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent7090") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 70.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent5070") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 50.0, 70.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent3050") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 30.0, 50.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent1030") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 10.0, 30.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent010") {
    cut->addCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->addCut(VarManager::kCentFT0C, 0.0, 10.0);
    return cut;
  }

  // ---------------------------------------------------
  // Barrel track kine cuts
  if (nameStr == "negTrack") {
    cut->addCut(VarManager::kCharge, -99., 0.);
    return cut;
  }

  if (nameStr == "posTrack") {
    cut->addCut(VarManager::kCharge, 0., 99.);
    return cut;
  }

  if (nameStr == "posEtaSel") {
    cut->addCut(VarManager::kEta, 0., 0.8);
    return cut;
  }

  if (nameStr == "negEtaSel") {
    cut->addCut(VarManager::kEta, -0.8, 0.);
    return cut;
  }

  if (nameStr == "etaSel") {
    cut->addCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "pt02Sel") {
    cut->addCut(VarManager::kPt, 0.2, 20.0);
    return cut;
  }

  if (nameStr == "pt04Sel") {
    cut->addCut(VarManager::kPt, 0.4, 20.0);
    return cut;
  }

  if (nameStr == "openEtaSel") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "insideTPCsector") {
    cut->addCut(VarManager::kTrackIsInsideTPCModule, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "rho0Kine") {
    cut->addCut(VarManager::kPt, 0.1, 1000.0);
    cut->addCut(VarManager::kEta, -1.1, 1.1);
    return cut;
  }

  if (nameStr == "pionQuality") {
    cut->addCut(VarManager::kTPCncls, 50.0, 1000.);
    return cut;
  }

  if (nameStr == "primaryVertexContributor") {
    cut->addCut(VarManager::kPVContributor, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "jpsiStandardKine") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine2") {
    cut->addCut(VarManager::kPt, 0.9, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine3") {
    cut->addCut(VarManager::kP, 1.2, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiPIDcalibKine_posEta") {
    cut->addCut(VarManager::kPin, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, 0.0, 0.9);
    return cut;
  }

  if (nameStr == "jpsiPIDcalibKine_negEta") {
    cut->addCut(VarManager::kPin, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.0);
    return cut;
  }

  if (nameStr == "jpsiStandardKine4") {
    cut->addCut(VarManager::kP, 1.5, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine5") {
    cut->addCut(VarManager::kP, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiKineSkimmed") {
    cut->addCut(VarManager::kPt, 0.7, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "lmeePrefilterKine") {
    cut->addCut(VarManager::kPt, 0., 20.0);
    cut->addCut(VarManager::kEta, -1.2, 1.2);
    return cut;
  }

  if (nameStr == "lmeeStandardKine") {
    cut->addCut(VarManager::kPt, 0.2, 20.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "lmeeStandardKine_pt04") {
    cut->addCut(VarManager::kPt, 0.4, 20.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "lmeeLowBKine") {
    cut->addCut(VarManager::kPt, 0.075, 20.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "dalitzStandardKine") {
    cut->addCut(VarManager::kPt, 0.15, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "PIDStandardKine") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "PIDStandardKine2") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kPt, 0.1, 1000.0);
    return cut;
  }

  if (nameStr == "PIDStandardKine3") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kPt, 0.5, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow04") {
    cut->addCut(VarManager::kPt, 0.4, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow03") {
    cut->addCut(VarManager::kPt, 0.3, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow02") {
    cut->addCut(VarManager::kPt, 0.2, 1000.0);
    return cut;
  }

  // -----------------------------------------------
  // Barrel track quality cuts

  // ---------------------------------------------------
  // MC generated particle acceptance cuts

  if (nameStr == "rapidity08") {
    cut->addCut(VarManager::kMCY, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "rapidity09") {
    cut->addCut(VarManager::kMCY, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_pp13600") {
    cut->addCut(VarManager::kMCY, -0.8, 0.8);
    cut->addCut(VarManager::kMCPt1, 1.0, 1000.0);
    cut->addCut(VarManager::kMCPt2, 1.0, 1000.0);
    cut->addCut(VarManager::kMCEta1, -0.8, 0.8);
    cut->addCut(VarManager::kMCEta2, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "acceptance_pp5360") {
    cut->addCut(VarManager::kMCY, -0.9, 0.9);
    cut->addCut(VarManager::kMCPt1, 1.0, 1000.0);
    cut->addCut(VarManager::kMCPt2, 1.0, 1000.0);
    cut->addCut(VarManager::kMCEta1, -0.9, 0.9);
    cut->addCut(VarManager::kMCEta2, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_PbPb5360") {
    cut->addCut(VarManager::kMCY, -0.9, 0.9);
    cut->addCut(VarManager::kMCP1, 1.0, 1000.0);
    cut->addCut(VarManager::kMCP2, 1.0, 1000.0);
    cut->addCut(VarManager::kMCEta1, -0.9, 0.9);
    cut->addCut(VarManager::kMCEta2, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_PbPb5360_y08") {
    cut->addCut(VarManager::kMCY, -0.8, 0.8);
    cut->addCut(VarManager::kMCP1, 1.0, 1000.0);
    cut->addCut(VarManager::kMCP2, 1.0, 1000.0);
    cut->addCut(VarManager::kMCEta1, -0.8, 0.8);
    cut->addCut(VarManager::kMCEta2, -0.8, 0.8);
    return cut;
  }

  // ---------------------------------------------------
  // MC generated particle acceptance cuts

  // Run 2 only

  if (nameStr == "highPtHadron") {
    cut->addCut(VarManager::kPt, 4.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->addCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  if (nameStr == "TightGlobalTrack") {
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->addCut(VarManager::kTPCncls, 60.0, 170.);
    cut->addCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQuality") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->addCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityBenchmark") {
    cut->addCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->addCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  // Run 2 or run 3

  if (nameStr == "jpsi_trackCut_debug") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 90., 159);
    cut->addCut(VarManager::kITSncls, 2.5, 7.5);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_noITSCuts_debug") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 90., 159);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug2") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 90., 159);
    cut->addCut(VarManager::kITSncls, 2.5, 7.5);
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kTrackDCAxy, -1, 1);
    cut->addCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug3") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 90., 159);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCnclsCR, 70., 159);
    cut->addCut(VarManager::kITSncls, 2.5, 7.5);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kTrackDCAxy, -0.05, 0.05);
    cut->addCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug4") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 90., 159);
    cut->addCut(VarManager::kITSncls, 2.5, 7.5);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTrackDCAxy, -1, 1);
    cut->addCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug5") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70., 159);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug6") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 120., 159);
    cut->addCut(VarManager::kTPCnclsCR, 140., 159);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "lmee_trackCut_debug") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 80., 159);
    cut->addCut(VarManager::kITSncls, 0.0, 7.5);
    return cut;
  }

  if (nameStr == "lmee_skimming_cuts") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCnclsCR, 60.0, 161.);
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
    return cut;
  }

  if (nameStr == "trackQuality_compareDQEMframework") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kITSncls, 4.5, 7.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->addCut(VarManager::kTPCnCRoverFindCls, 0.8, 1e+10);
    cut->addCut(VarManager::kTPCncls, 90.0, 170.);
    return cut;
  }

  if ((nameStr == "TightGlobalTrackRun3") || (nameStr == "lmeeQCTrackCuts")) {
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kITSncls, 4.5, 7.5);
    cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->addCut(VarManager::kTPCncls, 90.0, 170.);
    return cut;
  }

  std::vector<TString> vecTypetrack;
  vecTypetrack.emplace_back("");                  // default TightGlobalTrackRun3 as above
  vecTypetrack.emplace_back("_7ITSncls");         // default TightGlobalTrackRun3 but with 7 ITS clusters
  vecTypetrack.emplace_back("_ITS");              // Ask only for ITS requirements
  vecTypetrack.emplace_back("_ITSalone");         // Ask only for ITS requirements + ITSalone (no TPC matching)
  vecTypetrack.emplace_back("_TPC");              // Ask only for TPC requirements
  vecTypetrack.emplace_back("_TPCalone");         // Ask only for TPC requirements + TPCalone (no ITS matching)
  vecTypetrack.emplace_back("_TPCnoTRD");         // Ask only for TPC requirements no TRD matching
  vecTypetrack.emplace_back("_TPCstrongncls");    // default TightGlobalTrackRun3 but with 130 TPC clusters
  vecTypetrack.emplace_back("_ITSanyfirsttwo");   // default TightGlobalTrackRun3 but with a cluster in any of the first two layers
  vecTypetrack.emplace_back("_ITSanyfirstthree"); // default TightGlobalTrackRun3 but with a cluster in any of the first three layers

  // loop to define PID cuts with and without post calibration
  for (size_t icase = 1; icase < vecTypetrack.size(); icase++) {
    if (nameStr == Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())) {
      if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 6.5, 7.5);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 4.5, 7.5);
      } else if (icase == 3) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 4.5, 7.5);
        cut->addCut(VarManager::kHasTPC, -0.5, 0.5);
      } else if (icase == 4) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
      } else if (icase == 5) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
        cut->addCut(VarManager::kHasITS, -0.5, 0.5);
      } else if (icase == 6) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
        cut->addCut(VarManager::kHasTRD, -0.5, 0.5);
      } else if (icase == 7) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 4.5, 7.5);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 130.0, 170.);
      } else if (icase == 8) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 4.5, 7.5);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
      } else {
        cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
        cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->addCut(VarManager::kITSncls, 4.5, 7.5);
        cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->addCut(VarManager::kTPCncls, 90.0, 170.);
      }
      return cut;
    }
  }

  if (nameStr == "LooseGlobalTrackRun3") {
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.0, 6.0);
    cut->addCut(VarManager::kITSncls, 3.5, 7.5);
    cut->addCut(VarManager::kTPCncls, 70.0, 170.);
    return cut;
  }

  if (nameStr == "TightGlobalTrackRun3_strongTPC") {
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kITSncls, 4.5, 7.5);
    cut->addCut(VarManager::kTPCnclsCR, 140.0, 161.);
    cut->addCut(VarManager::kTPCncls, 120.0, 170.);
    return cut;
  }

  if (nameStr == "TightTPCTrackRun3") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
    return cut;
  }

  if (nameStr == "TightTPCTrack") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "SPDfirst") {
    cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "noTPC") {
    cut->addCut(VarManager::kHasTPC, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "SPDany") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "ITSiball") {
    cut->addCut(VarManager::kIsITSibAll, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "ITSibany") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    return cut;
  }

  // List of 30 variations in ITS and TPC parameters
  std::vector<double> cutVarITSchi2 = {6., 6., 5., 4., 4., 6., 6., 5., 4., 5., 4., 5., 6., 5., 6., 5., 6., 5., 5., 4., 6., 4., 6., 5., 6., 4., 4., 6., 4., 5.};
  std::vector<double> cutVarTPCchi2 = {5., 5., 4., 3., 5., 4., 5., 3., 5., 4., 5., 3., 3., 5., 4., 5., 3., 5., 5., 5., 3., 5., 5., 4., 3., 4., 5., 5., 5., 3.};
  std::vector<double> cutVarITSnCl = {4.5, 5.5, 5.5, 4.5, 6.5, 4.5, 4.5, 4.5, 4.5, 5.5, 4.5, 5.5, 5.5, 5.5, 5.5, 4.5, 5.5, 6.5, 5.5, 4.5, 4.5, 5.5, 5.5, 5.5, 6.5, 5.5, 4.5, 4.5, 6.5, 6.5};
  std::vector<double> cutVarTPCnClsCR = {90., 80., 80., 80., 90., 80., 70., 90., 70., 80., 70., 90., 90., 70., 90., 90., 70., 80., 90., 80., 80., 90., 70., 70., 70., 80., 90., 70., 70., 80.};
  std::vector<double> cutVarTPCnCls = {80., 100., 80., 90., 90., 80., 80., 80., 80., 90., 100., 100., 80., 80., 80., 80., 100., 90., 100., 90., 90., 100., 100., 80., 100., 90., 90., 100., 90., 90.};

  for (unsigned int i = 0; i < cutVarITSchi2.size(); i++) {
    if (nameStr == Form("lmeeCutVarTrackCuts%i", i)) {
      cut->addCut(VarManager::kIsSPDfirst, 0.5, 1.5);
      cut->addCut(VarManager::kITSchi2, 0.0, cutVarITSchi2.at(i));
      cut->addCut(VarManager::kTPCchi2, 0.0, cutVarTPCchi2.at(i));
      cut->addCut(VarManager::kITSncls, cutVarITSnCl.at(i), 7.5);
      cut->addCut(VarManager::kTPCnclsCR, cutVarTPCnClsCR.at(i), 161.);
      cut->addCut(VarManager::kTPCncls, cutVarTPCnCls.at(i), 170.);
      return cut;
    }
  }

  if (nameStr == "electronStandardQualityForO2MCdebug") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug2") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug3") {
    cut->addCut(VarManager::kITSncls, 0.5, 10);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug4") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityITSOnly") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQualitybAnyITSOnly") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly2") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 100, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly3") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 120, 161.);
    return cut;
  }

  if (nameStr == "NoelectronStandardQualityTPCOnly") {
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0, true, VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kTPCncls, 60, 161);
    cut->addCut(VarManager::kTrackDCAxy, -1.5, 1.5);
    cut->addCut(VarManager::kTrackDCAz, -1.5, 1.5);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed2") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kTPCncls, 60, 161);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed3") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCnclsCR, 70, 161);
    cut->addCut(VarManager::kTPCncls, 70, 161);
    return cut;
  }

  if (nameStr == "electronTrackQuality_Maolin") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kITSchi2, 0.0, 15.0);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 70, 161.);
    cut->addCut(VarManager::kTPCnclsCR, 70, 161);
    cut->addCut(VarManager::kTrackDCAxy, -2.0, 2.0);
    cut->addCut(VarManager::kTrackDCAz, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "pionQualityCut1") {
    cut->addCut(VarManager::kPt, 0.15, 1000.0);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCncls, 70, 161);
    return cut;
  }

  if (nameStr == "pionQualityCut2") {
    cut->addCut(VarManager::kPt, 0.15, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCncls, 90, 161);
    cut->addCut(VarManager::kTPCnclsCR, 70, 161);
    return cut;
  }

  if (nameStr == "protonPVcut") {
    cut->addCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->addCut(VarManager::kTrackDCAz, -0.15, 0.15);
    cut->addCut(VarManager::kPt, 0.4, 3);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kTPCncls, 80, 161.);
    return cut;
  }

  if (nameStr == "pidbasic") {
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    cut->addCut(VarManager::kTPCncls, 60, 161.);
    return cut;
  }

  if (nameStr == "standardPrimaryTrack") {
    cut->addCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->addCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "trackDCA1cm") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->addCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->addCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "electronPrimary_dca3sigma" || nameStr == "electronPrimary_dca7sigma") {
    std::shared_ptr<TF1> fDCAxyresLow = std::make_shared<TF1>("fDCAxyresLow", "[0] + [1] * pow(x, -[2])", 0.1, 1000.);
    std::shared_ptr<TF1> fDCAzresLow = std::make_shared<TF1>("fDCAzresLow", "[0] + [1] * pow(x, -[2])", 0.1, 1000.);
    std::shared_ptr<TF1> fDCAxyresUp = std::make_shared<TF1>("fDCAxyresUp", "[0] + [1] * pow(x, -[2])", 0.1, 1000.);
    std::shared_ptr<TF1> fDCAzresUp = std::make_shared<TF1>("fDCAzresUp", "[0] + [1] * pow(x, -[2])", 0.1, 1000.);

    if (nameStr == "electronPrimary_dca3sigma") {
      // DCAxy and DCAz 3 sigma cut. DCA resolution vs pt extracted from fits of Users/m/mfaggin/test/inputsTrackTuner/pp2024/pass1_minBias/vsPhi (used for the track tuner)
      // we add in addition a term for the misalignment of the mean of the distribution, which seems to be at most 20 mum for DCAxy and 10 mum for DCAz
      fDCAxyresLow->SetParameters(-3 * 8.7e-4 - 20e-4, -3 * 25.4e-4, 0.79); // res is 8.7 + 25.4/pt^0.79 mum
      fDCAzresLow->SetParameters(-3 * 9.4e-4 - 10e-4, -3 * 26.5e-4, 0.79);  // res is 9.4 + 26.5/pt^0.79 mum
      fDCAxyresUp->SetParameters(3 * 8.7e-4 + 20e-4, 3 * 25.4e-4, 0.79);    // res is 8.7 + 25.4/pt^0.79 mum
      fDCAzresUp->SetParameters(3 * 9.4e-4 + 10e-4, 3 * 26.5e-4, 0.79);     // res is 9.4 + 26.5/pt^0.79 mum
      cut->addCut(VarManager::kTrackDCAxy, fDCAxyresLow, fDCAxyresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->addCut(VarManager::kTrackDCAz, fDCAzresLow, fDCAzresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->addCut(VarManager::kIsITSibFirst, 0.5, 1.5);
      cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
      cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
      cut->addCut(VarManager::kTPCncls, 100, 161.);
      return cut;
    }

    if (nameStr == "electronPrimary_dca7sigma") {
      // DCAxy and DCAz 7 sigma cut
      fDCAxyresLow->SetParameters(-7 * 8.7e-4 - 20e-4, -7 * 25.4e-4, 0.79);
      fDCAzresLow->SetParameters(-7 * 9.4e-4 - 10e-4, -7 * 26.5e-4, 0.79);
      fDCAxyresUp->SetParameters(7 * 8.7e-4 + 20e-4, 7 * 25.4e-4, 0.79);
      fDCAzresUp->SetParameters(7 * 9.4e-4 + 10e-4, 7 * 26.5e-4, 0.79);
      cut->addCut(VarManager::kTrackDCAxy, fDCAxyresLow, fDCAxyresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->addCut(VarManager::kTrackDCAz, fDCAzresLow, fDCAzresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->addCut(VarManager::kIsITSibFirst, 0.5, 1.5);
      cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
      cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
      cut->addCut(VarManager::kTPCncls, 100, 161.);
      return cut;
    }
  }

  if (nameStr == "dcaCut1_ionut") {
    cut->addCut(VarManager::kTrackDCAxy, -0.5, 0.5);
    cut->addCut(VarManager::kTrackDCAz, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "trackQuality_ionut") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCncls, 70, 161);
    cut->addCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->addCut(VarManager::kTPCchi2, 0.0, 2.0);
    return cut;
  }

  if (nameStr == "trackQualityTight_ionut") {
    cut->addCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->addCut(VarManager::kTPCncls, 100, 161);
    cut->addCut(VarManager::kITSchi2, 0.0, 3.0);
    cut->addCut(VarManager::kTPCchi2, 0.0, 2.0);
    cut->addCut(VarManager::kITSncls, 5.0, 8.0);
    return cut;
  }

  if (nameStr == "kineJpsiEle_ionut") {
    cut->addCut(VarManager::kP, 1.0, 15.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "pidJpsiEle0_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 1.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kPin, 4.0, 150.0);
    cut->addCut(VarManager::kTPCnSigmaEl, 98.1, 98.11, false, VarManager::kPin, 0.0, 1.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle1_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 1.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kPin, 4.0, 150.0);
    cut->addCut(VarManager::kTPCnSigmaEl, 98.1, 98.11, false, VarManager::kPin, 0.0, 1.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle2_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle3_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, -0.5, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle4_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, 0.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle5_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, 0.5, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle6_ionut") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    cut->addCut(VarManager::kTOFnSigmaEl, -1.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle7_ionut") {
    cut->addCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsiEle8_ionut") {
    cut->addCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -1.5, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsiEle9_ionut") {
    cut->addCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsi_TPCpion0") {
    cut->addCut(VarManager::kTPCnSigmaPi, 4.0, 1000.0);
    return cut;
  }

  if (nameStr == "pidJpsi_noTOF_prot") {
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 1000.0, false, VarManager::kHasTOF, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "pidJpsi_beta") {
    cut->addCut(VarManager::kTOFbeta, 0.98, 1.02, false, VarManager::kHasTOF, 0.5, 1.5);
    return cut;
  }

  // Magnus cuts ----------------------------------------------------------

  if (nameStr == "pidJpsi_magnus_ele1") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_ele2") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_ele3") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_prot1") {
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_prot2") {
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_pion1") {
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_pion2") {
    cut->addCut(VarManager::kTPCnSigmaPi, 3.5, 1000.0);
    return cut;
  }

  // ----------------------------------------------------------------------------------

  if (nameStr == "standardPrimaryTrackDCAz") {
    cut->addCut(VarManager::kTrackDCAxy, -3.0, 3.0);
    cut->addCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "standardPrimaryTrackDCA") {
    cut->addCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->addCut(VarManager::kTrackDCAz, -0.15, 0.15);
    return cut;
  }

  if (nameStr == "PrimaryTrack_looseDCA") {
    cut->addCut(VarManager::kTrackDCAxy, -3.0, 3.0);
    cut->addCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "tightPrimaryTrack") {
    cut->addCut(VarManager::kTrackDCAsigXY, -3.0, 3.0);
    cut->addCut(VarManager::kTrackDCAsigZ, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCA05") {
    cut->addCut(VarManager::kTrackDCAsigXY, -0.5, 0.5);
    cut->addCut(VarManager::kTrackDCAsigZ, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCAz") {
    cut->addCut(VarManager::kTrackDCAz, -0.3, 0.3);
    return cut;
  }

  if (nameStr == "hasTOF") {
    cut->addCut(VarManager::kHasTOF, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "noTOF") {
    cut->addCut(VarManager::kHasTOF, -0.5, 0.5);
    return cut;
  }

  // -----------------------------------------------------
  // V0 and Dalitz legs selections

  for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("dalitzLeg%d", i)) {
      cut->addCut(VarManager::kIsDalitzLeg + i - 1, -0.5, 0.5, true);
      return cut;
    }

    if (nameStr == Form("notDalitzLeg%d", i)) {
      cut->addCut(VarManager::kIsDalitzLeg + i - 1, -0.5, 0.5);
      return cut;
    }
  }

  if (nameStr == "pidcalib_ele") {
    cut->addCut(VarManager::kIsLegFromGamma, -0.5, 0.5, true);
    return cut;
  }

  if (nameStr == "pidcalib_pion") {
    cut->addCut(VarManager::kIsLegFromK0S, -0.5, 0.5, true);
    return cut;
  }

  if (nameStr == "pidcalib_proton") {
    cut->addCut(VarManager::kIsProtonFromLambdaAndAntiLambda, -0.5, 0.5, true);
    return cut;
  }

  if (nameStr == "pidcalib_kaon") {
    cut->addCut(VarManager::kTOFnSigmaKa, -2.0, 2.0);
    cut->addCut(VarManager::kTOFnSigmaPi, -2.0, 2.0, true);
    cut->addCut(VarManager::kITSncls, 1.5, 7.5);
    return cut;
  }

  // ------------------------------------------------
  // Barrel PID cuts
  if (nameStr == "jpsi_TPCPID_debug1") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug2") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug3") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 3.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 3.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug4") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug5") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug5_noCorr") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 999);
    return cut;
  }
  if (nameStr == "electronPIDLooseSkimmed") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 999);
    return cut;
  }

  if (nameStr == "electronPIDLooseSkimmed2") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 1.5, 999, false, VarManager::kPin, 3.0, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 1.5, 999, false, VarManager::kPin, 3.0, 999);
    return cut;
  }

  if (nameStr == "electronPIDLooseSkimmed3") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kPin, 0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug6") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug7") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.5, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug8") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 3.0, false, VarManager::kPin, 0.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0, false, VarManager::kPin, 3.0, 999.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kPin, 0.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kPin, 5.0, 999.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 3.0, 999.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug9") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.5, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 1.0, 999, false, VarManager::kPin, 3.0, 999.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug10") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.5, 2.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 4.0, 999);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 999);
    return cut;
  }

  if (nameStr == "pidCut_lowP_Corr") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0, false, VarManager::kP, 0.0, 5.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kP, 0.0, 5.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kP, 0.0, 5.0);
    return cut;
  }

  if (nameStr == "EleInclusion_highP_Corr") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1.0, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "EleInclusion_highP2_Corr") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -0.5, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "PionExclusion_highP_Corr") {
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 2.0, 999, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "pidCut_lowP") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0, false, VarManager::kP, 0.0, 5.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kP, 0.0, 5.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kP, 0.0, 5.0);
    return cut;
  }

  if (nameStr == "EleInclusion_highP") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "EleInclusion_highP2") {
    cut->addCut(VarManager::kTPCnSigmaEl, -0.5, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "PionExclusion_highP") {
    cut->addCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }

  if (nameStr == "lmee_TPCPID_debug1") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
    return cut;
  }

  if (nameStr == "electronPID1" || nameStr == "electronPID1shiftUp" || nameStr == "electronPID1shiftDown" || nameStr == "electronPID2" || nameStr == "electronPID3") {
    std::shared_ptr<TF1> cutLow1 = std::make_shared<TF1>("cutLow1", "pol1", 0., 10.);
    if (nameStr == "electronPID1") {
      cutLow1->SetParameters(130., -40.0);
      cut->addCut(VarManager::kTPCsignal, 70., 100.);
      cut->addCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID1shiftUp") {
      cut->addCut(VarManager::kTPCsignal, 70. - 0.85, 100. - 0.85);
      cutLow1->SetParameters(130. - 0.85, -40.0);
      cut->addCut(VarManager::kTPCsignal, cutLow1, 100.0 - 0.85, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID1shiftDown") {
      cut->addCut(VarManager::kTPCsignal, 70.0 + 0.85, 100.0 + 0.85);
      cutLow1->SetParameters(130. + 0.85, -40.0);
      cut->addCut(VarManager::kTPCsignal, cutLow1, 100.0 + 0.85, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID2") {
      cutLow1->SetParameters(130., -40.0);
      cut->addCut(VarManager::kTPCsignal, 73., 100.);
      cut->addCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID3") {
      cutLow1->SetParameters(130., -40.0);
      cut->addCut(VarManager::kTPCsignal, 60., 110.);
      cut->addCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }
  }

  if (nameStr == "electronPIDnsigma") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0, 3000.0);
    return cut;
  }

  if (nameStr == "lmee_commonDQEM_PID_TPC") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->addCut(VarManager::kTPCnSigmaEl, -2.5, 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->addCut(VarManager::kTPCnSigmaPi, -1e12, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
    cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->addCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "lmee_commonDQEM_PID_TOF") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->addCut(VarManager::kTPCnSigmaEl, -2.5, 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
    cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
    return cut;
  }

  std::vector<TString> vecPIDcase;
  vecPIDcase.emplace_back("");              // without post calibration
  vecPIDcase.emplace_back("_Corr");         // case of using post calibrated PID spectra
  vecPIDcase.emplace_back("_CorrWithKaon"); // case of using post calibrated PID spectra with also the kaons

  // loop to define TPC PID cuts with and without post calibration
  for (size_t icase = 0; icase < vecPIDcase.size(); icase++) {
    if (nameStr == Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_loose", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TPCnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    // List of nSigma values for lower and upper edge of El, Pi, Ka, Pr, TPC PID selections/rejection
    std::vector<double> cutVarTPCnSigmaElLow = {-4., -4., -4., -2., -3., -2., -3., -4., -2., -4., -3., -3., -3., -4., -3., -3., -4., -3., -4., -4., -3., -4., -3., -3., -2., -4., -4., -3., -4., -2};
    std::vector<double> cutVarTPCnSigmaElUp = {2., 3., 2., 2., 2., 2., 2., 4., 2., 2., 4., 3., 3., 2., 4., 2., 2., 4., 3., 4., 2., 4., 2., 3., 2., 3., 4., 2., 3., 2};
    std::vector<double> cutVarTPCnSigmaPiLow = {-3., -2., -3., -4., -4., -3., -4., -2., -2., -2., -3., -3., -2., -2., -4., -3., -3., -2., -3., -2., -4., -2., -4., -4., -3., -3., -3., -2., -4., -4};
    std::vector<double> cutVarTPCnSigmaPiUp = {3., 3., 4., 4., 3., 2., 4., 4., 3., 4., 4., 3., 4., 4., 3., 2., 4., 2., 4., 2., 3., 4., 2., 2., 3., 2., 3., 4., 2., 4};
    std::vector<double> cutVarTPCnSigmaKaLow = {-4., -2., -2., -2., -4., -3., -2., -4., -3., -3., -4., -2., -3., -3., -4., -2., -4., -2., -3., -4., -4., -2., -2., -3., -2., -2., -3., -3., -2., -4};
    std::vector<double> cutVarTPCnSigmaKaUp = {4., 3., 2., 3., 4., 3., 4., 4., 4., 4., 4., 4., 4., 4., 2., 4., 4., 2., 2., 4., 3., 3., 2., 4., 2., 4., 3., 3., 3., 3};
    std::vector<double> cutVarTPCnSigmaPrLow = {-4., -2., -2., -3., -4., -4., -3., -2., -2., -4., -4., -2., -3., -4., -2., -3., -3., -2., -3., -3., -2., -2., -2., -2., -2., -3., -2., -3., -3., -3};
    std::vector<double> cutVarTPCnSigmaPrUp = {2., 2., 3., 2., 3., 3., 3., 2., 4., 3., 3., 4., 4., 3., 4., 4., 3., 4., 2., 3., 4., 4., 3., 4., 3., 2., 3., 3., 2., 3};

    for (unsigned int i = 0; i < cutVarTPCnSigmaElLow.size(); i++) {
      if (nameStr == Form("electronPID_TPCnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
          cut->addCut(VarManager::kTPCnSigmaEl, cutVarTPCnSigmaElLow.at(i), cutVarTPCnSigmaElUp.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPi, cutVarTPCnSigmaPiLow.at(i), cutVarTPCnSigmaPiUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaKa, cutVarTPCnSigmaKaLow.at(i), cutVarTPCnSigmaKaUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPr, cutVarTPCnSigmaPrLow.at(i), cutVarTPCnSigmaPrUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
          cut->addCut(VarManager::kTPCnSigmaEl_Corr, cutVarTPCnSigmaElLow.at(i), cutVarTPCnSigmaElUp.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPi_Corr, cutVarTPCnSigmaPiLow.at(i), cutVarTPCnSigmaPiUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaKa, cutVarTPCnSigmaKaLow.at(i), cutVarTPCnSigmaKaUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPr_Corr, cutVarTPCnSigmaPrLow.at(i), cutVarTPCnSigmaPrUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
          cut->addCut(VarManager::kTPCnSigmaEl_Corr, cutVarTPCnSigmaElLow.at(i), cutVarTPCnSigmaElUp.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPi_Corr, cutVarTPCnSigmaPiLow.at(i), cutVarTPCnSigmaPiUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaKa_Corr, cutVarTPCnSigmaKaLow.at(i), cutVarTPCnSigmaKaUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPr_Corr, cutVarTPCnSigmaPrLow.at(i), cutVarTPCnSigmaPrUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        }
        return cut;
      }
    }

    if (nameStr == Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_lowB_TPC%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCloosenopkrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) {
        cut->addCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCPbPbnopkrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCloose%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }
  }

  if (nameStr == "electronPIDnsigmaOpen") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.0, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaVeryVeryLoose") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaVeryLoose") {
    cut->addCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.4, 1e+10, false);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaLoose") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaMedium") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.7, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.7, 3000.0);
    return cut;
  }
  if (nameStr == "electronPIDnsigmaMedium_withLargeTOFPID") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 2.7, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 2.7, 3000.0);
    cut->addCut(VarManager::kTOFnSigmaEl, -5.0, 5.0);
    return cut;
  }
  if (nameStr == "electronPIDnsigmaSkewed") {
    cut->addCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaSkewed_2") {
    cut->addCut(VarManager::kTPCnSigmaEl, -0.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRej") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -3.0, 3.0, true);
    cut->addCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true);
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRej_Corr") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3.0, 3.0, true);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true);
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRejLoose") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, -2.0, 2.0, true);
    cut->addCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->addCut(VarManager::kTPCnSigmaPi, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRejLoose_Corr") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, -2.0, 2.0, true);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaEMu") {
    cut->addCut(VarManager::kTPCnSigmaEl, -1.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaKa, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma") {
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "kaonRejNsigma") {
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma2") {
    cut->addCut(VarManager::kTPCnSigmaKa, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "kaonPID_TPCnTOF") {
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    cut->addCut(VarManager::kTOFnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma700") {
    cut->addCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    cut->addCut(VarManager::kPin, 0.0, 0.7);
    return cut;
  }

  if (nameStr == "AssocKine") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "electronPIDworseRes") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0 * 0.8, 3000.0); // emulates a 20% degradation in PID resolution
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0 * 0.8, 3000.0); // proton and pion rejections are effectively relaxed by 20%
    return cut;
  }

  if (nameStr == "electronPIDshift") {
    cut->addCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaPr, 3.0 - 0.2, 3000.0);
    cut->addCut(VarManager::kTPCnSigmaPi, 3.0 - 0.2, 3000.0);
    return cut;
  }

  if (nameStr == "pionPIDnsigma") {
    cut->addCut(VarManager::kTPCnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "pionPID_TPCnTOF") {
    cut->addCut(VarManager::kTPCnSigmaPi, -3.0, 3.0);
    cut->addCut(VarManager::kTOFnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "protonPID_TPCnTOF") {
    cut->addCut(VarManager::kTPCnSigmaPr, -3.0, 3.0);
    cut->addCut(VarManager::kTOFnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "protonPID_TPCnTOF2") {
    cut->addCut(VarManager::kTPCnSigmaPr, -2.5, 2.5);
    return cut;
  }

  if (nameStr == "tpc_pion_rejection") {
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+[1]*x", 0, 10);
    f1maxPi->SetParameters(85, -50);
    cut->addCut(VarManager::kTPCsignal, 70, f1maxPi, true, VarManager::kPin, 0.0, 0.4, false);
    return cut;
  }

  if (nameStr == "tpc_pion_band_rejection") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+[1]*log(x)", 0, 10);
    f1minPi->SetParameters(-115, -90);
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+[1]*log(x)", 0, 10);
    f1maxPi->SetParameters(-70, -90);
    cut->addCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.05, 0.3, false);
    return cut;
  }

  if (nameStr == "tpc_pion_muon_band_rejection") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPi->SetParameters(37, -18, 5.5);
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPi->SetParameters(67, -17, 5.9);
    cut->addCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_pion_rejection_highp") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+[1]*x", 0, 10);
    f1minPi->SetParameters(65, 4.);
    cut->addCut(VarManager::kTPCsignal, f1minPi, 110., false, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_kaon_rejection") {
    std::shared_ptr<TF1> f1minKa = std::make_shared<TF1>("f1minKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1minKa->SetParameters(37, -4, 5.6);
    std::shared_ptr<TF1> f1maxKa = std::make_shared<TF1>("f1maxKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxKa->SetParameters(60, -4.1, 6.);
    cut->addCut(VarManager::kTPCsignal, f1minKa, f1maxKa, true, VarManager::kPin, 0.0, 10.0, false);
    return cut;
  }

  if (nameStr == "tpc_proton_rejection") {
    std::shared_ptr<TF1> f1minPr = std::make_shared<TF1>("f1minPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPr->SetParameters(37, -2.6, 6.1);
    std::shared_ptr<TF1> f1maxPr = std::make_shared<TF1>("f1maxPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPr->SetParameters(60, -2.4, 6.2);
    cut->addCut(VarManager::kTPCsignal, f1minPr, f1maxPr, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_electron") {
    cut->addCut(VarManager::kTPCsignal, 60, 110, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron") {
    cut->addCut(VarManager::kTOFbeta, 0.99, 1.01, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron_sigma") {
    cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron_sigma_2") {
    cut->addCut(VarManager::kTOFnSigmaEl, -3., 3.);
    return cut;
  }

  if (nameStr == "tof_electron_loose") {
    cut->addCut(VarManager::kTOFbeta, 0.95, 1.05, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  // loop to define TOF PID cuts with and without post calibration
  for (size_t icase = 0; icase < vecPIDcase.size(); icase++) {
    if (nameStr == Form("electronPID_TOFnsigma%s_loose", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TOFnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    // List of nSigma values for lower and upper edge of TPC: El, Pi and TOF: El PID selections/rejection
    std::vector<double> cutVarTPCnSigmaElLow = {-4., -4., -4., -2., -3., -2., -3., -4., -2., -4., -3., -3., -3., -4., -3., -3., -4., -3., -4., -4., -3., -4., -3., -3., -2., -4., -4., -3., -4., -2};
    std::vector<double> cutVarTPCnSigmaElUp = {2., 3., 2., 2., 2., 2., 2., 4., 2., 2., 4., 3., 3., 2., 4., 2., 2., 4., 3., 4., 2., 4., 2., 3., 2., 3., 4., 2., 3., 2};
    std::vector<double> cutVarTPCnSigmaPiLow = {-3., -2., -3., -4., -4., -3., -4., -2., -2., -2., -3., -3., -2., -2., -4., -3., -3., -2., -3., -2., -4., -2., -4., -4., -3., -3., -3., -2., -4., -4};
    std::vector<double> cutVarTPCnSigmaPiUp = {3., 3., 4., 4., 3., 2., 4., 4., 3., 4., 4., 3., 4., 4., 3., 2., 4., 2., 4., 2., 3., 4., 2., 2., 3., 2., 3., 4., 2., 4};
    std::vector<double> cutVarTOFnSigmaElLow = {-4., -2., -4., -4., -3., -2., -4., -4., -4., -2., -2., -4., -3., -3., -4., -4., -4., -2., -4., -4., -2., -2., -3., -4., -4., -2., -4., -2., -3., -3};
    std::vector<double> cutVarTOFnSigmaElUp = {4., 2., 4., 2., 4., 3., 2., 3., 3., 3., 4., 3., 2., 3., 4., 3., 3., 3., 4., 4., 2., 2., 2., 3., 3., 3., 2., 3., 2., 4};

    for (unsigned int i = 0; i < cutVarTOFnSigmaElLow.size(); i++) {
      if (nameStr == Form("electronPID_TOFnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
          cut->addCut(VarManager::kTPCnSigmaEl, cutVarTPCnSigmaElLow.at(i), cutVarTPCnSigmaElUp.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPi, cutVarTPCnSigmaPiLow.at(i), cutVarTPCnSigmaPiUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTOFnSigmaEl, cutVarTOFnSigmaElLow.at(i), cutVarTOFnSigmaElUp.at(i), false, VarManager::kPin, 0.3, 1e+10, false);
        } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
          cut->addCut(VarManager::kTPCnSigmaEl_Corr, cutVarTPCnSigmaElLow.at(i), cutVarTPCnSigmaElUp.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTPCnSigmaPi_Corr, cutVarTPCnSigmaPiLow.at(i), cutVarTPCnSigmaPiUp.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->addCut(VarManager::kTOFnSigmaEl, cutVarTOFnSigmaElLow.at(i), cutVarTOFnSigmaElUp.at(i), false, VarManager::kPin, 0.3, 1e+10, false);
        }
        return cut;
      }
    }

    if (nameStr == Form("electronPID_TPC_TOFnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TOFnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->addCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->addCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFreq%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_lowB_TOF%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaMu, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 3., true, VarManager::kPin, 0.0, 4.0, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 4.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 3., true, VarManager::kPin, 0.0, 4.0, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 4.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 2.0, false);
        cut->addCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 2.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->addCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 2.0, false);
        cut->addCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 2.0, 1e+10, false);
        cut->addCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts
  if (nameStr == "GlobalMuonTrack") {
    cut->addCut(VarManager::kMuonTrackType, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "MFTMCH") {
    cut->addCut(VarManager::kMuonTrackType, 1.5, 2.5);
    return cut;
  }

  if (nameStr == "MCHMID") {
    cut->addCut(VarManager::kMuonTrackType, 2.5, 3.5);
    return cut;
  }

  if (nameStr == "MCHStandalone") {
    cut->addCut(VarManager::kMuonTrackType, 3.5, 4.5);
    return cut;
  }

  if (nameStr == "muonEtaCut") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    return cut;
  }

  if (nameStr == "muonMinimalCuts") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    return cut;
  }

  if (nameStr == "muonMinimalCuts10SigmaPDCA") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 990.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 540.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    return cut;
  }

  if (nameStr == "muonQualityCuts") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonQualityCuts5SigmaPDCA_Run3") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 500.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 335.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonQualityCuts10SigmaPDCA") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 990.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 540.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "matchedQualityCuts") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "matchedQualityCutsMFTeta") {
    cut->addCut(VarManager::kEta, -3.6, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "muonQualityCutsMatchingOnly") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonLowPt") {
    cut->addCut(VarManager::kPt, 0.5, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt2") {
    cut->addCut(VarManager::kPt, 0.7, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt3") {
    cut->addCut(VarManager::kPt, 0.8, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt4") {
    cut->addCut(VarManager::kPt, 0.9, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt5") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt6") {
    cut->addCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt1p5") {
    cut->addCut(VarManager::kPt, 1.5, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt1p7") {
    cut->addCut(VarManager::kPt, 1.7, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt") {
    cut->addCut(VarManager::kPt, 3.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt2") {
    cut->addCut(VarManager::kPt, 4.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt3") {
    cut->addCut(VarManager::kPt, 6.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt4") {
    cut->addCut(VarManager::kPt, 8.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt5") {
    cut->addCut(VarManager::kPt, 10.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt6") {
    cut->addCut(VarManager::kPt, 20.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonTightQualityCutsForTests") {
    cut->addCut(VarManager::kEta, -4.0, -2.5);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 20.0, 60.0);
    cut->addCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    return cut;
  }

  if (nameStr == "mchTrack") {
    cut->addCut(VarManager::kMuonTrackType, 3.5, 4.5);
    return cut;
  }

  if (nameStr == "matchedMchMid") {
    cut->addCut(VarManager::kMuonTrackType, 2.5, 3.5);
    return cut;
  }

  if (nameStr == "muonMchAndMchMid") {
    cut->addCut(VarManager::kMuonTrackType, 2.5, 4.5); // selecting both MCH and MCH-MID tracks
    return cut;
  }

  if (nameStr == "matchedFwd") {
    cut->addCut(VarManager::kMuonTrackType, 1.5, 2.5);
    return cut;
  }

  if (nameStr == "matchedGlobal") {
    cut->addCut(VarManager::kMuonTrackType, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut1") {
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0, 30);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut2") {
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0, 40);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut3") {
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0, 50);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut4") {
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0, 60);
    return cut;
  }

  // -----------------------------------------------------------------------------------------------
  // Pair cuts
  if (nameStr == "pairDalitz1") {
    cut->addCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->addCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    std::shared_ptr<TF1> fcutHigh = std::make_shared<TF1>("f1", "[0] - [0]/[1]*x", -1.5, 1.5);
    fcutHigh->SetParameters(0.6, 0.12);
    std::shared_ptr<TF1> fcutLow = std::make_shared<TF1>("f2", "-[0] + [0]/[1]*x", -1.5, 1.5);
    fcutLow->SetParameters(0.6, 0.12);
    cut->addCut(VarManager::kPsiPair, fcutLow, fcutHigh, true, VarManager::kDeltaPhiPair, 0, 0.12);
    return cut;
  }

  if (nameStr == "pairDalitz1Strong") {
    cut->addCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->addCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    cut->addCut(VarManager::kDeltaPhiPair, -1., 0.);
    return cut;
  }

  if (nameStr == "pairDalitz2") {
    cut->addCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->addCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    return cut;
  }

  if (nameStr == "pairDalitz3") {
    cut->addCut(VarManager::kMass, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "paira_prefilter1") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    return cut;
  }

  if (nameStr == "paira_prefilter2") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "paira_prefilter3") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "paira_prefilter4") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "paira_prefilter5") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "paira_prefilter6") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "paira_prefilter7") {
    cut->addCut(VarManager::kMass, 0.0, 0.06);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairb_prefilter1") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairb_prefilter2") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "pairb_prefilter3") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairb_prefilter4") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairb_prefilter5") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairb_prefilter6") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairb_prefilter7") {
    cut->addCut(VarManager::kMass, 0.0, 0.05);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairc_prefilter1") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    return cut;
  }

  if (nameStr == "pairc_prefilter2") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "pairc_prefilter3") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairc_prefilter4") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairc_prefilter5") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairc_prefilter6") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairc_prefilter7") {
    cut->addCut(VarManager::kMass, 0.0, 0.04);
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "paird_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "paire_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairf_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairg_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairh_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairi_prefilter1") {
    cut->addCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairNoCut") {
    cut->addCut(VarManager::kMass, 0.0, 1000.0);
    return cut;
  }

  if (nameStr == "DipionMassCut1") {
    cut->addCut(VarManager::kMass, 0.5, 1.0);
    return cut;
  }

  if (nameStr == "DipionMassCut2") {
    cut->addCut(VarManager::kMass, 0.0, 1.0);
    return cut;
  }

  if (nameStr == "pairMassLow1") {
    cut->addCut(VarManager::kMass, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow2") {
    cut->addCut(VarManager::kMass, 1.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow3") {
    cut->addCut(VarManager::kMass, 1.6, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow4") {
    cut->addCut(VarManager::kMass, 1.7, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow5") {
    cut->addCut(VarManager::kMass, 1.8, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow6") {
    cut->addCut(VarManager::kMass, 1.85, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow7") {
    cut->addCut(VarManager::kMass, 1.9, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow8") {
    cut->addCut(VarManager::kMass, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow9") {
    cut->addCut(VarManager::kMass, 2.2, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow10") {
    cut->addCut(VarManager::kMass, 2.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow11") {
    cut->addCut(VarManager::kMass, 3.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow12") {
    cut->addCut(VarManager::kMass, 3.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMass1to2") {
    cut->addCut(VarManager::kMass, 1., 2.);
    return cut;
  }

  if (nameStr == "pairMassIMR") {
    cut->addCut(VarManager::kMass, 1.1, 2.7);
    return cut;
  }

  if (nameStr == "pairMass1_5to2_7") {
    cut->addCut(VarManager::kMass, 1.5, 2.7);
    return cut;
  }

  if (nameStr == "pairMass1_3to3_5") {
    cut->addCut(VarManager::kMass, 1.3, 3.5);
    return cut;
  }

  if (nameStr == "pairMass1_3") {
    cut->addCut(VarManager::kMass, 1.3, 1000.0);
    return cut;
  }

  if (nameStr == "pairMass1_5to3_5") {
    cut->addCut(VarManager::kMass, 1.5, 3.5);
    return cut;
  }

  if (nameStr == "pairD0") {
    cut->addCut(VarManager::kMass, 1.814, 1.914);
    return cut;
  }

  if (nameStr == "pairD0Topology") {
    std::vector<double> ptBoundsD0 = {0.0, 2.0, 4.0, 6.0, 10.0, 30.0};
    std::vector<double> cosThetaPBoundsD0 = {0.9995, 0.9990, 0.9985, 0.9980, 0.9975};

    for (size_t thetaPBin = 1; thetaPBin < ptBoundsD0.size(); ++thetaPBin) {
      cut->addCut(VarManager::kCosPointingAngle, cosThetaPBoundsD0[thetaPBin - 1], 1.0, false, VarManager::kPt, ptBoundsD0[thetaPBin - 1], ptBoundsD0[thetaPBin]);
    }

    std::vector<double> daughterPtLowerBounds = {0.2, 0.3, 0.5, 0.7, 1.0};
    std::vector<double> dcaProductMax = {-0.020, -0.015, -0.010, -0.005, -0.002};
    std::vector<double> dcaDaughtersMax = {0.0050, 0.0040, 0.0030, 0.0025, 0.0020};

    for (size_t pTBin = 1; pTBin < ptBoundsD0.size(); ++pTBin) {
      double daughterLowerBound = daughterPtLowerBounds[pTBin - 1];
      double dcaProdUpperBound = dcaProductMax[pTBin - 1];
      double dcaDaughters = dcaDaughtersMax[pTBin - 1];

      double pTLowerValue = ptBoundsD0[pTBin - 1];
      double pTUpperValue = ptBoundsD0[pTBin];

      cut->addCut(VarManager::kPt1, daughterLowerBound, 1000.0, false, VarManager::kPt, pTLowerValue, pTUpperValue);
      cut->addCut(VarManager::kPt2, daughterLowerBound, 1000.0, false, VarManager::kPt, pTLowerValue, pTUpperValue);
      cut->addCut(VarManager::kDCATrackProd, -1e10, dcaProdUpperBound, false, VarManager::kPt, pTLowerValue, pTUpperValue);
      cut->addCut(VarManager::kDCADaughters, 0.0f, dcaDaughters, false, VarManager::kPt, pTLowerValue, pTUpperValue);
    }

    cut->addCut(VarManager::kCosThetaHE, -0.8f, 0.8f);
    cut->addCut(VarManager::kRap, -2.0, 2.0);

    return cut;
  }

  if (nameStr == "pairLambdaCTopology") {
    std::vector<double> pTBoundsLambdaC = {0.0, 2.0, 4.0, 6.0, 10.0, 30.0};

    std::vector<double> cosThetaPBoundsLambdaC = {0.9950, 0.9970, 0.9980, 0.9985, 0.9990};
    std::vector<double> daughterPtLowerBounds = {0.15, 0.25, 0.35, 0.50, 0.70};
    std::vector<double> lXYMinLambdaC = {0.0010, 0.0015, 0.0020, 0.0025, 0.0030};

    for (size_t pTBin = 1; pTBin < pTBoundsLambdaC.size(); ++pTBin) {
      double lowerBound = pTBoundsLambdaC[pTBin - 1];
      double upperBound = pTBoundsLambdaC[pTBin];

      double thetaP = cosThetaPBoundsLambdaC[pTBin - 1];
      double daughterPt = daughterPtLowerBounds[pTBin - 1];
      double lxyMin = lXYMinLambdaC[pTBin - 1];

      cut->addCut(VarManager::kCosPointingAngle, thetaP, 1.0, false, VarManager::kPt, lowerBound, upperBound);

      cut->addCut(VarManager::kPt1, daughterPt, 1000.0, false, VarManager::kPt, lowerBound, upperBound);
      cut->addCut(VarManager::kPt2, daughterPt, 1000.0, false, VarManager::kPt, lowerBound, upperBound);
      cut->addCut(VarManager::kPt3, daughterPt, 1000.0, false, VarManager::kPt, lowerBound, upperBound);

      cut->addCut(VarManager::kVertexingLxy, lxyMin, 100.0, false, VarManager::kPt, lowerBound, upperBound);
    }

    cut->addCut(VarManager::kRap, -2.0, 2.0);

    return cut;
  }

  if (nameStr == "pairJpsi") {
    cut->addCut(VarManager::kMass, 2.8, 3.3);
    return cut;
  }

  if (nameStr == "pairJpsi2") {
    cut->addCut(VarManager::kMass, 2.72, 3.2);
    return cut;
  }

  if (nameStr == "pairJpsi3") {
    cut->addCut(VarManager::kMass, 2.92, 3.14);
    return cut;
  }

  if (nameStr == "pairPsi2S") {
    cut->addCut(VarManager::kMass, 3.4, 3.9);
    return cut;
  }

  if (nameStr == "pairUpsilon") {
    cut->addCut(VarManager::kMass, 8.0, 11.0);
    return cut;
  }

  if (nameStr == "pairX3872") {
    cut->addCut(VarManager::kQ, 0.0, 0.3);
    return cut;
  }

  if (nameStr == "pairX3872_2") {
    cut->addCut(VarManager::kQuadDefaultDileptonMass, 3.0, 5.0);
    cut->addCut(VarManager::kQ, 0.0, 0.5);
    cut->addCut(VarManager::kDeltaR, 0.0, 5.0);
    cut->addCut(VarManager::kQuadPt, 5.0, 40.0);
    return cut;
  }

  if (nameStr == "pairX3872_3") {
    cut->addCut(VarManager::kQuadDefaultDileptonMass, 3.0, 5.0);
    cut->addCut(VarManager::kQ, 0.0, 0.5);
    cut->addCut(VarManager::kDeltaR, 0.0, 5.0);
    cut->addCut(VarManager::kQuadPt, 0.0, 1000.0);
    return cut;
  }

  std::shared_ptr<TF1> fDeltaR2High = std::make_shared<TF1>("fDeltaR2High", "[0]+[1]*x", 0.0, 5.0);
  fDeltaR2High->SetParameters(3.0, -0.4);
  if (nameStr == "pairX3872_minitree") {
    cut->addCut(VarManager::kRap, -0.8, 0.8);
    cut->addCut(VarManager::kQuadDefaultDileptonMass, 3.0, 5.0);
    cut->addCut(VarManager::kDeltaR2, 0.0, 1.0, false, VarManager::kPt, 5.0, 1000.0);
    cut->addCut(VarManager::kDeltaR2, 0.0, fDeltaR2High, false, VarManager::kPt, 0.0, 5.0);
    cut->addCut(VarManager::kQuadPt, 0.0, 1000.0);
    cut->addCut(VarManager::kVertexingProcCode, 0.5, 2.5);
    cut->addCut(VarManager::kVertexingQuadProcCode, 0.5, 2.5);
    return cut;
  }

  if (nameStr == "pairPtLow1") {
    cut->addCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairPtLow2") {
    cut->addCut(VarManager::kPt, 5.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairPtLow3") {
    cut->addCut(VarManager::kPt, 0, 0.15);
    return cut;
  }

  if (nameStr == "pairPtLow4") {
    cut->addCut(VarManager::kPt, 0, 10.0);
    return cut;
  }

  if (nameStr == "pairPtLow5") {
    cut->addCut(VarManager::kPt, 0.8, 1000.0);
    return cut;
  }

  if (nameStr == "pairRapidityForward") {
    cut->addCut(VarManager::kRap, 2.5, 4.0);
    return cut;
  }

  if (nameStr == "pairDCA") {
    cut->addCut(VarManager::kQuadDCAabsXY, .0, .50);
    return cut;
  }

  if (nameStr == "singleDCA") {
    cut->addCut(VarManager::kTrackDCAsigXY, 0.0, 5.);
    return cut;
  }

  if (nameStr == "pairPhiV") {
    cut->addCut(VarManager::kPairPhiv, 2., 3.2);
    return cut;
  }

  if (nameStr == "excludePairPhiV") {
    cut->addCut(VarManager::kPairPhiv, 2., 3.2, true);
    return cut;
  }

  if (nameStr == "pairLowMass") {
    cut->addCut(VarManager::kMass, 0., 0.1);
    return cut;
  }

  if (nameStr == "excludePairLowMass") {
    cut->addCut(VarManager::kMass, 0., 0.1, true);
    return cut;
  }

  if (nameStr == "pairLxyzProjected3sigma") {
    cut->addCut(VarManager::kVertexingLxyzProjected, 0.015, 10.);
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1") {
    cut->addCut(VarManager::kVertexingTauxyzProjected, 0.0005, 10.);
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1sigma") {
    cut->addCut(VarManager::kVertexingTauxyzProjected, 0.003, 10.);
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaLambdacCand") {
    std::shared_ptr<TF1> f1minLxyProjected = std::make_shared<TF1>("f1minLxyProjected", "[0]+[1]*x", 0., 20.);
    f1minLxyProjected->SetParameters(0.0065, -0.00023);
    cut->addCut(VarManager::kVertexingLxyProjected, f1minLxyProjected, 1., false, VarManager::kPt, 0., 20.);
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaDplusCand") {
    cut->addCut(VarManager::kVertexingLxyProjected, 0.009, 10.);
    return cut;
  }

  if (nameStr == "pairCosPointingPos") {
    cut->addCut(VarManager::kCosPointingAngle, 0.9, 1000.);
    return cut;
  }

  if (nameStr == "pairCosPointingNeg90") {
    cut->addCut(VarManager::kCosPointingAngle, -1000., -0.9);
    return cut;
  }

  if (nameStr == "pairCosPointingNeg85") {
    cut->addCut(VarManager::kCosPointingAngle, -1000., -0.85);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and pair selection in order or optimize the trigger
  // trigger selection cuts

  if (nameStr == "electronStandardQualityTriggerTest") {
    cut->addCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->addCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->addCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->addCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->addCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->addCut(VarManager::kTPCncls, 50.0, 161.);
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts") {
    cut->addCut(VarManager::kEta, -4.5, -2.0);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->addCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts") {
    cut->addCut(VarManager::kEta, -4.5, -2.0);
    cut->addCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->addCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->addCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->addCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->addCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->addCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest1") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest2") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest3") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest4") {
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->addCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->addCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->addCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  //---------------------------------------------------------------
  // ALICE 3 Analysis Cuts

  if (nameStr == "alice3KineSkim") {
    cut->addCut(VarManager::kPt, 0.05, 1000.0);
    cut->addCut(VarManager::kEta, -4.0, 4.0);
    return cut;
  }

  if (nameStr == "alice3StandardKine") {
    cut->addCut(VarManager::kPt, 0.1, 1000.0);
    cut->addCut(VarManager::kEta, -2.5, 2.5); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineTOFAcceptance") {
    cut->addCut(VarManager::kPt, 0.1, 1000.0);
    cut->addCut(VarManager::kEta, -2.0, 2.0); // TOF acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineRICHAcceptance") {
    cut->addCut(VarManager::kPt, 0.1, 1000.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8); // RICH acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3StandardKineD0") {
    cut->addCut(VarManager::kPt, 0.7, 1000.0);
    cut->addCut(VarManager::kEta, -2.5, 2.5); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3TOFEta") {
    cut->addCut(VarManager::kEta, -2.0, 2.0); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3RICHEta") {
    cut->addCut(VarManager::kEta, -0.8, 0.8); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineTOFAcceptanceD0") {
    cut->addCut(VarManager::kPt, 0.7, 1000.0);
    cut->addCut(VarManager::kEta, -2.0, 2.0); // TOF acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineRICHAcceptanceD0") {
    cut->addCut(VarManager::kPt, 0.7, 1000.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8); // TOF acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3LambdaCKine") {
    cut->addCut(VarManager::kPt, 0.2, 1000.0);
    cut->addCut(VarManager::kEta, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "alice3JpsiKine") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -2.5, 2.5);
    return cut;
  }

  if (nameStr == "alice3JpsiKineTOFAcceptance") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "alice3JpsiKineRICHAcceptance") {
    cut->addCut(VarManager::kPt, 1.0, 1000.0);
    cut->addCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "alice3TrackQuality") {
    cut->addCut(VarManager::kIsReconstructed, 0.5, 1.5);
    cut->addCut(VarManager::kNSiliconHits, 5.0, 12.0);
    return cut;
  }

  if (nameStr == "alice3TrackQualityHF") {
    cut->addCut(VarManager::kIsReconstructed, 0.5, 1.5);
    cut->addCut(VarManager::kNSiliconHits, 5.0, 12.0);
    cut->addCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->addCut(VarManager::kTrackDCAz, -0.1, 0.1);
    return cut;
  }

  if (nameStr == "alice3OTPIDEl") {
    cut->addCut(VarManager::kOTnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDPi") {
    cut->addCut(VarManager::kOTnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDKa") {
    cut->addCut(VarManager::kOTnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDPr") {
    cut->addCut(VarManager::kOTnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDEl") {
    cut->addCut(VarManager::kInnerTOFnSigmaEl, -3.0, 3.0, false, VarManager::kP, 0.0, 1.5);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPi") {
    cut->addCut(VarManager::kInnerTOFnSigmaPi, -3.0, 3.0, false, VarManager::kP, 0.0, 1.5);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDKa") {
    cut->addCut(VarManager::kInnerTOFnSigmaKa, -3.0, 3.0, false, VarManager::kP, 0.0, 1.5);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPr") {
    cut->addCut(VarManager::kInnerTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 1.5);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDEl") {
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0, false, VarManager::kP, 0.0, 2.5);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPi") {
    cut->addCut(VarManager::kOuterTOFnSigmaPi, -3.0, 3.0, false, VarManager::kP, 0.0, 2.5);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDKa") {
    cut->addCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0, false, VarManager::kP, 0.0, 2.5);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPr") {
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 2.5);
    return cut;
  }

  if (nameStr == "alice3KaonTOFRICH") {
    cut->addCut(VarManager::kInnerTOFnSigmaKa, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.057, 2.2);
    cut->addCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.3, 2.2);
    cut->addCut(VarManager::kTOFEventTimeErr, 0.0, 30.0, false, VarManager::kP, 0.0, 2.2);

    cut->addCut(VarManager::kOuterTOFnSigmaPi, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPi, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPr, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKEl, 3.0, 100.0);

    cut->addCut(VarManager::kEta, -0.8, 0.8, false, VarManager::kP, 2.2, 1000.0);
    cut->addCut(VarManager::kRICHnSigmaKa, -3.0, 3.0, false, VarManager::kHasRICHSigKa, 0.5, 1.5, false, VarManager::kEta, -0.8, 0.8);

    cut->addCut(VarManager::kRICHnSigmaPi, -3.0, 3.0, true, VarManager::kRICHSeparationKPi, 3.0, 100.0, false, VarManager::kHasRICHSigPi, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaPr, -3.0, 3.0, true, VarManager::kRICHSeparationKPr, 3.0, 100.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaEl, -3.0, 3.0, true, VarManager::kRICHSeparationKEl, 3.0, 100.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3PionTOFRICH") {
    cut->addCut(VarManager::kInnerTOFnSigmaPi, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.057, 2.2);
    cut->addCut(VarManager::kOuterTOFnSigmaPi, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.3, 2.2);
    cut->addCut(VarManager::kTOFEventTimeErr, 0.0, 30.0, false, VarManager::kP, 0.0, 2.2);

    cut->addCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPi, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, true, VarManager::kOuterTOFSeparationPiPr, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0, true, VarManager::kOuterTOFSeparationPiEl, 3.0, 100.0);

    cut->addCut(VarManager::kEta, -0.8, 0.8, false, VarManager::kP, 2.2, 1000.0);
    cut->addCut(VarManager::kRICHnSigmaPi, -3.0, 3.0, false, VarManager::kHasRICHSigPi, 0.5, 1.5, false, VarManager::kEta, -0.8, 0.8);

    cut->addCut(VarManager::kRICHnSigmaKa, -3.0, 3.0, true, VarManager::kRICHSeparationKPi, 3.0, 100.0, false, VarManager::kHasRICHSigKa, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaPr, -3.0, 3.0, true, VarManager::kRICHSeparationPiPr, 3.0, 100.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaEl, -3.0, 3.0, true, VarManager::kRICHSeparationPiEl, 3.0, 100.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3ProtonTOFRICH") {
    cut->addCut(VarManager::kInnerTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.057, 2.2);
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.3, 2.2);
    cut->addCut(VarManager::kTOFEventTimeErr, 0.0, 30.0, false, VarManager::kP, 0.0, 2.2);

    cut->addCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPi, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, true, VarManager::kOuterTOFSeparationPiPr, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0, true, VarManager::kOuterTOFSeparationPiEl, 3.0, 100.0);

    cut->addCut(VarManager::kEta, -0.8, 0.8, false, VarManager::kP, 2.2, 1000.0);
    cut->addCut(VarManager::kRICHnSigmaPi, -3.0, 3.0, false, VarManager::kHasRICHSigPi, 0.5, 1.5, false, VarManager::kEta, -0.8, 0.8);

    cut->addCut(VarManager::kRICHnSigmaKa, -3.0, 3.0, true, VarManager::kRICHSeparationKPi, 3.0, 100.0, false, VarManager::kHasRICHSigKa, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaPr, -3.0, 3.0, true, VarManager::kRICHSeparationPiPr, 3.0, 100.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaEl, -3.0, 3.0, true, VarManager::kRICHSeparationPiEl, 3.0, 100.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5);

    return cut;
  }

  if (nameStr == "alice3ElectronTOFRICH") {
    cut->addCut(VarManager::kInnerTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.057, 2.2);
    cut->addCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0, false, VarManager::kP, 0.0, 2.2, false, VarManager::kPt, 0.3, 2.2);
    cut->addCut(VarManager::kTOFEventTimeErr, 0.0, 30.0, false, VarManager::kP, 0.0, 2.2);

    cut->addCut(VarManager::kOuterTOFnSigmaPi, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPi, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKPr, 3.0, 100.0);
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0, true, VarManager::kOuterTOFSeparationKEl, 3.0, 100.0);

    cut->addCut(VarManager::kEta, -0.8, 0.8, false, VarManager::kP, 2.2, 1000.0);
    cut->addCut(VarManager::kRICHnSigmaPr, -3.0, 3.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5, false, VarManager::kEta, -0.8, 0.8);

    cut->addCut(VarManager::kRICHnSigmaPi, -3.0, 3.0, true, VarManager::kRICHSeparationKPi, 3.0, 100.0, false, VarManager::kHasRICHSigPi, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaKa, -3.0, 3.0, true, VarManager::kRICHSeparationKPr, 3.0, 100.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5);
    cut->addCut(VarManager::kRICHnSigmaEl, -3.0, 3.0, true, VarManager::kRICHSeparationKEl, 3.0, 100.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3CharmoniumPID") {
    cut->addCut(VarManager::kOuterTOFnSigmaEl, -2.0, 3.0, false, VarManager::kP, 0.0, 1.2);
    cut->addCut(VarManager::kRICHnSigmaEl, -2.0, 3.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5, false, VarManager::kP, 0.7, 1000.0);
    return cut;
  }

  if (nameStr == "alice3RICHPIDEl") {
    cut->addCut(VarManager::kRICHnSigmaEl, -3.0, 3.0, false, VarManager::kP, 0.7, 1000.0, false, VarManager::kHasRICHSigEl, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDPi") {
    cut->addCut(VarManager::kRICHnSigmaPi, -3.0, 3.0, false, VarManager::kP, 0.57, 1000.0, false, VarManager::kHasRICHSigPi, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDKa") {
    cut->addCut(VarManager::kRICHnSigmaKa, -3.0, 3.0, false, VarManager::kP, 2.0, 1000.0, false, VarManager::kHasRICHSigKa, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDPr") {
    cut->addCut(VarManager::kRICHnSigmaPr, -3.0, 3.0, false, VarManager::kP, 3.8, 1000.0, false, VarManager::kHasRICHSigPr, 0.5, 1.5);
    return cut;
  }

  delete cut;
  LOGF(fatal, Form("Did not find cut %s", cutName));
  return nullptr;
}

//________________________________________________________________________________________________
std::vector<AnalysisCut*> o2::aod::dqcuts::getCutsFromJSON(const char* json)
{
  //
  // configure cuts using a json file
  //

  std::vector<AnalysisCut*> cuts;
  // AnalysisCut* cuts[100];
  LOG(info) << "========================================== interpreting JSON for analysis cuts";
  LOG(info) << "JSON string: " << json;

  //
  // Create a vector of AnalysisCuts from a JSON formatted string
  //   The JSON is expected to contain a list of objects, with each object containing the fields needed
  //    to define either an AnalysisCut or an AnalysisCompositeCut
  rapidjson::Document document;

  // Check that the json is parsed correctly
  rapidjson::ParseResult ok = document.Parse(json);
  if (!ok) {
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")"; // NOLINT(google-readability-casting)
    return cuts;
  }

  // The json is expected to contain a list of objects, each of which should provide the configuration of a cut
  for (rapidjson::Value::ConstMemberIterator it = document.MemberBegin(); it != document.MemberEnd(); it++) {

    const char* cutName = it->name.GetString();
    LOG(info) << "=================================================== Configuring cut " << cutName;
    const auto& cut = it->value;

    // Detect if this is an AnalysisCut or a composite cut
    bool isAnalysisCut = false;
    bool isAnalysisCompositeCut = false;
    if (cut.HasMember("type")) {
      TString typeStr = cut.FindMember("type")->value.GetString();
      if (typeStr.CompareTo("AnalysisCut") == 0) {
        LOG(debug) << "This is an AnalysisCut";
        isAnalysisCut = true;
      }
      if (typeStr.CompareTo("AnalysisCompositeCut") == 0) {
        LOG(debug) << "This is an AnalysisCompositeCut";
        isAnalysisCompositeCut = true;
      }
    }
    if (!(isAnalysisCut || isAnalysisCompositeCut)) {
      LOG(fatal) << "Member is neither an AnalysisCut or AnalysisCompositeCut";
      return cuts;
    }

    // Parse the object, construct the cut and add it to the return vector
    if (isAnalysisCut) {
      AnalysisCut* anaCut = parseJSONAnalysisCut(&cut, cutName);
      if (anaCut != nullptr) {
        cuts.push_back(anaCut);
      } else {
        LOG(fatal) << "Something went wrong in creating the AnalysisCut " << cutName;
        return cuts;
      }
    }
    if (isAnalysisCompositeCut) {
      AnalysisCompositeCut* anaCut = parseJSONAnalysisCompositeCut(&cut, cutName);
      if (anaCut != nullptr) {
        cuts.push_back(anaCut);
      } else {
        LOG(fatal) << "Something went wrong in creating the AnalysisCompositeCut " << cutName;
        return cuts;
      }
    }
  }

  return cuts;
}

//________________________________________________________________________________________________
template <typename T>
bool o2::aod::dqcuts::validateJSONAnalysisCut(T cut)
{
  //
  // Validate cut definition in JSON file
  //
  // The type field is compulsory
  if (!cut->HasMember("type")) {
    LOG(fatal) << "Missing type information";
    return false;
  }
  TString typeStr = cut->FindMember("type")->value.GetString();
  if (typeStr.CompareTo("AnalysisCut") != 0) {
    LOG(fatal) << "Type is not AnalysisCut";
    return false;
  }

  return true;
}

//________________________________________________________________________________________________
template <typename T>
bool o2::aod::dqcuts::validateJSONaddCut(T addcut, bool isSimple)
{
  //
  // Validate addCut definition in JSON file
  //

  // Check if this addCut is adding an analysis cut (if the mother is a composite cut)
  bool isAnalysisCut = false;
  if (addcut->HasMember("type")) {
    isAnalysisCut = true;
  }
  // check if this is adding a basic variable range cut
  bool isBasicCut = false;
  if (addcut->HasMember("var") && addcut->HasMember("cutLow") && addcut->HasMember("cutHigh")) {
    isBasicCut = true;
  }

  // if neither of the two option is true, then something is wrong
  if (!(isBasicCut || isAnalysisCut)) {
    LOG(fatal) << "This is neither adding an AnalysisCut nor a basic variable range cut";
    return false;
  }
  if (isSimple && isAnalysisCut) {
    LOG(fatal) << "One cannot call addCut with an AnalysisCut in this case";
    return false;
  }
  if (!isSimple && isAnalysisCut) {
    return true;
  }

  // check that the variable to cut on is a valid one (implemented in the VarManager)
  const char* var = addcut->FindMember("var")->value.GetString();
  if (VarManager::fgVarNamesMap.find(var) == VarManager::fgVarNamesMap.end()) {
    LOG(fatal) << "Bad cut variable (" << var << ") specified for this addCut";
    return false;
  }

  // check whether the specified cut low and high are numbers or functions
  bool cutLowIsNumber = addcut->FindMember("cutLow")->value.IsNumber();
  bool cutHighIsNumber = addcut->FindMember("cutHigh")->value.IsNumber();
  if (!cutLowIsNumber) {
    auto& cutLow = addcut->FindMember("cutLow")->value;
    if (!cutLow.HasMember("funcName") || !cutLow.HasMember("funcBody") ||
        !cutLow.HasMember("xLow") || !cutLow.HasMember("xHigh")) {
      LOG(fatal) << "Missing fields for the cutLow TF1 definition";
      return false;
    }
  }
  if (!cutHighIsNumber) {
    auto& cutHigh = addcut->FindMember("cutHigh")->value;
    if (!cutHigh.HasMember("funcName") || !cutHigh.HasMember("funcBody") ||
        !cutHigh.HasMember("xLow") || !cutHigh.HasMember("xHigh")) {
      LOG(fatal) << "Missing fields for the cutLow TF1 definition";
      return false;
    }
  }
  if (!cutHighIsNumber || !cutLowIsNumber) {
    if (!addcut->HasMember("dependentVar") || !addcut->HasMember("depCutLow") || !addcut->HasMember("depCutHigh")) {
      LOG(fatal) << "For cutLow or cutHigh as a TF1, the definition of the dependentVar and range are also required";
      return false;
    }
    const char* depVar = addcut->FindMember("dependentVar")->value.GetString();
    if (VarManager::fgVarNamesMap.find(depVar) == VarManager::fgVarNamesMap.end()) {
      LOG(fatal) << "Bad cut variable (" << depVar << ") specified for the dependentVar";
      return false;
    }
  }
  if (addcut->HasMember("dependentVar1") && cutLowIsNumber && cutHighIsNumber) {
    const char* depVar1 = addcut->FindMember("dependentVar1")->value.GetString();
    if (VarManager::fgVarNamesMap.find(depVar1) == VarManager::fgVarNamesMap.end()) {
      LOG(fatal) << "Bad cut variable (" << depVar1 << ") specified for the dependentVar";
      return false;
    }
    if (!addcut->HasMember("depCut2Low") || !addcut->HasMember("depCut2High")) {
      LOG(fatal) << "dependentVar2 specified, but not its range";
      return false;
    }
  }
  if (addcut->HasMember("dependentVar2")) {
    const char* depVar2 = addcut->FindMember("dependentVar2")->value.GetString();
    if (VarManager::fgVarNamesMap.find(depVar2) == VarManager::fgVarNamesMap.end()) {
      LOG(fatal) << "Bad cut variable (" << depVar2 << ") specified for the dependentVar2";
      return false;
    }
    if (!addcut->HasMember("depCut2Low") || !addcut->HasMember("depCut2High")) {
      LOG(fatal) << "dependentVar2 specified, but not its range";
      return false;
    }
  }

  return true;
}

//_______________________________________________________________________________________________
template <typename T>
AnalysisCut* o2::aod::dqcuts::parseJSONAnalysisCut(T cut, const char* cutName)
{

  // Parse the json object and build an AnalysisCut
  if (!validateJSONAnalysisCut(cut)) {
    LOG(fatal) << "AnalysisCut not properly defined in the JSON file. Skipping it";
    return nullptr;
  }

  // If the analysis cut has the field "library", its just loaded from the preexisting cuts in the library and return
  if (cut->HasMember("library")) {
    return getAnalysisCut(cut->FindMember("library")->value.GetString());
  }

  // construct the AnalysisCut object and add the addCuts
  auto* retCut = new AnalysisCut(cutName, cut->HasMember("library") ? cut->FindMember("title")->value.GetString() : "");

  // loop over all the members for this cut and configure the addCut objects
  for (rapidjson::Value::ConstMemberIterator it = cut->MemberBegin(); it != cut->MemberEnd(); it++) {

    TString itName = it->name.GetString();

    if (itName.Contains("AddCut")) {

      LOG(debug) << "Parsing " << itName.Data();
      const auto& addcut = it->value;
      if (!validateJSONaddCut(&addcut, true)) {
        LOG(fatal) << "AddCut statement not properly defined";
        return nullptr;
      }

      const char* var = addcut.FindMember("var")->value.GetString();
      LOG(info) << "var " << var;
      bool cutLowIsNumber = addcut.FindMember("cutLow")->value.IsNumber();
      LOG(info) << "cutLowIsNumber " << cutLowIsNumber;
      bool cutHighIsNumber = addcut.FindMember("cutHigh")->value.IsNumber();
      LOG(info) << "cutHighIsNumber " << cutHighIsNumber;

      bool exclude = (addcut.HasMember("exclude") ? addcut.FindMember("exclude")->value.GetBool() : false);
      LOG(info) << "exclude " << exclude;
      const char* dependentVar = (addcut.HasMember("dependentVar") ? addcut.FindMember("dependentVar")->value.GetString() : "kNothing");
      LOG(info) << "dependentVar " << dependentVar;
      double depCutLow = (addcut.HasMember("depCutLow") ? addcut.FindMember("depCutLow")->value.GetDouble() : 0.0);
      LOG(info) << "depCutLow " << depCutLow;
      double depCutHigh = (addcut.HasMember("depCutHigh") ? addcut.FindMember("depCutHigh")->value.GetDouble() : 10.0);
      LOG(info) << "depCutHigh " << depCutHigh;
      bool depCutExclude = (addcut.HasMember("depCutExclude") ? addcut.FindMember("depCutExclude")->value.GetBool() : false);
      LOG(info) << "depCutExclude " << depCutExclude;
      const char* dependentVar2 = (addcut.HasMember("dependentVar2") ? addcut.FindMember("dependentVar2")->value.GetString() : "kNothing");
      LOG(info) << "dependentVar2 " << dependentVar2;
      double depCut2Low = (addcut.HasMember("depCut2Low") ? addcut.FindMember("depCut2Low")->value.GetDouble() : 0.0);
      LOG(info) << "depCut2Low " << depCut2Low;
      double depCut2High = (addcut.HasMember("depCut2High") ? addcut.FindMember("depCut2High")->value.GetDouble() : 10.0);
      LOG(info) << "depCut2High " << depCut2High;
      bool depCut2Exclude = (addcut.HasMember("depCut2Exclude") ? addcut.FindMember("depCut2Exclude")->value.GetBool() : false);
      LOG(info) << "depCut2Exclude " << depCut2Exclude;

      std::shared_ptr<TF1> cutLowFunc = nullptr;
      std::shared_ptr<TF1> cutHighFunc = nullptr;
      double cutLowNumber = 0.0;
      double cutHighNumber = 0.0;
      if (cutLowIsNumber) {
        cutLowNumber = addcut.FindMember("cutLow")->value.GetDouble();
        LOG(info) << "cutLowNumber " << cutLowNumber;
      } else {
        auto& cutLow = addcut.FindMember("cutLow")->value;
        cutLowFunc = std::make_shared<TF1>(cutLow.FindMember("funcName")->value.GetString(), cutLow.FindMember("funcBody")->value.GetString(),
                                           cutLow.FindMember("xLow")->value.GetDouble(), cutLow.FindMember("xHigh")->value.GetDouble());
        LOG(info) << "cutLowFunc " << cutLow.FindMember("funcName")->value.GetString() << ", " << cutLow.FindMember("funcBody")->value.GetString()
                  << ", " << cutLow.FindMember("xLow")->value.GetDouble() << ", " << cutLow.FindMember("xHigh")->value.GetDouble();
      }
      if (cutHighIsNumber) {
        cutHighNumber = addcut.FindMember("cutHigh")->value.GetDouble();
        LOG(info) << "cutHighNumber " << cutHighNumber;
      } else {
        auto& cutHigh = addcut.FindMember("cutHigh")->value;
        cutHighFunc = std::make_shared<TF1>(cutHigh.FindMember("funcName")->value.GetString(), cutHigh.FindMember("funcBody")->value.GetString(),
                                            cutHigh.FindMember("xLow")->value.GetDouble(), cutHigh.FindMember("xHigh")->value.GetDouble());
        LOG(info) << "cutHighFunc " << cutHigh.FindMember("funcName")->value.GetString() << ", " << cutHigh.FindMember("funcBody")->value.GetString()
                  << ", " << cutHigh.FindMember("xLow")->value.GetDouble() << ", " << cutHigh.FindMember("xHigh")->value.GetDouble();
      }
      if (cutLowIsNumber) {
        if (cutHighIsNumber) {
          retCut->addCut(VarManager::fgVarNamesMap[var], cutLowNumber, cutHighNumber, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        } else {
          retCut->addCut(VarManager::fgVarNamesMap[var], cutLowNumber, cutHighFunc, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        }
      } else {
        if (cutHighIsNumber) {
          retCut->addCut(VarManager::fgVarNamesMap[var], cutLowFunc, cutHighNumber, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        } else {
          retCut->addCut(VarManager::fgVarNamesMap[var], cutLowFunc, cutHighFunc, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        }
      }
    }
  }

  return retCut;
}

//_______________________________________________________________________________________________
template <typename T>
bool o2::aod::dqcuts::validateJSONAnalysisCompositeCut(T cut)
{
  //
  // Validate composite cut definition in JSON file
  //
  if (!cut->HasMember("type")) {
    LOG(fatal) << "Missing type field";
    return false;
  }
  if (!(cut->HasMember("library") || cut->HasMember("useAND"))) {
    LOG(fatal) << "Either library or useAND fields are required in an AnalysisCompositeCut definition";
  }
  TString typeStr = cut->FindMember("type")->value.GetString();
  if (typeStr.CompareTo("AnalysisCompositeCut") != 0) {
    LOG(fatal) << "Type is not AnalysisCompositeCut";
    return false;
  }

  return true;
}

//_______________________________________________________________________________________________
template <typename T>
AnalysisCompositeCut* o2::aod::dqcuts::parseJSONAnalysisCompositeCut(T cut, const char* cutName)
{

  // Configure an AnalysisCompositeCut
  if (!validateJSONAnalysisCompositeCut(cut)) {
    LOG(fatal) << "Composite Cut not properly defined in the JSON file. Skipping it";
    return nullptr;
  }

  if (cut->HasMember("library")) {
    return getCompositeCut(cut->FindMember("library")->value.GetString());
  }

  auto* retCut = new AnalysisCompositeCut(cutName, cut->HasMember("library") ? cut->FindMember("title")->value.GetString() : "", cut->FindMember("useAND")->value.GetBool());

  // Loop to find addCut objects
  for (rapidjson::Value::ConstMemberIterator it = cut->MemberBegin(); it != cut->MemberEnd(); it++) {

    TString itName = it->name.GetString();

    if (itName.Contains("AddCut")) {

      LOG(debug) << "Parsing " << itName.Data();
      const auto& addcut = it->value;
      if (!validateJSONaddCut(&addcut, false)) {
        LOG(fatal) << "AddCut statement not properly defined";
        return nullptr;
      }

      // For an AnalysisCompositeCut, one can call addCut with either an AnalysisCut or another AnalysisCompositeCut
      bool isAnalysisCut = false;
      bool isAnalysisCompositeCut = false;
      TString typeStr = addcut.FindMember("type")->value.GetString();
      if (typeStr.CompareTo("AnalysisCut") == 0) {
        isAnalysisCut = true;
      }
      if (typeStr.CompareTo("AnalysisCompositeCut") == 0) {
        isAnalysisCompositeCut = true;
      }

      // Add an AnalysisCut
      if (isAnalysisCut) {
        AnalysisCut* cutMember = parseJSONAnalysisCut(&addcut, itName.Data());
        if (cutMember != nullptr) {
          retCut->addCut(cutMember);
        } else {
          return nullptr;
        }
      }
      // Add an AnalysisCompositeCut
      if (isAnalysisCompositeCut) {
        AnalysisCompositeCut* cutMember = parseJSONAnalysisCompositeCut(&addcut, itName.Data());
        if (cutMember != nullptr) {
          retCut->addCut(cutMember);
        } else {
          return nullptr;
        }
      }
    }
  }

  return retCut;
}

//________________________________________________________________________________________________
o2::aod::dqmlcuts::BdtScoreConfig o2::aod::dqmlcuts::getBdtScoreCutsAndConfigFromJSON(const char* json)
{
  LOG(info) << "========================================== interpreting JSON for ML analysis cuts";
  if (!json) {
    LOG(fatal) << "JSON config string is null!";
    return {};
  }
  LOG(info) << "JSON string: " << json;

  rapidjson::Document document;

  // Check that the json is parsed correctly
  rapidjson::ParseResult ok = document.Parse(json);
  if (!ok) {
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")"; // NOLINT(google-readability-casting)
    return {};
  }

  for (auto it = document.MemberBegin(); it != document.MemberEnd(); ++it) {
    const auto& obj = it->value;

    // Classification type
    if (!obj.HasMember("type")) {
      LOG(fatal) << "Missing type (Binary/MultiClass)";
      return {};
    }
    TString typeStr = obj["type"].GetString();
    if (typeStr != "Binary" && typeStr != "MultiClass") {
      LOG(fatal) << "Unsupported classification type: " << typeStr;
      return {};
    }

    // Input features
    if (!obj.HasMember("inputFeatures") || !obj["inputFeatures"].IsArray()) {
      LOG(fatal) << "Missing inputFeatures member or array";
      return {};
    }
    std::vector<std::string> namesInputFeatures;
    for (const auto& feature : obj["inputFeatures"].GetArray()) {
      namesInputFeatures.emplace_back(feature.GetString());
      LOG(debug) << "Input features: " << feature.GetString();
    }

    // Model files
    if (!obj.HasMember("modelFiles") || !obj["modelFiles"].IsArray()) {
      LOG(fatal) << "Missing modelFiles member or array";
      return {};
    }
    std::vector<std::string> onnxFileNames;
    for (const auto& model : obj["modelFiles"].GetArray()) {
      onnxFileNames.emplace_back(model.GetString());
      LOG(debug) << "Model Files: " << model.GetString() << " ";
    }

    // Centrality estimation type
    if (!obj.HasMember("cent") || !obj["cent"].IsString()) {
      LOG(fatal) << "Missing cent member";
      return {};
    }
    std::string cent = obj["cent"].GetString();
    LOG(debug) << "Centrality type: " << cent;
    if (cent != "kCentFT0C" && cent != "kCentFT0A" && cent != "kCentFT0M") {
      LOG(fatal) << "Unsupported centrality type: " << cent;
      return {};
    }

    // Cut storage
    std::vector<std::pair<double, double>> centBins;
    std::vector<std::pair<double, double>> ptBins;
    std::vector<std::vector<double>> cutsMl;
    std::vector<int> cutDirs;
    std::vector<std::string> labelsFlatBin;
    bool cutDirsFilled = false;

    for (auto centMember = obj.MemberBegin(); centMember != obj.MemberEnd(); ++centMember) {
      TString centKey = centMember->name.GetString();
      if (!centKey.Contains("AddCentCut")) {
        continue;
      }

      const auto& centCut = centMember->value;

      // Centrality info
      if (!centCut.HasMember("centMin") || !centCut.HasMember("centMax")) {
        LOG(fatal) << "Missing centMin/centMax in " << centKey;
        return {};
      }
      double centMin = centCut["centMin"].GetDouble();
      double centMax = centCut["centMax"].GetDouble();

      for (auto ptMember = centCut.MemberBegin(); ptMember != centCut.MemberEnd(); ++ptMember) {
        TString ptKey = ptMember->name.GetString();
        if (!ptKey.Contains("AddPtCut")) {
          continue;
        }

        const auto& ptCut = ptMember->value;

        // Pt info
        if (!ptCut.HasMember("pTMin") || !ptCut.HasMember("pTMax")) {
          LOG(fatal) << "Missing pTMin/pTMax in " << ptKey;
          return {};
        }

        double ptMin = ptCut["pTMin"].GetDouble();
        double ptMax = ptCut["pTMax"].GetDouble();

        std::vector<double> binCuts;
        bool exclude = false;

        for (auto mlMember = ptCut.MemberBegin(); mlMember != ptCut.MemberEnd(); ++mlMember) {
          TString mlKey = mlMember->name.GetString();
          if (!mlKey.Contains("AddMLCut")) {
            continue;
          }

          const auto& mlcut = mlMember->value;

          if (!mlcut.HasMember("cut")) {
            LOG(fatal) << "Missing cut (score) in " << mlKey;
            return {};
          }

          double cutVal = mlcut["cut"].GetDouble();
          exclude = mlcut.HasMember("exclude") ? mlcut["exclude"].GetBool() : false;

          binCuts.push_back(cutVal);

          if (!cutDirsFilled) {
            cutDirs.push_back(exclude ? 0 : 1);
          }
        }

        cutDirsFilled = true;

        centBins.emplace_back(centMin, centMax);
        ptBins.emplace_back(ptMin, ptMax);
        cutsMl.push_back(binCuts);
        labelsFlatBin.push_back(Form("%s_cent%.0f_%.0f_pt%.1f_%.1f", cent.c_str(), centMin, centMax, ptMin, ptMax));
        LOG(info) << "Added cut for " << Form("%s_cent%.0f_%.0f_pt%.1f_%.1f", cent.c_str(), centMin, centMax, ptMin, ptMax) << " with cuts: [";
        std::string msg = "";
        for (size_t i = 0; i < binCuts.size(); ++i) {
          msg += std::to_string(binCuts[i]);
          if (i != binCuts.size() - 1) {
            msg += ", ";
          }
        }
        msg += "] and direction: ";
        msg += (exclude ? "CutGreater" : "CutSmaller");
        LOG(info) << msg;
      }
    }

    if (cutDirs.size() != cutsMl[0].size()) {
      LOG(fatal) << "Mismatch the cut size and direction size: cutsMl[0].size() = " << cutsMl[0].size()
                 << ", cutsMl[0].size() = " << cutDirs.size();
      return {};
    }

    std::vector<std::string> labelsClass;
    labelsClass.reserve(cutsMl[0].size());
    for (size_t j = 0; j < cutsMl[0].size(); ++j) {
      labelsClass.push_back(Form("score class %d", static_cast<int>(j)));
    }

    size_t nFlatBins = cutsMl.size();
    std::vector<double> binsMl(nFlatBins + 1);
    std::iota(binsMl.begin(), binsMl.end(), 0);

    // Binary
    if (typeStr == "Binary") {
      dqmlcuts::BinaryBdtScoreConfig binaryCfg;
      binaryCfg.inputFeatures = namesInputFeatures;
      binaryCfg.onnxFiles = onnxFileNames;
      binaryCfg.binsCent = centBins;
      binaryCfg.binsPt = ptBins;
      binaryCfg.binsMl = binsMl;
      binaryCfg.cutDirs = cutDirs;
      binaryCfg.centType = cent;
      binaryCfg.cutsMl = makeLabeledCutsMl(cutsMl, labelsFlatBin, labelsClass);

      return binaryCfg;
    }

    // MultiClass
    if (typeStr == "MultiClass") {
      dqmlcuts::MultiClassBdtScoreConfig multiCfg;
      multiCfg.inputFeatures = namesInputFeatures;
      multiCfg.onnxFiles = onnxFileNames;
      multiCfg.binsCent = centBins;
      multiCfg.binsPt = ptBins;
      multiCfg.binsMl = binsMl;
      multiCfg.cutDirs = cutDirs;
      multiCfg.centType = cent;
      multiCfg.cutsMl = makeLabeledCutsMl(cutsMl, labelsFlatBin, labelsClass);

      return multiCfg;
    }
  }

  return {};
}

o2::framework::LabeledArray<double> o2::aod::dqmlcuts::makeLabeledCutsMl(const std::vector<std::vector<double>>& cuts,
                                                                         const std::vector<std::string>& labelsflatBin,
                                                                         const std::vector<std::string>& labelsClass)
{
  const size_t nRows = cuts.size();
  const size_t nCols = cuts.empty() ? 0 : cuts[0].size();
  std::vector<double> flat;

  for (const auto& row : cuts) {
    flat.insert(flat.end(), row.begin(), row.end());
  }

  o2::framework::Array2D<double> arr(flat.data(), nRows, nCols);
  return o2::framework::LabeledArray<double>(arr, labelsflatBin, labelsClass);
}

int o2::aod::dqmlcuts::getMlBinIndex(double cent, double pt,
                                     const std::vector<std::pair<double, double>>& binsCent,
                                     const std::vector<std::pair<double, double>>& binsPt)
{
  LOG(debug) << "Searching for Ml bin index for cent: " << cent << ", pt: " << pt;
  for (size_t i = 0; i < binsCent.size(); ++i) {
    if (cent >= binsCent[i].first && cent < binsCent[i].second && pt >= binsPt[i].first && pt < binsPt[i].second) {
      LOG(debug) << " - Found at index: " << i;
      return static_cast<int>(i);
    }
  }
  return -1; // not found
}
