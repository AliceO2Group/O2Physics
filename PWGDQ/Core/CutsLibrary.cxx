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
/// \file CutsLibrary.cxxx
/// \brief Library with cuts to use

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

#include <RtypesCore.h>

#include <cstddef>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

AnalysisCompositeCut* o2::aod::dqcuts::GetCompositeCut(const char* cutName)
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
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));
    return cut;
  }
  if (nameStr == "Electron2023") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5_noCorr"));
    return cut;
  }
  if (nameStr == "Electron2025_1") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -3.0, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_2") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -3.0, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_3") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_4") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_5") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.25, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "Electron2025_4_ldong") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kP, 1.0, 1000.0);
    kineCut->AddCut(VarManager::kEta, -0.9, 0.9);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    qualityCuts->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    qualityCuts->AddCut(VarManager::kTPCncls, 70, 161.);
    qualityCuts->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);
    qualityCuts->AddCut(VarManager::kTrackDCAxy, -0.5, 0.5);

    auto* pidCuts = new AnalysisCut("pidCuts", "pid cuts");
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -2.5, 4.0, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.7, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPi, 2.7, 999, false, VarManager::kPin, 5.0, 1000.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0.0, 5.0);
    pidCuts->AddCut(VarManager::kTPCnSigmaPr, 2.7, 999, false, VarManager::kPin, 5.0, 1000.0);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(pidCuts);
    return cut;
  }

  if (nameStr == "LowMassElectron2023") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("lmee_pp_502TeV_TOFloose_pionrej"));
    return cut;
  }
  if (nameStr == "MuonLow2022") {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }
  if (nameStr == "MuonHigh2022") {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }
  if (nameStr == "MuonLow2023") {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }
  if (nameStr == "MuonHigh2023") {
    cut->AddCut(GetAnalysisCut("muonHighPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }
  if (nameStr == "ElectronForEMu") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaLoose"));
    return cut;
  }
  if (nameStr == "MuonForEMu") {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }
  // ///////////////////////////////////////////////
  //           End of Cuts for CEFP               //
  // ///////////////////////////////////////////////

  if (nameStr == "jpsiO2MCdebugCuts") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1"));
    return cut;
  }

  if (nameStr == "jpsiBenchmarkCuts") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityBenchmark"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts2_prefiltered1") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    cut->AddCut(GetAnalysisCut("notDalitzLeg1"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts3") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts4") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaLoose"));
    return cut;
  }

  if (nameStr == "electronSelection1_ionut") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }
  if (nameStr == "jpsiKineDcaQualitynoPID") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    return cut;
  }
  if (nameStr == "jpsiKineDcaQualitynoPID2") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    return cut;
  }
  if (nameStr == "electronSelection1_ionut_withTOFPID") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium_withLargeTOFPID"));
    return cut;
  }
  if (nameStr == "electronSelection1_idstoreh") { // same as electronSelection1_ionut, but with kIsSPDAny -> kIsITSibAny
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "electronSelection1pos_ionut") {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }
  if (nameStr == "electronSelection1neg_ionut") {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    return cut;
  }

  if (nameStr == "electronSelection2_ionut") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    cut->AddCut(GetAnalysisCut("insideTPCsector"));
    return cut;
  }
  if (nameStr == "electronSelection2pos_ionut") {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    cut->AddCut(GetAnalysisCut("insideTPCsector"));
    return cut;
  }
  if (nameStr == "electronSelection2neg_ionut") {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaMedium"));
    cut->AddCut(GetAnalysisCut("insideTPCsector"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts4_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts5") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryLoose"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts6") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryVeryLoose"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_noCorr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5_noCorr"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr_2") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine2"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts7_Corr_3") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine3"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts8") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts9") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftDown"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts10_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts10_Corr_Amb") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    cut->AddCut(GetAnalysisCut("ambiguousTrack")); // IsAmbiguous
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts11_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug3")); // cut on 1 ITS cluster
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts12") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryLoose"));     // with 3 sigma El TOF
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine4"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("pidCut_lowP_Corr"));

    auto* pidCut_highP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", kFALSE);
    pidCut_highP->AddCut(GetAnalysisCut("EleInclusion_highP_Corr"));
    pidCut_highP->AddCut(GetAnalysisCut("PionExclusion_highP_Corr"));
    cut->AddCut(pidCut_highP);
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine4"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("pidCut_lowP"));

    auto* pidCut_highP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", kFALSE);
    pidCut_highP->AddCut(GetAnalysisCut("EleInclusion_highP"));
    pidCut_highP->AddCut(GetAnalysisCut("PionExclusion_highP"));
    cut->AddCut(pidCut_highP);
    return cut;
  }
  if (nameStr == "jpsiO2MCdebugCuts_Pdependent2_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine4"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("pidCut_lowP_Corr"));

    auto* pidCut_highP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", kFALSE);
    pidCut_highP->AddCut(GetAnalysisCut("EleInclusion_highP2_Corr"));
    pidCut_highP->AddCut(GetAnalysisCut("PionExclusion_highP_Corr"));
    cut->AddCut(pidCut_highP);
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts_Pdependent2") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine4"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("pidCut_lowP"));

    auto* pidCut_highP = new AnalysisCompositeCut("pidCut_highP", "pidCut_highP", kFALSE);
    pidCut_highP->AddCut(GetAnalysisCut("EleInclusion_highP2"));
    pidCut_highP->AddCut(GetAnalysisCut("PionExclusion_highP"));
    cut->AddCut(pidCut_highP);
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts1") { // please do not remove or modify, this is used for the common Skimmed tree production, (Xiaozhi Bai)
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed"));
    cut->AddCut(GetAnalysisCut("electronPIDLooseSkimmed"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts1") { // please do not remove or modify, this is used for the common Skimmed tree production, (Xiaozhi Bai)
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed"));
    cut->AddCut(GetAnalysisCut("electronPIDLooseSkimmed"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts2") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed"));
    cut->AddCut(GetAnalysisCut("electronPIDLooseSkimmed2"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts3") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed2"));
    cut->AddCut(GetAnalysisCut("electronPIDLooseSkimmed2"));
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts4") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed2"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug9")); // loose cut
    return cut;
  }

  if (nameStr == "JpsiPWGSkimmedCuts5") {
    cut->AddCut(GetAnalysisCut("electronTrackQualitySkimmed3"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug8"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut") {
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine3"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut_posEta") {
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    cut->AddCut(GetAnalysisCut("jpsiPIDcalibKine_posEta"));
    return cut;
  }

  if (nameStr == "pidElectron_ionut_negEta") {
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    cut->AddCut(GetAnalysisCut("jpsiPIDcalibKine_negEta"));
    return cut;
  }

  if (nameStr == "pidPion_ionut") {
    cut->AddCut(GetAnalysisCut("pidcalib_pion"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine3"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts13_Corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA")); // with DCA cut
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts14") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed"));
    return cut;
  }

  if (nameStr == "jpsiO2MCdebugCuts14andDCA") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "emu_electronCuts") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tof") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed"));
    cut->AddCut(GetAnalysisCut("tof_electron_sigma_2"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tightTPC") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed_2"));
    return cut;
  }

  if (nameStr == "emu_electronCuts_tof_tightTPC") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaSkewed_2"));
    cut->AddCut(GetAnalysisCut("tof_electron_sigma_2"));
    return cut;
  }

  if (nameStr == "jpsiKineAndQuality") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "jpsiPID1") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1"));
    return cut;
  }

  if (nameStr == "jpsiPID2") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID2"));
    return cut;
  }

  if (nameStr == "jpsiPIDnsigma") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDCut1") {
    cut->AddCut(GetAnalysisCut("pionQualityCut1"));
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDCut2") {
    cut->AddCut(GetAnalysisCut("pionQualityCut2"));
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "PIDCalibElectron") {
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "PIDCalibPion") {
    cut->AddCut(GetAnalysisCut("pidcalib_pion"));
    return cut;
  }

  if (nameStr == "PIDCalibKaon") {
    cut->AddCut(GetAnalysisCut("pidcalib_kaon"));
    return cut;
  }

  if (nameStr == "PIDCalibProton") {
    cut->AddCut(GetAnalysisCut("pidcalib_proton"));
    return cut;
  }

  if (nameStr == "PIDCalib_basic") {
    cut->AddCut(GetAnalysisCut("pidbasic"));
    return cut;
  }

  if (nameStr == "PIDefficiency_wPID") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "PIDefficiency_woPID") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "highPtHadron") {
    cut->AddCut(GetAnalysisCut("highPtHadron"));
    return cut;
  }

  if (nameStr == "rho0Cuts") {
    cut->AddCut(GetAnalysisCut("rho0Kine"));
    cut->AddCut(GetAnalysisCut("pionQuality"));
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "rho0Kine") {
    cut->AddCut(GetAnalysisCut("rho0Kine"));
    return cut;
  }

  if (nameStr == "openEtaSel") {
    cut->AddCut(GetAnalysisCut("openEtaSel"));
    return cut;
  }

  if (nameStr == "hasTOF") {
    cut->AddCut(GetAnalysisCut("hasTOF"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts1") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("SPDany"));
    cut->AddCut(GetAnalysisCut("openEtaSel"));
    cut->AddCut(GetAnalysisCut("pionQuality"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts2") {
    cut->AddCut(GetAnalysisCut("muonLowPt3"));
    cut->AddCut(GetAnalysisCut("ITSiball"));
    cut->AddCut(GetAnalysisCut("openEtaSel"));
    cut->AddCut(GetAnalysisCut("pionQuality"));
    return cut;
  }

  if (nameStr == "singleGapTrackCuts3") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine2"));
    cut->AddCut(GetAnalysisCut("SPDany"));

    auto* cut_notpc = new AnalysisCompositeCut("NoTPC", "NoTPC", kTRUE);
    cut_notpc->AddCut(GetAnalysisCut("noTPC"));

    auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
    cut_tpcpid->AddCut(GetAnalysisCut("pionQuality"));

    auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
    cut_OR->AddCut(cut_notpc);
    cut_OR->AddCut(cut_tpcpid);
    cut->AddCut(cut_OR);
    return cut;
  }

  if (nameStr == "singleGapTrackCuts4") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine2"));
    cut->AddCut(GetAnalysisCut("ITSibany"));

    auto* cut_notpc = new AnalysisCompositeCut("NoTPC", "NoTPC", kTRUE);
    cut_notpc->AddCut(GetAnalysisCut("noTPC"));

    auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
    cut_tpcpid->AddCut(GetAnalysisCut("pionQuality"));

    auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
    cut_OR->AddCut(cut_notpc);
    cut_OR->AddCut(cut_tpcpid);
    cut->AddCut(cut_OR);
    return cut;
  }

  if (nameStr == "PIDCalib") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("pidcalib_ele"));
    return cut;
  }

  if (nameStr == "NoPID") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine2")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly2"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    return cut;
  }

  if (nameStr == "KineCutOnly") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "KineCutOnly2") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine2")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "KineCutOnly3") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine3")); // standard kine cuts usually are applied via Filter in the task
    return cut;
  }

  if (nameStr == "kaonPID") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  if (nameStr == "kaonPID2") {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma2"));
    return cut;
  }

  if (nameStr == "kaonPID3") {
    cut->AddCut(GetAnalysisCut("AssocKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID3_withDCA") { // same as kaonPID3 but with cut on DCA and SPDAny->ITSAny
    cut->AddCut(GetAnalysisCut("AssocKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug4"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID4") {
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "kaonPID5") {
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma"));
    return cut;
  }

  if (nameStr == "kaonPID6") {
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma700"));
    return cut;
  }

  if (nameStr == "kaonPIDTPCTOForTPC") {
    auto* cut_tpctof_nSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", kTRUE);
    cut_tpctof_nSigma->AddCut(GetAnalysisCut("hasTOF"));
    cut_tpctof_nSigma->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("noTOF"));
    cut_tpc_nSigma->AddCut(GetAnalysisCut("kaonPIDnsigma"));

    auto* cut_pid_OR = new AnalysisCompositeCut("kaon_nsigma", "kaon_nsigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpctof_nSigma);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "kaonPIDTPCTOForTPC700") {
    auto* cut_tpctof_nSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", kTRUE);
    cut_tpctof_nSigma->AddCut(GetAnalysisCut("hasTOF"));
    cut_tpctof_nSigma->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("noTOF"));
    cut_tpc_nSigma->AddCut(GetAnalysisCut("kaonPIDnsigma700"));

    auto* cut_pid_OR = new AnalysisCompositeCut("kaon_nsigma", "kaon_nsigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpctof_nSigma);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "kaonPosPID4") {
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "kaonPosPID4Pt05") {
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "kaonNegPID4") {
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "kaonNegPID4Pt05") {
    cut->AddCut(GetAnalysisCut("kaonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "pionPID") {
    cut->AddCut(GetAnalysisCut("pionPID_TPCnTOF"));
    return cut;
  }

  if (nameStr == "pionPID2") {
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    return cut;
  }

  if (nameStr == "pionPIDTPCTOForTPC") {
    auto* cut_tpctof_nSigma = new AnalysisCompositeCut("pid_TPCTOFnSigma", "pid_TPCTOFnSigma", kTRUE);
    cut_tpctof_nSigma->AddCut(GetAnalysisCut("pionPID_TPCnTOF"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("pionPIDnsigma"));

    auto* cut_pid_OR = new AnalysisCompositeCut("pion_nsigma", "pion_nsigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpctof_nSigma);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "pionPosPID") {
    cut->AddCut(GetAnalysisCut("pionPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "pionPosPID2") {
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "pionNegPID") {
    cut->AddCut(GetAnalysisCut("pionPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "pionNegPID2") {
    cut->AddCut(GetAnalysisCut("pionPIDnsigma"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "protonPosPID") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "protonPosPIDPt05") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "protonNegPID") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "protonNegPIDPt05") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "protonPIDPV") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF2"));
    cut->AddCut(GetAnalysisCut("protonPVcut"));
    return cut;
  }

  if (nameStr == "protonPIDPV2") {
    cut->AddCut(GetAnalysisCut("protonPID_TPCnTOF2"));
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCAz") {
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "posPrimaryTrack_DCAz") {
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negPrimaryTrack_DCAz") {
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posStandardPrimaryTrackDCA") {
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negStandardPrimaryTrackDCA") {
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posTrack") {
    cut->AddCut(GetAnalysisCut("posTrack"));
    return cut;
  }

  if (nameStr == "negTrack") {
    cut->AddCut(GetAnalysisCut("negTrack"));
    return cut;
  }

  if (nameStr == "posTrackKaonRej") {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("kaonRejNsigma"));
    return cut;
  }

  if (nameStr == "negTrackKaonRej") {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("kaonRejNsigma"));
    return cut;
  }

  if (nameStr == "pTLow05") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    return cut;
  }

  if (nameStr == "pTLow04") {
    cut->AddCut(GetAnalysisCut("pTLow04"));
    return cut;
  }

  if (nameStr == "pTLow03") {
    cut->AddCut(GetAnalysisCut("pTLow03"));
    return cut;
  }

  if (nameStr == "pTLow02") {
    cut->AddCut(GetAnalysisCut("pTLow02"));
    return cut;
  }

  if (nameStr == "pTLow05DCAzHigh03") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "pTLow04DCAzHigh03") {
    cut->AddCut(GetAnalysisCut("pTLow04"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  if (nameStr == "pTLow03DCAzHigh03") {
    cut->AddCut(GetAnalysisCut("pTLow03"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_DCAz"));
    return cut;
  }

  // NOTE Below there are several TPC pid cuts used for studies of the Run3 TPC post PID calib.
  if (nameStr == "Jpsi_TPCPost_calib_debug1") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug2") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_noITSCuts_debug2") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug3") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug3"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_debug4") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }
  if (nameStr == "Jpsi_TPCPost_calib_noITSCuts_debug4") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug6") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug2"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug6"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug7") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug2"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug7"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug8") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug5"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug8"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug9") {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug4"));
    cut->AddCut(GetAnalysisCut("electronPIDLooseSkimmed3"));
    return cut;
  }

  if (nameStr == "Jpsi_TPCPost_calib_debug10") {
    cut->AddCut(GetAnalysisCut("jpsiKineSkimmed"));
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug6"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug10"));
    return cut;
  }

  if (nameStr == "LMee_TPCPost_calib_debug1") {
    cut->AddCut(GetAnalysisCut("lmee_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("lmee_TPCPID_debug1"));
    return cut;
  }

  if (nameStr == "ITSalone_prefilter") {
    cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "ITSalonebAny_prefilter") {
    cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "TPCalone_prefilter") {
    cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "ITSTPC_prefilter") {
    cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  for (int iCut = 0; iCut < 10; iCut++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("jpsiEleSel%d_ionut", iCut)) {
      cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
      cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
      cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
      cut->AddCut(GetAnalysisCut(Form("pidJpsiEle%d_ionut", iCut)));
      return cut;
    }

    if (nameStr == Form("jpsiEleSelTight%d_ionut", iCut)) {
      cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
      cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
      cut->AddCut(GetAnalysisCut("trackQualityTight_ionut"));
      cut->AddCut(GetAnalysisCut(Form("pidJpsiEle%d_ionut", iCut)));
      return cut;
    }
  }

  // Magnus composite cuts -----------------------------------------------------------------------------------------------------------------

  auto* magnus_PID111 = new AnalysisCompositeCut("magnus_PID111", "");
  magnus_PID111->AddCut(GetAnalysisCut("pidJpsi_magnus_ele1"));
  magnus_PID111->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID111->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization111") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID111);
    return cut;
  }

  auto* magnus_PID211 = new AnalysisCompositeCut("magnus_PID211", "");
  magnus_PID211->AddCut(GetAnalysisCut("pidJpsi_magnus_ele2"));
  magnus_PID211->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID211->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization211") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID211);
    return cut;
  }

  auto* magnus_PID311 = new AnalysisCompositeCut("magnus_PID311", "");
  magnus_PID311->AddCut(GetAnalysisCut("pidJpsi_magnus_ele3"));
  magnus_PID311->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID311->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization311") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID311);
    return cut;
  }

  auto* magnus_PID121 = new AnalysisCompositeCut("magnus_PID121", "");
  magnus_PID121->AddCut(GetAnalysisCut("pidJpsi_magnus_ele1"));
  magnus_PID121->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID121->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization121") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID121);
    return cut;
  }

  auto* magnus_PID112 = new AnalysisCompositeCut("magnus_PID112", "");
  magnus_PID112->AddCut(GetAnalysisCut("pidJpsi_magnus_ele1"));
  magnus_PID112->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID112->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization112") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID112);
    return cut;
  }

  auto* magnus_PID122 = new AnalysisCompositeCut("magnus_PID122", "");
  magnus_PID122->AddCut(GetAnalysisCut("pidJpsi_magnus_ele1"));
  magnus_PID122->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID122->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization122") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID122);
    return cut;
  }

  auto* magnus_PID222 = new AnalysisCompositeCut("magnus_PID222", "");
  magnus_PID222->AddCut(GetAnalysisCut("pidJpsi_magnus_ele2"));
  magnus_PID222->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID222->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization222") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID222);
    return cut;
  }

  auto* magnus_PID212 = new AnalysisCompositeCut("magnus_PID212", "");
  magnus_PID212->AddCut(GetAnalysisCut("pidJpsi_magnus_ele2"));
  magnus_PID212->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID212->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization212") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID212);
    return cut;
  }

  auto* magnus_PID221 = new AnalysisCompositeCut("magnus_PID221", "");
  magnus_PID221->AddCut(GetAnalysisCut("pidJpsi_magnus_ele2"));
  magnus_PID221->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID221->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization221") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID221);
    return cut;
  }

  auto* magnus_PID321 = new AnalysisCompositeCut("magnus_PID321", "");
  magnus_PID321->AddCut(GetAnalysisCut("pidJpsi_magnus_ele3"));
  magnus_PID321->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID321->AddCut(GetAnalysisCut("pidJpsi_magnus_prot1"));
  if (nameStr == "MagnussOptimization321") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID321);
    return cut;
  }

  auto* magnus_PID312 = new AnalysisCompositeCut("magnus_PID312", "");
  magnus_PID312->AddCut(GetAnalysisCut("pidJpsi_magnus_ele3"));
  magnus_PID312->AddCut(GetAnalysisCut("pidJpsi_magnus_pion1"));
  magnus_PID312->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization312") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID312);
    return cut;
  }

  auto* magnus_PID322 = new AnalysisCompositeCut("magnus_PID322", "");
  magnus_PID322->AddCut(GetAnalysisCut("pidJpsi_magnus_ele1"));
  magnus_PID322->AddCut(GetAnalysisCut("pidJpsi_magnus_pion2"));
  magnus_PID322->AddCut(GetAnalysisCut("pidJpsi_magnus_prot2"));
  if (nameStr == "MagnussOptimization322") {
    cut->AddCut(GetAnalysisCut("kineJpsiEle_ionut"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("trackQuality_ionut"));
    cut->AddCut(magnus_PID322);
    return cut;
  }
  //-------------------------------------------------------------------------------------------------------

  //---------------------------------------------------------------
  // Cuts for the selection of legs from dalitz decay
  //
  if (nameStr == "DalitzCut1") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (nameStr == "DalitzCut2") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (nameStr == "DalitzCut2_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut3") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (nameStr == "DalitzCut3_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut1SPDfirst") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (nameStr == "DalitzCut1SPDfirst_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut2SPDfirst") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (nameStr == "DalitzCut2SPDfirst_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (nameStr == "DalitzCut3SPDfirst") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (nameStr == "DalitzCut3SPDfirst_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  if (nameStr == "Dalitz_WithTOF_SPDfirst") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose"));

    auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma"));

    auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "Dalitz_WithTOF_SPDfirst_Corr") {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));

    auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_Corr"));

    auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("dalitzSelected%d", i)) {
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }
  }

  if (nameStr == "electronPrimaryTag0") {
    // with tight 3 sigma DCA cut for selecting primary electrons
    cut->AddCut(GetAnalysisCut("electronPID_TPCnsigma_loose")); // 3 sigma inclusion, 3sigma rejection
    cut->AddCut(GetAnalysisCut("electronPrimary_dca3sigma"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryTag1") {
    // with 7 sigma DCA cut for selecting primary electrons
    cut->AddCut(GetAnalysisCut("electronPID_TPCnsigma_loose")); // 3 sigma inclusion, 3sigma rejection
    cut->AddCut(GetAnalysisCut("electronPrimary_dca7sigma"));
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_TPC") {
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_ITS") {
    cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA"));
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "electronPrimaryProbe_ITSTPC") {
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA"));
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    return cut;
  }

  if (nameStr == "jpsiPIDworseRes") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDworseRes"));
    return cut;
  }

  if (nameStr == "jpsiPIDshift") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPIDshift"));
    return cut;
  }

  if (nameStr == "jpsiPID1shiftUp") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftUp"));
    return cut;
  }

  if (nameStr == "jpsiPID1shiftDown") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQuality"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    cut->AddCut(GetAnalysisCut("electronPID1shiftDown"));
    return cut;
  }
  // -------------------------------------------------------------------------------------------------
  //
  // Q vector contributor cut
  //
  if (nameStr == "selTPCCentral") {
    auto* kineCut = new AnalysisCut("kineCut", "kine cut");
    kineCut->AddCut(VarManager::kEta, -0.8, 0.8);
    kineCut->AddCut(VarManager::kPt, 0.15, 5);

    auto* qualityCuts = new AnalysisCut("qualityCuts", "quality cuts");
    qualityCuts->AddCut(VarManager::kTPCchi2, 0., 4.);
    qualityCuts->AddCut(VarManager::kTPCnCRoverFindCls, 0.8, 1.);
    qualityCuts->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    qualityCuts->AddCut(VarManager::kITSchi2, 0., 36.);

    auto* dcaCuts = new AnalysisCut("dcaCuts", "DCA cuts");
    std::shared_ptr<TF1> f1dcaxyHigh = std::make_shared<TF1>("f1dcaxy", "[0]+[1]/pow(x,[2])", 0., 10.);
    f1dcaxyHigh->SetParameters(0.0105, 0.035, 1.1);
    std::shared_ptr<TF1> f1dcaxyLow = std::make_shared<TF1>("f1dcaxy_low", "[0]+[1]/pow(x,[2])", 0., 10.);
    f1dcaxyLow->SetParameters(-0.0105, -0.035, 1.1);
    dcaCuts->AddCut(VarManager::kTrackDCAxy, f1dcaxyLow, f1dcaxyHigh);
    dcaCuts->AddCut(VarManager::kTrackDCAz, -2., 2.);

    cut->AddCut(kineCut);
    cut->AddCut(qualityCuts);
    cut->AddCut(dcaCuts);
  }

  // -------------------------------------------------------------------------------------------------
  //
  // LMee cuts
  // List of cuts used for low mass dielectron analyses
  //
  // Skimming cuts:
  if (nameStr == "lmee_skimming") {
    cut->AddCut(GetAnalysisCut("lmee_skimming_cuts"));
    return cut;
  }

  // LMee Run2 PID cuts

  if (nameStr == "lmeePID_TPChadrejTOFrec") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    auto* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));

    auto* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron"));

    auto* cut_pid_OR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_hadrej);
    cut_pid_OR->AddCut(cut_tof_rec);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "lmeePID_TPChadrej") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    auto* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));
    cut->AddCut(cut_tpc_hadrej);
    return cut;
  }

  if (nameStr == "lmee_eNSigmaRun2") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma"));

    auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma"));

    auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "lmeePID_TOFrec") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));

    auto* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron"));

    cut->AddCut(cut_tof_rec);
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun2") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun2_lowPt") {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun2_lowPt") {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun2") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  // LMee Run3 PID cuts

  if (nameStr == "lmeePID_TPChadrejTOFrecRun3") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

    auto* cut_tpc_hadrej = new AnalysisCompositeCut("pid_TPChadrej", "pid_TPChadrej", kTRUE);
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_muon_band_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_pion_rejection_highp"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_kaon_rejection"));
    cut_tpc_hadrej->AddCut(GetAnalysisCut("tpc_proton_rejection"));

    auto* cut_tof_rec = new AnalysisCompositeCut("pid_tof_rec", "pid_tof_rec", kTRUE);
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_electron"));
    cut_tof_rec->AddCut(GetAnalysisCut("tof_electron_loose"));
    cut_tof_rec->AddCut(GetAnalysisCut("tpc_pion_rejection_highp"));

    auto* cut_pid_OR = new AnalysisCompositeCut("pid_TPChadrejTOFrec", "pid_TPChadrejTOFrec", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_hadrej);
    cut_pid_OR->AddCut(cut_tof_rec);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (nameStr == "lmee_Run3_TPCelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_TOFelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("tof_electron_sigma"));
    return cut;
  }

  if (nameStr == "lmee_Run3_posTrack_TPCelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_posTrack_TOFelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("tof_electron_sigma"));
    return cut;
  }

  if (nameStr == "lmee_Run3_negTrack_TPCelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("lmee_pp_502TeV_TPCPbPbnopkrej"));
    return cut;
  }

  if (nameStr == "lmee_Run3_negTrack_TOFelectron") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
    cut->AddCut(GetAnalysisCut("tof_electron_sigma"));
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
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee%s_PbPb_selection_pt04", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee%s_TrackCuts_Resol", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("openEtaSel")); // No pt cut and wider eta cut to produce resolution maps
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 4 cuts to separate pos & neg tracks in pos & neg eta range
    if (nameStr == Form("lmee_posTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("posEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("posEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_posTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("negEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("negEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 2 cuts to separate pos & neg tracks
    if (nameStr == Form("lmee_posTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("etaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    if (nameStr == Form("lmee_negTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("etaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      return cut;
    }

    // 4 cuts to separate pos & neg tracks in pos & neg eta range low B field
    if (nameStr == Form("lmee_lowB_posTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("posEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_posEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("posEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_posTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("negEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_negEta_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("negEtaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    // 2 cuts to separate pos & neg tracks in low B field
    if (nameStr == Form("lmee_lowB_posTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("posTrack"));
      cut->AddCut(GetAnalysisCut("etaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
      return cut;
    }

    if (nameStr == Form("lmee_lowB_negTrack_selection%s", vecTypetrack.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("negTrack"));
      cut->AddCut(GetAnalysisCut("etaSel"));
      cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrack.at(icase).Data())));
      cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
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
      cut->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITSTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      cut->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITS_ifTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_notpc = new AnalysisCompositeCut("NoTPC", "NoTPC", kTRUE);
      cut_notpc->AddCut(GetAnalysisCut("noTPC"));

      auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
      cut_tpcpid->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut_tpcpid->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
      cut_OR->AddCut(cut_notpc);
      cut_OR->AddCut(cut_tpcpid);
      cut->AddCut(cut_OR);
      return cut;
    }

    if (nameStr == Form("ITS_ifTPCStandard_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_notpcstandard = new AnalysisCompositeCut("NoTPCstandard", "NoTPCstandard", kTRUE);
      cut_notpcstandard->AddCut(GetAnalysisCut("NoelectronStandardQualityTPCOnly"));

      auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
      cut_tpcpid->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut_tpcpid->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
      cut_OR->AddCut(cut_notpcstandard);
      cut_OR->AddCut(cut_tpcpid);
      cut->AddCut(cut_OR);
      return cut;
    }

    if (nameStr == Form("ITSTPCbAny_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      cut->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("ITSbAny_ifTPC_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_notpc = new AnalysisCompositeCut("NoTPC", "NoTPC", kTRUE);
      cut_notpc->AddCut(GetAnalysisCut("noTPC"));

      auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
      cut_tpcpid->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut_tpcpid->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
      cut_OR->AddCut(cut_notpc);
      cut_OR->AddCut(cut_tpcpid);
      cut->AddCut(cut_OR);
      return cut;
    }

    if (nameStr == Form("ITSbAny_ifTPCStandard_TPCPID%s_prefilter", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeePrefilterKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualitybAnyITSOnly"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_notpcstandard = new AnalysisCompositeCut("NoTPCstandard", "NoTPCstandard", kTRUE);
      cut_notpcstandard->AddCut(GetAnalysisCut("NoelectronStandardQualityTPCOnly"));

      auto* cut_tpcpid = new AnalysisCompositeCut("pid_TPC", "pid_TPC", kTRUE);
      cut_tpcpid->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut_tpcpid->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));

      auto* cut_OR = new AnalysisCompositeCut("OR", "OR", kFALSE);
      cut_OR->AddCut(cut_notpcstandard);
      cut_OR->AddCut(cut_tpcpid);
      cut->AddCut(cut_OR);
      return cut;
    }

    for (unsigned int i = 0; i < 30; i++) { // o2-linter: disable=magic-number (number of cuts)
      if (nameStr == Form("ElSelCutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeCutVarTrackCuts%i", i)));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }
    }

    for (size_t jcase = 0; jcase < vecTypetrackWithPID.size(); jcase++) {
      // All previous cut with TightGlobalTrackRun3
      if (nameStr == Form("ITSTPC%s_TPCPIDalone%s_PbPb", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPIDOnly%s", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_loose", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_loose", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_loose", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongHadRej", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TOFreq", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TPC_TOFnsigma%s", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_Resol", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("openEtaSel")); // No pt cut and wider eta cut to produce resolution maps
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_tightNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_tightNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_eNSigmaRun3%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz")); // to reject looper using DCAz

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_eNSigmaRun3%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz")); // to reject looper using DCAz

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3%s_tightNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      // 8 cuts for QC
      if (nameStr == Form("lmee%s_NSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt02Sel"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt02Sel"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt04Sel"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_NSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt04Sel"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      // 6 cuts for QC
      if (nameStr == Form("lmee%s_posTOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_posTOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt04Sel"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_TOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("pt04Sel"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_negTOFreqRun3_posEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_negTOFreqRun3_negEta%s_strongNSigEPbPb_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));
        cut->AddCut(GetAnalysisCut(Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())));
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TPC_PID", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        cut->AddCut(cut_tpc_nSigma);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_TOF_PID", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        cut->AddCut(cut_tof_nSigma);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_strongNSigE_DCA05", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_DCA05"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID
      if (nameStr == Form("lmee%s_posNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_posNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_negNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID for low B field
      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_posEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_negEta%s_strongNSigE", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID for low B field with bad TOF rejection
      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_posEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_posEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("posEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_posNSigmaRun3_negEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("posTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_lowB_negNSigmaRun3_negEta%s_strongNSigE_rejBadTOF", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("negTrack"));
        cut->AddCut(GetAnalysisCut("negEtaSel"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      // some older cuts
      if (nameStr == Form("lmee%s_pp502TeV_PID%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
        if (nameStr == Form("lmee%s_pp502TeV_PID%s_UsePrefilter%d", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data(), i)) {
          cut->AddCut(GetAnalysisCut(Form("notDalitzLeg%d", i)));
          cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
          cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
          cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

          auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
          cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())));

          auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
          cut_tof_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())));

          auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
          cut_pid_OR->AddCut(cut_tpc_nSigma);
          cut_pid_OR->AddCut(cut_tof_nSigma);
          cut->AddCut(cut_pid_OR);
          return cut;
        }
      }

      if (nameStr == Form("lmee%s_pp502TeV_lowB_PID%s", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCAz")); // DCAz to reject loopers

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_lowB_TPC%s", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_lowB_TOF%s", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      if (nameStr == Form("lmee%s_eNSigmaRun3%s_pt04", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data())) {
        cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
        cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
        cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

        auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
        cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
        cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

        auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
        cut_pid_OR->AddCut(cut_tpc_nSigma);
        cut_pid_OR->AddCut(cut_tof_nSigma);
        cut->AddCut(cut_pid_OR);
        return cut;
      }

      for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
        if (nameStr == Form("lmee%s_eNSigmaRun3%s_UsePrefilter%d", vecTypetrackWithPID.at(jcase).Data(), vecPIDcase.at(icase).Data(), i)) {
          cut->AddCut(GetAnalysisCut(Form("notDalitzLeg%d", i)));
          cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
          cut->AddCut(GetAnalysisCut(Form("lmeeQCTrackCuts%s", vecTypetrackWithPID.at(jcase).Data())));
          cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

          auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
          cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

          auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
          cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

          auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
          cut_pid_OR->AddCut(cut_tpc_nSigma);
          cut_pid_OR->AddCut(cut_tof_nSigma);
          cut->AddCut(cut_pid_OR);
          return cut;
        }
      }
    }

    if (nameStr == Form("lmee_eNSigmaRun3%s_strongTPC", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3_strongTPC"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
      cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())));

      auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
      cut_tof_nSigma->AddCut(GetAnalysisCut(Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())));

      auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
      cut_pid_OR->AddCut(cut_tpc_nSigma);
      cut_pid_OR->AddCut(cut_tof_nSigma);
      cut->AddCut(cut_pid_OR);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_PID%s", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
      cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TPCloose%s", vecPIDcase.at(icase).Data())));

      auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
      cut_tof_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())));

      auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
      cut_pid_OR->AddCut(cut_tpc_nSigma);
      cut_pid_OR->AddCut(cut_tof_nSigma);
      cut->AddCut(cut_pid_OR);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtestb_PID%s", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));

      auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
      cut_tpc_nSigma->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TPCloosenopkrej%s", vecPIDcase.at(icase).Data())));

      cut->AddCut(cut_tpc_nSigma);
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF%s", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      cut->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF_pionrej%s", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
      cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
      cut->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())));
      return cut;
    }

    if (nameStr == Form("lmee_skimmingtesta_TOF_pionrej_noDCA%s", vecPIDcase.at(icase).Data())) {
      cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
      cut->AddCut(GetAnalysisCut("LooseGlobalTrackRun3"));
      cut->AddCut(GetAnalysisCut(Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())));
      return cut;
    }
  }

  if (nameStr == "testCut_chic") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine5"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun3") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_GlobalTrackRun3_lowPt") {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun3_lowPt") {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "lmee_TPCTrackRun3") {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCA"));
    return cut;
  }

  if (nameStr == "trackCut_compareDQEMframework") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("trackQuality_compareDQEMframework"));
    cut->AddCut(GetAnalysisCut("trackDCA1cm"));
    auto* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("lmee_commonDQEM_PID_TPC"));

    auto* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("lmee_commonDQEM_PID_TOF"));

    auto* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    return cut;
  }

  if (nameStr == "jpsi_debug_TPCTOF3_rejBadTOF") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine5"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly3"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("dcaCut1_ionut"));
    cut->AddCut(GetAnalysisCut("pidJpsi_TPCpion0"));
    cut->AddCut(GetAnalysisCut("pidJpsi_beta"));
    cut->AddCut(GetAnalysisCut("pidJpsi_noTOF_prot"));
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // lmee pair cuts

  if (nameStr == "pairPhiV") {
    auto* cut_pairPhiV = new AnalysisCompositeCut("cut_pairPhiV", "cut_pairPhiV", kTRUE);
    cut_pairPhiV->AddCut(GetAnalysisCut("pairLowMass"));
    cut_pairPhiV->AddCut(GetAnalysisCut("pairPhiV"));
    cut->AddCut(cut_pairPhiV);
    return cut;
  }

  if (nameStr == "excludePairPhiV") {
    auto* cut_pairlowPhiV = new AnalysisCompositeCut("cut_pairlowPhiV", "cut_pairlowPhiV", kFALSE);
    cut_pairlowPhiV->AddCut(GetAnalysisCut("excludePairLowMass"));
    cut_pairlowPhiV->AddCut(GetAnalysisCut("excludePairPhiV"));
    cut->AddCut(cut_pairlowPhiV);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts

  if (nameStr == "muonQualityCutsMatchingOnly") {
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonMinimalCuts") {
    cut->AddCut(GetAnalysisCut("muonMinimalCuts"));
    return cut;
  }

  if (nameStr == "muonMinimalCuts10SigmaPDCA") {
    cut->AddCut(GetAnalysisCut("muonMinimalCuts10SigmaPDCA"));
    return cut;
  }

  if (nameStr == "muonQualityCuts") {
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "matchedQualityCuts") {
    cut->AddCut(GetAnalysisCut("matchedQualityCuts"));
    return cut;
  }

  if (nameStr == "matchedQualityCutsMFTeta") {
    cut->AddCut(GetAnalysisCut("matchedQualityCutsMFTeta"));
    return cut;
  }

  if (nameStr == "muonQualityCuts5SigmaPDCA_Run3") {
    cut->AddCut(GetAnalysisCut("muonQualityCuts5SigmaPDCA_Run3"));
    return cut;
  }

  if (nameStr == "muonLowPt5SigmaPDCA_Run3") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts5SigmaPDCA_Run3"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonQualityCuts10SigmaPDCA_MCHMID") {
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt10SigmaPDCA") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt210SigmaPDCA") {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt510SigmaPDCA") {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt610SigmaPDCA") {
    cut->AddCut(GetAnalysisCut("muonLowPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts10SigmaPDCA"));
    cut->AddCut(GetAnalysisCut("MCHMID"));
    return cut;
  }

  if (nameStr == "muonLowPt") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt2") {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt3") {
    cut->AddCut(GetAnalysisCut("muonLowPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt4") {
    cut->AddCut(GetAnalysisCut("muonLowPt4"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt5") {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPt6") {
    cut->AddCut(GetAnalysisCut("muonLowPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly2") {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly3") {
    cut->AddCut(GetAnalysisCut("muonLowPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly4") {
    cut->AddCut(GetAnalysisCut("muonLowPt4"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly5") {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonLowPtMatchingOnly6") {
    cut->AddCut(GetAnalysisCut("muonLowPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPt") {
    cut->AddCut(GetAnalysisCut("muonHighPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt2") {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt3") {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt4") {
    cut->AddCut(GetAnalysisCut("muonHighPt4"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt5") {
    cut->AddCut(GetAnalysisCut("muonHighPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPt6") {
    cut->AddCut(GetAnalysisCut("muonHighPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly2") {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly3") {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonTightQualityCutsForTests") {
    cut->AddCut(GetAnalysisCut("muonTightQualityCutsForTests"));
    return cut;
  }

  if (nameStr == "mchTrack") {
    cut->AddCut(GetAnalysisCut("mchTrack"));
    return cut;
  }

  if (nameStr == "matchedMchMid") {
    cut->AddCut(GetAnalysisCut("matchedMchMid"));
    return cut;
  }

  if (nameStr == "matchedFwd") {
    cut->AddCut(GetAnalysisCut("matchedFwd"));
    return cut;
  }

  if (nameStr == "matchedGlobal") {
    cut->AddCut(GetAnalysisCut("matchedGlobal"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut1") {
    cut->AddCut(GetAnalysisCut("Chi2MCHMFTCut1"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut2") {
    cut->AddCut(GetAnalysisCut("Chi2MCHMFTCut2"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut3") {
    cut->AddCut(GetAnalysisCut("Chi2MCHMFTCut3"));
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut4") {
    cut->AddCut(GetAnalysisCut("Chi2MCHMFTCut4"));
    return cut;
  }

  if (nameStr == "muonQualityCutsMUONStandalone") {
    cut->AddCut(GetAnalysisCut("matchedMchMid"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (nameStr == "muonQualityCutsGlobal") {
    cut->AddCut(GetAnalysisCut("matchedGlobal"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }
  // -----------------------------------------------------------
  // Pair cuts
  if (nameStr == "pairNoCut") {
    cut->AddCut(GetAnalysisCut("pairNoCut"));
    return cut;
  }

  if (nameStr == "pairMassLow1") {
    cut->AddCut(GetAnalysisCut("pairMassLow1"));
    return cut;
  }

  if (nameStr == "pairMassLow2") {
    cut->AddCut(GetAnalysisCut("pairMassLow2"));
    return cut;
  }

  if (nameStr == "pairPtLow3") {
    cut->AddCut(GetAnalysisCut("pairPtLow3"));
    return cut;
  }

  if (nameStr == "pairPtLow4") {
    cut->AddCut(GetAnalysisCut("pairPtLow4"));
    return cut;
  }

  if (nameStr == "pairPtLow5") {
    cut->AddCut(GetAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairMassLow3") {
    cut->AddCut(GetAnalysisCut("pairMassLow3"));
    return cut;
  }

  if (nameStr == "pairMassLow4") {
    cut->AddCut(GetAnalysisCut("pairMassLow4"));
    return cut;
  }

  if (nameStr == "pairMassLow5") {
    cut->AddCut(GetAnalysisCut("pairMassLow5"));
    return cut;
  }

  if (nameStr == "pairMassLow6") {
    cut->AddCut(GetAnalysisCut("pairMassLow6"));
    return cut;
  }

  if (nameStr == "pairMassLow7") {
    cut->AddCut(GetAnalysisCut("pairMassLow7"));
    return cut;
  }

  if (nameStr == "pairMassLow8") {
    cut->AddCut(GetAnalysisCut("pairMassLow8"));
    return cut;
  }

  if (nameStr == "pairMassLow9") {
    cut->AddCut(GetAnalysisCut("pairMassLow9"));
    return cut;
  }

  if (nameStr == "pairMassLow10") {
    cut->AddCut(GetAnalysisCut("pairMassLow10"));
    return cut;
  }

  if (nameStr == "pairMassLow11") {
    cut->AddCut(GetAnalysisCut("pairMassLow11"));
    return cut;
  }

  if (nameStr == "pairMassLow12") {
    cut->AddCut(GetAnalysisCut("pairMassLow12"));
    return cut;
  }

  if (nameStr == "pairMass1to2") {
    cut->AddCut(GetAnalysisCut("pairMass1to2"));
    return cut;
  }

  if (nameStr == "pairMassIMR") {
    cut->AddCut(GetAnalysisCut("pairMassIMR"));
    return cut;
  }

  if (nameStr == "pairMass1_5to2_7") {
    cut->AddCut(GetAnalysisCut("pairMass1_5to2_7"));
    return cut;
  }

  if (nameStr == "pairMass1_3to3_5") {
    cut->AddCut(GetAnalysisCut("pairMass1_3to3_5"));
    return cut;
  }

  if (nameStr == "pairMass1_3") {
    cut->AddCut(GetAnalysisCut("pairMass1_3"));
    return cut;
  }

  if (nameStr == "pairMass1_5to3_5") {
    cut->AddCut(GetAnalysisCut("pairMass1_5to3_5"));
    return cut;
  }

  if (nameStr == "pairDalitz1") {
    cut->AddCut(GetAnalysisCut("pairDalitz1"));
    return cut;
  }

  if (nameStr == "pairDalitz1Strong") {
    cut->AddCut(GetAnalysisCut("pairDalitz1Strong"));
    return cut;
  }

  if (nameStr == "pairDalitz2") {
    cut->AddCut(GetAnalysisCut("pairDalitz2"));
    return cut;
  }

  if (nameStr == "pairDalitz3") {
    cut->AddCut(GetAnalysisCut("pairDalitz3"));
    return cut;
  }

  if (nameStr == "paira_prefilter1") {
    cut->AddCut(GetAnalysisCut("paira_prefilter1"));
    return cut;
  }

  if (nameStr == "paira_prefilter2") {
    cut->AddCut(GetAnalysisCut("paira_prefilter2"));
    return cut;
  }

  if (nameStr == "paira_prefilter3") {
    cut->AddCut(GetAnalysisCut("paira_prefilter3"));
    return cut;
  }

  if (nameStr == "paira_prefilter4") {
    cut->AddCut(GetAnalysisCut("paira_prefilter4"));
    return cut;
  }

  if (nameStr == "paira_prefilter5") {
    cut->AddCut(GetAnalysisCut("paira_prefilter5"));
    return cut;
  }

  if (nameStr == "paira_prefilter6") {
    cut->AddCut(GetAnalysisCut("paira_prefilter6"));
    return cut;
  }

  if (nameStr == "paira_prefilter7") {
    cut->AddCut(GetAnalysisCut("paira_prefilter7"));
    return cut;
  }

  if (nameStr == "pairb_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter1"));
    return cut;
  }

  if (nameStr == "pairb_prefilter2") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter2"));
    return cut;
  }

  if (nameStr == "pairb_prefilter3") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter3"));
    return cut;
  }

  if (nameStr == "pairb_prefilter4") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter4"));
    return cut;
  }

  if (nameStr == "pairb_prefilter5") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter5"));
    return cut;
  }

  if (nameStr == "pairb_prefilter6") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter6"));
    return cut;
  }

  if (nameStr == "pairb_prefilter7") {
    cut->AddCut(GetAnalysisCut("pairb_prefilter7"));
    return cut;
  }

  if (nameStr == "pairc_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter1"));
    return cut;
  }

  if (nameStr == "pairc_prefilter2") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter2"));
    return cut;
  }

  if (nameStr == "pairc_prefilter3") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter3"));
    return cut;
  }

  if (nameStr == "pairc_prefilter4") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter4"));
    return cut;
  }

  if (nameStr == "pairc_prefilter5") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter5"));
    return cut;
  }

  if (nameStr == "pairc_prefilter6") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter6"));
    return cut;
  }

  if (nameStr == "pairc_prefilter7") {
    cut->AddCut(GetAnalysisCut("pairc_prefilter7"));
    return cut;
  }

  if (nameStr == "paird_prefilter1") {
    cut->AddCut(GetAnalysisCut("paird_prefilter1"));
    return cut;
  }

  if (nameStr == "paire_prefilter1") {
    cut->AddCut(GetAnalysisCut("paire_prefilter1"));
    return cut;
  }

  if (nameStr == "pairf_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairf_prefilter1"));
    return cut;
  }

  if (nameStr == "pairg_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairg_prefilter1"));
    return cut;
  }

  if (nameStr == "pairh_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairh_prefilter1"));
    return cut;
  }

  if (nameStr == "pairi_prefilter1") {
    cut->AddCut(GetAnalysisCut("pairi_prefilter1"));
    return cut;
  }

  if (nameStr == "pairJpsi") {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    return cut;
  }

  if (nameStr == "pairJpsi2") {
    cut->AddCut(GetAnalysisCut("pairJpsi2"));
    return cut;
  }

  if (nameStr == "pairJpsi3") {
    cut->AddCut(GetAnalysisCut("pairJpsi3"));
    return cut;
  }

  if (nameStr == "pairPsi2S") {
    cut->AddCut(GetAnalysisCut("pairPsi2S"));
    return cut;
  }

  if (nameStr == "pairUpsilon") {
    cut->AddCut(GetAnalysisCut("pairUpsilon"));
    return cut;
  }

  if (nameStr == "pairX3872Cut1") {
    cut->AddCut(GetAnalysisCut("pairX3872"));
    return cut;
  }

  if (nameStr == "pairX3872Cut2") {
    cut->AddCut(GetAnalysisCut("pairX3872_2"));
    return cut;
  }

  if (nameStr == "pairX3872Cut3") {
    cut->AddCut(GetAnalysisCut("pairX3872_3"));
    return cut;
  }

  if (nameStr == "DipionPairCut1") {
    cut->AddCut(GetAnalysisCut("DipionMassCut1"));
    return cut;
  }

  if (nameStr == "DipionPairCut2") {
    cut->AddCut(GetAnalysisCut("DipionMassCut2"));
    return cut;
  }

  if (nameStr == "pairRapidityForward") {
    cut->AddCut(GetAnalysisCut("pairRapidityForward"));
    return cut;
  }

  if (nameStr == "pairJpsiLowPt1") {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    cut->AddCut(GetAnalysisCut("pairPtLow1"));
    return cut;
  }

  if (nameStr == "pairJpsiLowPt2") {
    cut->AddCut(GetAnalysisCut("pairJpsi"));
    cut->AddCut(GetAnalysisCut("pairPtLow2"));
    return cut;
  }

  if (nameStr == "pairCoherentRho0") {
    cut->AddCut(GetAnalysisCut("pairPtLow3"));
    return cut;
  }

  if (nameStr == "pairD0") {
    cut->AddCut(GetAnalysisCut("pairD0"));
    return cut;
  }

  if (nameStr == "pairD0HighPt1") {
    cut->AddCut(GetAnalysisCut("pairLxyzProjected3sigma"));
    cut->AddCut(GetAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairD0HighPt2") {
    cut->AddCut(GetAnalysisCut("pairTauxyzProjected1"));
    cut->AddCut(GetAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairD0HighPt3") {
    cut->AddCut(GetAnalysisCut("pairTauxyzProjected1sigma"));
    cut->AddCut(GetAnalysisCut("pairPtLow5"));
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1") {
    cut->AddCut(GetAnalysisCut("pairTauxyzProjected1"));
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaLambdacCand") {
    cut->AddCut(GetAnalysisCut("pairLxyProjected3sigmaLambdacCand"));
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaDplusCand") {
    cut->AddCut(GetAnalysisCut("pairLxyProjected3sigmaDplusCand"));
    return cut;
  }

  if (nameStr == "pairCosPointingPos") {
    cut->AddCut(GetAnalysisCut("pairCosPointingPos"));
    return cut;
  }

  if (nameStr == "pairCosPointingNeg90") {
    cut->AddCut(GetAnalysisCut("pairCosPointingNeg90"));
    return cut;
  }

  if (nameStr == "pairCosPointingNeg85") {
    cut->AddCut(GetAnalysisCut("pairCosPointingNeg85"));
    return cut;
  }

  if (nameStr == "pairTauxyzProjectedCosPointing1") {
    cut->AddCut(GetAnalysisCut("pairCosPointingNeg"));
    cut->AddCut(GetAnalysisCut("pairTauxyzProjected1"));
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and in order or optimize the trigger
  // trigger selection cuts

  if (nameStr == "jpsiO2TriggerTestCuts_LooseNsigma") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_LooseNsigma_corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));
    return cut;
  }
  if (nameStr == "jpsiO2TriggerTestCuts_MediumNsigma") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_MediumNsigma_corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (nameStr == "jpsiO2TriggerTestCuts_TightNsigma") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TightNsigma_corr") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID1") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest1"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID2") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest2"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID3") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest3"));
    return cut;
  }

  if (nameStr == "jpsiO2TriggerTestCuts_TPCPID4") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest4"));
    return cut;
  }

  if (nameStr == "emu_electron_test") {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronTrackQuality_Maolin"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaEMu"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts") {
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_LowPt") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_HighPt2") {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts_HighPt3") {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly2") {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonHighPtMatchingOnly3") {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts") {
    cut->AddCut(GetAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts_LowPt") {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  //---------------------------------------------------------------
  // ALICE 3 studies composite cuts

  if (nameStr == "alice3StandardKine") {
    cut->AddCut(GetAnalysisCut("alice3StandardKine"));
    return cut;
  }

  if (nameStr == "alice3KineSkim") {
    cut->AddCut(GetAnalysisCut("alice3KineSkim"));
    return cut;
  }

  if (nameStr == "alice3TrackQuality") {
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TrackQualityTightDCA") {
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrack") {
    cut->AddCut(GetAnalysisCut("alice3StandardKine"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackTOFAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3KineTOFAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3StandardTrackRICHAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3KineRICHAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrack") {
    cut->AddCut(GetAnalysisCut("alice3StandardKine"));
    cut->AddCut(GetAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrackTOFAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3KineTOFAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3TightDCATrackRICHAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3KineRICHAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQualityTightDCA"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDEl") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPi") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDKa") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPr") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDEl") {
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPi") {
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDKa") {
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPr") {
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDEl") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDPi") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPi"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPi"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDKa") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDKa"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDKa"));
    return cut;
  }

  if (nameStr == "alice3FullTOFPIDPr") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPr"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPr"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDEl") {
    cut->AddCut(GetAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDPi") {
    cut->AddCut(GetAnalysisCut("alice3RICHPIDPi"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDKa") {
    cut->AddCut(GetAnalysisCut("alice3RICHPIDKa"));
    return cut;
  }

  if (nameStr == "alice3RICHPIDPr") {
    cut->AddCut(GetAnalysisCut("alice3RICHPIDPr"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDEl") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDPi") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPi"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPi"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDPi"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDKa") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDKa"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDKa"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDKa"));
    return cut;
  }

  if (nameStr == "alice3RICHTOFPIDPr") {
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDPr"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDPr"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDPr"));
    return cut;
  }

  if (nameStr == "alice3DielectronPID") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKine"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFOnly") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDRICHOnly") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFOnly") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKine"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDTOFAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKineTOFAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  if (nameStr == "alice3DielectronPIDRICHAcceptance") {
    cut->AddCut(GetAnalysisCut("alice3JpsiKineRICHAcceptance"));
    cut->AddCut(GetAnalysisCut("alice3TrackQuality"));
    cut->AddCut(GetAnalysisCut("alice3iTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3oTOFPIDEl"));
    cut->AddCut(GetAnalysisCut("alice3RICHPIDEl"));
    return cut;
  }

  delete cut;
  LOGF(fatal, Form("Did not find cut %s. Returning nullptr", cutName));
  return nullptr;
}

AnalysisCut* o2::aod::dqcuts::GetAnalysisCut(const char* cutName)
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
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventIsTVXTriggered") {
    cut->AddCut(VarManager::kIsTVXTriggered, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandard") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    return cut;
  }

  if (nameStr == "eventStandardtest") {
    cut->AddCut(VarManager::kVtxZ, -30.0, 30.0);
    return cut;
  }

  if (nameStr == "eventSel8") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    return cut;
  }
  if (nameStr == "eventStandardSel8") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    return cut;
  }
  if (nameStr == "eventStandardSel8WithITSROFRecomputedCut") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorderRecomputed, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBorder") { // Redundant w.r.t. eventStandardSel8, to be removed
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBNoITSROFB") { // Redundant w.r.t. eventStandardSel8, to be removed
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoTFBNoITSROFBrecomp") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorderRecomputed, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8NoSameBunch") {
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8TriggerZNAZNC") {
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTriggerZNAZNC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8TriggerZNAZNCNoPileUp") {
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTriggerZNAZNC, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSel8NoSameBunchGoodZvtx") {
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQuality") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8NoPileup") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityCent90") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityGoodITSLayersAll") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodITSLayersAll, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancy") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0., 1000);
    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityFirmTrackOccupancy") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0., 2000);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityLooseTrackOccupancy") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0., 5000);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancyCollInTime") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0., 500);
    cut->AddCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);

    return cut;
  }

  if (nameStr == "eventStandardSel8PbPbQualityTightTrackOccupancyCollInTime") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0., 1000);
    cut->AddCut(VarManager::kNoCollInTimeRangeStandard, 0.5, 1.5);

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
      cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
      cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
      cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
      cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
      cut->AddCut(VarManager::kTrackOccupancyInTimeRange, vecOccupancies[icase], vecOccupancies[icase + 1]);

      return cut;
    }
  }

  for (size_t icase = 0; icase < vecOccupancies.size() - 1; icase++) {
    if (nameStr == Form("eventStandardSel8PbPbQualityTrackOccupancySlice_0_%lu", icase)) {
      cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
      cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
      cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
      cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
      cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
      cut->AddCut(VarManager::kTrackOccupancyInTimeRange, 0, vecOccupancies[icase]);

      return cut;
    }
  }

  if (nameStr == "eventStandardSel8ppQuality") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8ppQualityNoVtxZ") {
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8multAnalysis") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventStandardSel8VtxQuality1") { // kIsSel8 = kIsTriggerTVX && kNoITSROFrameBorder && kNoTimeFrameBorder
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexITSTPC, 0.5, 1.5);
    cut->AddCut(VarManager::kIsVertexTOFmatched, 0.5, 1.5);
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
    cut->AddCut(VarManager::kVtxNcontribReal, fMultPVCutLow, fMultPVCutHigh, false, VarManager::kCentFT0C, 0.0, 100.0, false);
    cut->AddCut(VarManager::kMultA, fMultCutLow, fMultCutHigh, false, VarManager::kCentFT0C, 0.0, 100.0, false);
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventDimuonStandard") {
    cut->AddCut(VarManager::kIsMuonUnlikeLowPt7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventMuonStandard") {
    cut->AddCut(VarManager::kIsMuonSingleLowPt7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventTPCMultLow") {
    cut->AddCut(VarManager::kMultTPC, 0, 50);
    return cut;
  }

  if (nameStr == "eventExclusivePair") {
    cut->AddCut(VarManager::kVtxNcontrib, 2, 2);
    return cut;
  }

  if (nameStr == "eventVtxNContrib") {
    cut->AddCut(VarManager::kVtxNcontrib, 0, 10);
    return cut;
  }

  if (nameStr == "eventTPCMult3") {
    cut->AddCut(VarManager::kMultTPC, 3, 3);
    return cut;
  }

  if (nameStr == "int7vtxZ5") {
    cut->AddCut(VarManager::kVtxZ, -5.0, 5.0);
    cut->AddCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventDoubleGap") {
    cut->AddCut(VarManager::kIsDoubleGap, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGap") {
    cut->AddCut(VarManager::kIsSingleGap, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapA") {
    cut->AddCut(VarManager::kIsSingleGapA, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapAZDC") {
    cut->AddCut(VarManager::kIsSingleGapA, 0.5, 1.5);
    cut->AddCut(VarManager::kEnergyCommonZNA, -1000., 1.);
    cut->AddCut(VarManager::kEnergyCommonZNC, 1., 1000.);
    return cut;
  }

  if (nameStr == "eventSingleGapC") {
    cut->AddCut(VarManager::kIsSingleGapC, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapCZDC") {
    cut->AddCut(VarManager::kIsSingleGapC, 0.5, 1.5);
    cut->AddCut(VarManager::kEnergyCommonZNC, -1000., 1.);
    cut->AddCut(VarManager::kEnergyCommonZNA, 1., 1000.);
    return cut;
  }

  if (nameStr == "eventSingleGapACZDC") {
    auto* cutA = new AnalysisCompositeCut("singleGapAZDC", "singleGapAZDC", kTRUE);
    cutA->AddCut(GetAnalysisCut("eventSingleGapAZDC"));

    auto* cutC = new AnalysisCompositeCut("singleGapCZDC", "singleGapCZDC", kTRUE);
    cutC->AddCut(GetAnalysisCut("eventSingleGapCZDC"));

    auto* cutAorC = new AnalysisCompositeCut("singleGapACZDC", "singleGapACZDC", kFALSE);
    cutAorC->AddCut(cutA);
    cutAorC->AddCut(cutC);
    return cutAorC;
  }

  if (nameStr == "eventUPCMode") {
    cut->AddCut(VarManager::kIsITSUPCMode, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "eventSingleGapACZDC_UPCMode") {
    auto* cutA = new AnalysisCompositeCut("singleGapAZDC", "singleGapAZDC", kTRUE);
    cutA->AddCut(GetAnalysisCut("eventSingleGapAZDC"));
    cutA->AddCut(GetAnalysisCut("eventUPCMode"));

    auto* cutC = new AnalysisCompositeCut("singleGapCZDC", "singleGapCZDC", kTRUE);
    cutC->AddCut(GetAnalysisCut("eventSingleGapCZDC"));
    cutC->AddCut(GetAnalysisCut("eventUPCMode"));

    auto* cutAorC = new AnalysisCompositeCut("singleGapACZDC", "singleGapACZDC", kFALSE);
    cutAorC->AddCut(cutA);
    cutAorC->AddCut(cutC);
    return cutAorC;
  }

  if (nameStr == "eventXn0nTime") {
    cut->AddCut(VarManager::kTimeZNA, -2.0, 2.0);
    cut->AddCut(VarManager::kTimeZNC, -2.0, 2.0, true);
    return cut;
  }

  if (nameStr == "event0nXnTime") {
    cut->AddCut(VarManager::kTimeZNA, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTimeZNC, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "event0n0nTime") {
    cut->AddCut(VarManager::kTimeZNA, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTimeZNC, -2.0, 2.0, true);
    return cut;
  }

  if (nameStr == "eventXnXnTime") {
    cut->AddCut(VarManager::kTimeZNA, -2.0, 2.0);
    cut->AddCut(VarManager::kTimeZNC, -2.0, 2.0);
    return cut;
  }

  // Event cuts based on centrality
  if (nameStr == "eventStandardNoINT7Cent090") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent7090") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 70.0, 90.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent5070") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 50.0, 70.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent3050") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 30.0, 50.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent1030") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 10.0, 30.0);
    return cut;
  }

  if (nameStr == "eventStandardNoINT7Cent010") {
    cut->AddCut(VarManager::kVtxZ, -10.0, 10.0);
    cut->AddCut(VarManager::kCentFT0C, 0.0, 10.0);
    return cut;
  }

  // ---------------------------------------------------
  // Barrel track kine cuts
  if (nameStr == "negTrack") {
    cut->AddCut(VarManager::kCharge, -99., 0.);
    return cut;
  }

  if (nameStr == "posTrack") {
    cut->AddCut(VarManager::kCharge, 0., 99.);
    return cut;
  }

  if (nameStr == "posEtaSel") {
    cut->AddCut(VarManager::kEta, 0., 0.8);
    return cut;
  }

  if (nameStr == "negEtaSel") {
    cut->AddCut(VarManager::kEta, -0.8, 0.);
    return cut;
  }

  if (nameStr == "etaSel") {
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "pt02Sel") {
    cut->AddCut(VarManager::kPt, 0.2, 20.0);
    return cut;
  }

  if (nameStr == "pt04Sel") {
    cut->AddCut(VarManager::kPt, 0.4, 20.0);
    return cut;
  }

  if (nameStr == "openEtaSel") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "insideTPCsector") {
    cut->AddCut(VarManager::kTrackIsInsideTPCModule, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "rho0Kine") {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -1.1, 1.1);
    return cut;
  }

  if (nameStr == "pionQuality") {
    cut->AddCut(VarManager::kTPCncls, 50.0, 1000.);
    return cut;
  }

  if (nameStr == "primaryVertexContributor") {
    cut->AddCut(VarManager::kPVContributor, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "jpsiStandardKine") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine2") {
    cut->AddCut(VarManager::kPt, 0.9, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine3") {
    cut->AddCut(VarManager::kP, 1.2, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiPIDcalibKine_posEta") {
    cut->AddCut(VarManager::kPin, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, 0.0, 0.9);
    return cut;
  }

  if (nameStr == "jpsiPIDcalibKine_negEta") {
    cut->AddCut(VarManager::kPin, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.0);
    return cut;
  }

  if (nameStr == "jpsiStandardKine4") {
    cut->AddCut(VarManager::kP, 1.5, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiStandardKine5") {
    cut->AddCut(VarManager::kP, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiKineSkimmed") {
    cut->AddCut(VarManager::kPt, 0.7, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "jpsiKineSkimmed") {
    cut->AddCut(VarManager::kPt, 0.7, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "lmeePrefilterKine") {
    cut->AddCut(VarManager::kPt, 0., 20.0);
    cut->AddCut(VarManager::kEta, -1.2, 1.2);
    return cut;
  }

  if (nameStr == "lmeeStandardKine") {
    cut->AddCut(VarManager::kPt, 0.2, 20.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "lmeeStandardKine_pt04") {
    cut->AddCut(VarManager::kPt, 0.4, 20.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "lmeeLowBKine") {
    cut->AddCut(VarManager::kPt, 0.075, 20.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "dalitzStandardKine") {
    cut->AddCut(VarManager::kPt, 0.15, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "PIDStandardKine") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "PIDStandardKine2") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    return cut;
  }

  if (nameStr == "PIDStandardKine3") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kPt, 0.5, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow04") {
    cut->AddCut(VarManager::kPt, 0.4, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow03") {
    cut->AddCut(VarManager::kPt, 0.3, 1000.0);
    return cut;
  }

  if (nameStr == "pTLow02") {
    cut->AddCut(VarManager::kPt, 0.2, 1000.0);
    return cut;
  }

  // -----------------------------------------------
  // Barrel track quality cuts

  // ---------------------------------------------------
  // MC generated particle acceptance cuts

  if (nameStr == "rapidity08") {
    cut->AddCut(VarManager::kMCY, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "rapidity09") {
    cut->AddCut(VarManager::kMCY, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_pp13600") {
    cut->AddCut(VarManager::kMCY, -0.8, 0.8);
    cut->AddCut(VarManager::kMCPt1, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCPt2, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCEta1, -0.8, 0.8);
    cut->AddCut(VarManager::kMCEta2, -0.8, 0.8);
    return cut;
  }

  if (nameStr == "acceptance_pp5360") {
    cut->AddCut(VarManager::kMCY, -0.9, 0.9);
    cut->AddCut(VarManager::kMCPt1, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCPt2, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCEta1, -0.9, 0.9);
    cut->AddCut(VarManager::kMCEta2, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_PbPb5360") {
    cut->AddCut(VarManager::kMCY, -0.9, 0.9);
    cut->AddCut(VarManager::kMCP1, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCP2, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCEta1, -0.9, 0.9);
    cut->AddCut(VarManager::kMCEta2, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "acceptance_PbPb5360_y08") {
    cut->AddCut(VarManager::kMCY, -0.8, 0.8);
    cut->AddCut(VarManager::kMCP1, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCP2, 1.0, 1000.0);
    cut->AddCut(VarManager::kMCEta1, -0.8, 0.8);
    cut->AddCut(VarManager::kMCEta2, -0.8, 0.8);
    return cut;
  }

  // ---------------------------------------------------
  // MC generated particle acceptance cuts

  // Run 2 only

  if (nameStr == "highPtHadron") {
    cut->AddCut(VarManager::kPt, 4.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  if (nameStr == "TightGlobalTrack") {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kTPCncls, 60.0, 170.);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQuality") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityBenchmark") {
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 70.0, 161.);
    return cut;
  }

  // Run 2 or run 3

  if (nameStr == "jpsi_trackCut_debug") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    cut->AddCut(VarManager::kITSncls, 2.5, 7.5);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_noITSCuts_debug") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug2") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    cut->AddCut(VarManager::kITSncls, 2.5, 7.5);
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTrackDCAxy, -1, 1);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug3") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCnclsCR, 70., 159);
    cut->AddCut(VarManager::kITSncls, 2.5, 7.5);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTrackDCAxy, -0.05, 0.05);
    cut->AddCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug4") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    cut->AddCut(VarManager::kITSncls, 2.5, 7.5);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTrackDCAxy, -1, 1);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug5") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70., 159);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "jpsi_trackCut_debug6") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 120., 159);
    cut->AddCut(VarManager::kTPCnclsCR, 140., 159);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "lmee_trackCut_debug") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 80., 159);
    cut->AddCut(VarManager::kITSncls, 0.0, 7.5);
    return cut;
  }

  if (nameStr == "lmee_skimming_cuts") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 60.0, 161.);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
    return cut;
  }

  if (nameStr == "trackQuality_compareDQEMframework") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kTPCnCRoverFindCls, 0.8, 1e+10);
    cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
    return cut;
  }

  if ((nameStr == "TightGlobalTrackRun3") || (nameStr == "lmeeQCTrackCuts")) {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
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
        cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 6.5, 7.5);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
      } else if (icase == 3) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
        cut->AddCut(VarManager::kHasTPC, -0.5, 0.5);
      } else if (icase == 4) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
      } else if (icase == 5) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
        cut->AddCut(VarManager::kHasITS, -0.5, 0.5);
      } else if (icase == 6) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
        cut->AddCut(VarManager::kHasTRD, -0.5, 0.5);
      } else if (icase == 7) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 130.0, 170.);
      } else if (icase == 8) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
      } else {
        cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
        cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
        cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
        cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
        cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
        cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
      }
      return cut;
    }
  }

  if (nameStr == "LooseGlobalTrackRun3") {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 6.0);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    cut->AddCut(VarManager::kTPCncls, 70.0, 170.);
    return cut;
  }

  if (nameStr == "TightGlobalTrackRun3_strongTPC") {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
    cut->AddCut(VarManager::kTPCnclsCR, 140.0, 161.);
    cut->AddCut(VarManager::kTPCncls, 120.0, 170.);
    return cut;
  }

  if (nameStr == "TightTPCTrackRun3") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    return cut;
  }

  if (nameStr == "TightTPCTrack") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "SPDfirst") {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "noTPC") {
    cut->AddCut(VarManager::kHasTPC, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "SPDany") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "ITSiball") {
    cut->AddCut(VarManager::kIsITSibAll, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "ITSibany") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    return cut;
  }

  // List of 30 variations in ITS and TPC parameters
  std::vector<double> cutVar_ITSchi2 = {6., 6., 5., 4., 4., 6., 6., 5., 4., 5., 4., 5., 6., 5., 6., 5., 6., 5., 5., 4., 6., 4., 6., 5., 6., 4., 4., 6., 4., 5.};
  std::vector<double> cutVar_TPCchi2 = {5., 5., 4., 3., 5., 4., 5., 3., 5., 4., 5., 3., 3., 5., 4., 5., 3., 5., 5., 5., 3., 5., 5., 4., 3., 4., 5., 5., 5., 3.};
  std::vector<double> cutVar_ITSnCl = {4.5, 5.5, 5.5, 4.5, 6.5, 4.5, 4.5, 4.5, 4.5, 5.5, 4.5, 5.5, 5.5, 5.5, 5.5, 4.5, 5.5, 6.5, 5.5, 4.5, 4.5, 5.5, 5.5, 5.5, 6.5, 5.5, 4.5, 4.5, 6.5, 6.5};
  std::vector<double> cutVar_TPCnClsCR = {90., 80., 80., 80., 90., 80., 70., 90., 70., 80., 70., 90., 90., 70., 90., 90., 70., 80., 90., 80., 80., 90., 70., 70., 70., 80., 90., 70., 70., 80.};
  std::vector<double> cutVar_TPCnCls = {80., 100., 80., 90., 90., 80., 80., 80., 80., 90., 100., 100., 80., 80., 80., 80., 100., 90., 100., 90., 90., 100., 100., 80., 100., 90., 90., 100., 90., 90.};

  for (unsigned int i = 0; i < cutVar_ITSchi2.size(); i++) {
    if (nameStr == Form("lmeeCutVarTrackCuts%i", i)) {
      cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
      cut->AddCut(VarManager::kITSchi2, 0.0, cutVar_ITSchi2.at(i));
      cut->AddCut(VarManager::kTPCchi2, 0.0, cutVar_TPCchi2.at(i));
      cut->AddCut(VarManager::kITSncls, cutVar_ITSnCl.at(i), 7.5);
      cut->AddCut(VarManager::kTPCnclsCR, cutVar_TPCnClsCR.at(i), 161.);
      cut->AddCut(VarManager::kTPCncls, cutVar_TPCnCls.at(i), 170.);
      return cut;
    }
  }

  if (nameStr == "electronStandardQualityForO2MCdebug") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug2") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 100.0, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug3") {
    cut->AddCut(VarManager::kITSncls, 0.5, 10);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityForO2MCdebug4") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityITSOnly") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQualitybAnyITSOnly") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly2") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 100, 161.);
    return cut;
  }

  if (nameStr == "electronStandardQualityTPCOnly3") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 120, 161.);
    return cut;
  }

  if (nameStr == "NoelectronStandardQualityTPCOnly") {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0, true, VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCncls, 60, 161);
    cut->AddCut(VarManager::kTrackDCAxy, -1.5, 1.5);
    cut->AddCut(VarManager::kTrackDCAz, -1.5, 1.5);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed2") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCncls, 60, 161);
    return cut;
  }

  if (nameStr == "electronTrackQualitySkimmed3") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCnclsCR, 70, 161);
    cut->AddCut(VarManager::kTPCncls, 70, 161);
    return cut;
  }

  if (nameStr == "electronTrackQuality_Maolin") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 15.0);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    cut->AddCut(VarManager::kTPCnclsCR, 70, 161);
    cut->AddCut(VarManager::kTrackDCAxy, -2.0, 2.0);
    cut->AddCut(VarManager::kTrackDCAz, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "pionQualityCut1") {
    cut->AddCut(VarManager::kPt, 0.15, 1000.0);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCncls, 70, 161);
    return cut;
  }

  if (nameStr == "pionQualityCut2") {
    cut->AddCut(VarManager::kPt, 0.15, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCncls, 90, 161);
    cut->AddCut(VarManager::kTPCnclsCR, 70, 161);
    return cut;
  }

  if (nameStr == "protonPVcut") {
    cut->AddCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->AddCut(VarManager::kTrackDCAz, -0.15, 0.15);
    cut->AddCut(VarManager::kPt, 0.4, 3);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 80, 161.);
    return cut;
  }

  if (nameStr == "pidbasic") {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCncls, 60, 161.);
    return cut;
  }

  if (nameStr == "standardPrimaryTrack") {
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "trackDCA1cm") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->AddCut(VarManager::kTrackDCAz, -1.0, 1.0);
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
      cut->AddCut(VarManager::kTrackDCAxy, fDCAxyresLow, fDCAxyresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->AddCut(VarManager::kTrackDCAz, fDCAzresLow, fDCAzresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->AddCut(VarManager::kIsITSibFirst, 0.5, 1.5);
      cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
      cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
      cut->AddCut(VarManager::kTPCncls, 100, 161.);
      return cut;
    }

    if (nameStr == "electronPrimary_dca7sigma") {
      // DCAxy and DCAz 7 sigma cut
      fDCAxyresLow->SetParameters(-7 * 8.7e-4 - 20e-4, -7 * 25.4e-4, 0.79);
      fDCAzresLow->SetParameters(-7 * 9.4e-4 - 10e-4, -7 * 26.5e-4, 0.79);
      fDCAxyresUp->SetParameters(7 * 8.7e-4 + 20e-4, 7 * 25.4e-4, 0.79);
      fDCAzresUp->SetParameters(7 * 9.4e-4 + 10e-4, 7 * 26.5e-4, 0.79);
      cut->AddCut(VarManager::kTrackDCAxy, fDCAxyresLow, fDCAxyresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->AddCut(VarManager::kTrackDCAz, fDCAzresLow, fDCAzresUp, false, VarManager::kPt, 0.2, 1000.);
      cut->AddCut(VarManager::kIsITSibFirst, 0.5, 1.5);
      cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
      cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
      cut->AddCut(VarManager::kTPCncls, 100, 161.);
      return cut;
    }
  }

  if (nameStr == "dcaCut1_ionut") {
    cut->AddCut(VarManager::kTrackDCAxy, -0.5, 0.5);
    cut->AddCut(VarManager::kTrackDCAz, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "trackQuality_ionut") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCncls, 70, 161);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 2.0);
    return cut;
  }

  if (nameStr == "trackQualityTight_ionut") {
    cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCncls, 100, 161);
    cut->AddCut(VarManager::kITSchi2, 0.0, 3.0);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 2.0);
    cut->AddCut(VarManager::kITSncls, 5.0, 8.0);
    return cut;
  }

  if (nameStr == "kineJpsiEle_ionut") {
    cut->AddCut(VarManager::kP, 1.0, 15.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "pidJpsiEle0_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 4.0, false, VarManager::kPin, 1.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kPin, 4.0, 150.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, 98.1, 98.11, false, VarManager::kPin, 0.0, 1.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle1_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0, false, VarManager::kPin, 1.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kPin, 4.0, 150.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, 98.1, 98.11, false, VarManager::kPin, 0.0, 1.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle2_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle3_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -0.5, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle4_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle5_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, 0.5, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle6_ionut") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    cut->AddCut(VarManager::kTOFnSigmaEl, -1.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -4.0, 4.0, true);
    return cut;
  }

  if (nameStr == "pidJpsiEle7_ionut") {
    cut->AddCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsiEle8_ionut") {
    cut->AddCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.5, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsiEle9_ionut") {
    cut->AddCut(VarManager::kTOFnSigmaEl, -3.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 4.0);
    return cut;
  }

  if (nameStr == "pidJpsi_TPCpion0") {
    cut->AddCut(VarManager::kTPCnSigmaPi, 4.0, 1000.0);
    return cut;
  }

  if (nameStr == "pidJpsi_noTOF_prot") {
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 1000.0, false, VarManager::kHasTOF, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "pidJpsi_beta") {
    cut->AddCut(VarManager::kTOFbeta, 0.98, 1.02, false, VarManager::kHasTOF, 0.5, 1.5);
    return cut;
  }

  // Magnus cuts ----------------------------------------------------------

  if (nameStr == "pidJpsi_magnus_ele1") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_ele2") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_ele3") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_prot1") {
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_prot2") {
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_pion1") {
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 1000.0);
    return cut;
  }
  if (nameStr == "pidJpsi_magnus_pion2") {
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.5, 1000.0);
    return cut;
  }

  // ----------------------------------------------------------------------------------

  if (nameStr == "standardPrimaryTrackDCAz") {
    cut->AddCut(VarManager::kTrackDCAxy, -3.0, 3.0);
    cut->AddCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "standardPrimaryTrackDCA") {
    cut->AddCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->AddCut(VarManager::kTrackDCAz, -0.15, 0.15);
    return cut;
  }

  if (nameStr == "PrimaryTrack_looseDCA") {
    cut->AddCut(VarManager::kTrackDCAxy, -3.0, 3.0);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "tightPrimaryTrack") {
    cut->AddCut(VarManager::kTrackDCAsigXY, -3.0, 3.0);
    cut->AddCut(VarManager::kTrackDCAsigZ, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCA05") {
    cut->AddCut(VarManager::kTrackDCAsigXY, -0.5, 0.5);
    cut->AddCut(VarManager::kTrackDCAsigZ, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "PrimaryTrack_DCAz") {
    cut->AddCut(VarManager::kTrackDCAz, -0.3, 0.3);
    return cut;
  }

  if (nameStr == "hasTOF") {
    cut->AddCut(VarManager::kHasTOF, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "noTOF") {
    cut->AddCut(VarManager::kHasTOF, -0.5, 0.5);
    return cut;
  }

  // -----------------------------------------------------
  // V0 and Dalitz legs selections

  for (int i = 1; i <= 8; i++) { // o2-linter: disable=magic-number (number of cuts)
    if (nameStr == Form("dalitzLeg%d", i)) {
      cut->AddCut(VarManager::kIsDalitzLeg + i - 1, 0.5, 1.5);
      return cut;
    }

    if (nameStr == Form("notDalitzLeg%d", i)) {
      cut->AddCut(VarManager::kIsDalitzLeg + i - 1, -0.5, 0.5);
      return cut;
    }
  }

  if (nameStr == "pidcalib_ele") {
    cut->AddCut(VarManager::kIsLegFromGamma, 0.5, 1.5, false);
    return cut;
  }

  if (nameStr == "pidcalib_pion") {
    cut->AddCut(VarManager::kIsLegFromK0S, 0.5, 1.5, false);
    return cut;
  }

  if (nameStr == "pidcalib_proton") {
    cut->AddCut(VarManager::kIsProtonFromLambdaAndAntiLambda, 0.5, 1.5, false);
    return cut;
  }

  if (nameStr == "pidcalib_kaon") {
    cut->AddCut(VarManager::kTOFnSigmaKa, -2.0, 2.0);
    cut->AddCut(VarManager::kTOFnSigmaPi, -2.0, 2.0, true);
    cut->AddCut(VarManager::kITSncls, 1.5, 7.5);
    return cut;
  }

  // ------------------------------------------------
  // Barrel PID cuts
  if (nameStr == "jpsi_TPCPID_debug1") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug2") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug3") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug4") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug5") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug5_noCorr") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999);
    return cut;
  }
  if (nameStr == "electronPIDLooseSkimmed") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999);
    return cut;
  }

  if (nameStr == "electronPIDLooseSkimmed2") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 999, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 1.5, 999, false, VarManager::kPin, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 1.5, 999, false, VarManager::kPin, 3.0, 999);
    return cut;
  }

  if (nameStr == "electronPIDLooseSkimmed3") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 0, 3.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug6") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug7") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug8") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 3.0, false, VarManager::kPin, 0.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0, false, VarManager::kPin, 3.0, 999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kPin, 0.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kPin, 5.0, 999.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999, false, VarManager::kPin, 3.0, 999.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug9") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.5, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 1.0, 999, false, VarManager::kPin, 3.0, 999.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.0, 999);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_debug10") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.5, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 4.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 999);
    return cut;
  }

  if (nameStr == "pidCut_lowP_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0, false, VarManager::kP, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kP, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kP, 0.0, 5.0);
    return cut;
  }

  if (nameStr == "EleInclusion_highP_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1.0, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "EleInclusion_highP2_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -0.5, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "PionExclusion_highP_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.0, 999, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "pidCut_lowP") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0, false, VarManager::kP, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999, false, VarManager::kP, 0.0, 5.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 999, false, VarManager::kP, 0.0, 5.0);
    return cut;
  }

  if (nameStr == "EleInclusion_highP") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "EleInclusion_highP2") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -0.5, 4.0, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }
  if (nameStr == "PionExclusion_highP") {
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 999, false, VarManager::kP, 5.0, 999.0);
    return cut;
  }

  if (nameStr == "lmee_TPCPID_debug1") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
    return cut;
  }

  if (nameStr == "electronPID1" || nameStr == "electronPID1shiftUp" || nameStr == "electronPID1shiftDown" || nameStr == "electronPID2" || nameStr == "electronPID3") {
    std::shared_ptr<TF1> cutLow1 = std::make_shared<TF1>("cutLow1", "pol1", 0., 10.);
    if (nameStr == "electronPID1") {
      cutLow1->SetParameters(130., -40.0);
      cut->AddCut(VarManager::kTPCsignal, 70., 100.);
      cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID1shiftUp") {
      cut->AddCut(VarManager::kTPCsignal, 70. - 0.85, 100. - 0.85);
      cutLow1->SetParameters(130. - 0.85, -40.0);
      cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0 - 0.85, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID1shiftDown") {
      cut->AddCut(VarManager::kTPCsignal, 70.0 + 0.85, 100.0 + 0.85);
      cutLow1->SetParameters(130. + 0.85, -40.0);
      cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0 + 0.85, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID2") {
      cutLow1->SetParameters(130., -40.0);
      cut->AddCut(VarManager::kTPCsignal, 73., 100.);
      cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }

    if (nameStr == "electronPID3") {
      cutLow1->SetParameters(130., -40.0);
      cut->AddCut(VarManager::kTPCsignal, 60., 110.);
      cut->AddCut(VarManager::kTPCsignal, cutLow1, 100.0, false, VarManager::kPin, 0.5, 3.0);
      return cut;
    }
  }

  if (nameStr == "electronPIDnsigma") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 3000.0);
    return cut;
  }

  if (nameStr == "lmee_commonDQEM_PID_TPC") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.5, 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -1e12, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "lmee_commonDQEM_PID_TOF") { // cut setting to check least common factor between reduced data sets of PWGEM and PWGDQ
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.5, 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
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
        cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_loose", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TPCnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TPCnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TPCnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    // List of nSigma values for lower and upper edge of El, Pi, Ka, Pr, TPC PID selections/rejection
    std::vector<double> cutVar_TPCnSigmaEl_low = {-4., -4., -4., -2., -3., -2., -3., -4., -2., -4., -3., -3., -3., -4., -3., -3., -4., -3., -4., -4., -3., -4., -3., -3., -2., -4., -4., -3., -4., -2};
    std::vector<double> cutVar_TPCnSigmaEl_up = {2., 3., 2., 2., 2., 2., 2., 4., 2., 2., 4., 3., 3., 2., 4., 2., 2., 4., 3., 4., 2., 4., 2., 3., 2., 3., 4., 2., 3., 2};
    std::vector<double> cutVar_TPCnSigmaPi_low = {-3., -2., -3., -4., -4., -3., -4., -2., -2., -2., -3., -3., -2., -2., -4., -3., -3., -2., -3., -2., -4., -2., -4., -4., -3., -3., -3., -2., -4., -4};
    std::vector<double> cutVar_TPCnSigmaPi_up = {3., 3., 4., 4., 3., 2., 4., 4., 3., 4., 4., 3., 4., 4., 3., 2., 4., 2., 4., 2., 3., 4., 2., 2., 3., 2., 3., 4., 2., 4};
    std::vector<double> cutVar_TPCnSigmaKa_low = {-4., -2., -2., -2., -4., -3., -2., -4., -3., -3., -4., -2., -3., -3., -4., -2., -4., -2., -3., -4., -4., -2., -2., -3., -2., -2., -3., -3., -2., -4};
    std::vector<double> cutVar_TPCnSigmaKa_up = {4., 3., 2., 3., 4., 3., 4., 4., 4., 4., 4., 4., 4., 4., 2., 4., 4., 2., 2., 4., 3., 3., 2., 4., 2., 4., 3., 3., 3., 3};
    std::vector<double> cutVar_TPCnSigmaPr_low = {-4., -2., -2., -3., -4., -4., -3., -2., -2., -4., -4., -2., -3., -4., -2., -3., -3., -2., -3., -3., -2., -2., -2., -2., -2., -3., -2., -3., -3., -3};
    std::vector<double> cutVar_TPCnSigmaPr_up = {2., 2., 3., 2., 3., 3., 3., 2., 4., 3., 3., 4., 4., 3., 4., 4., 3., 4., 2., 3., 4., 4., 3., 4., 3., 2., 3., 3., 2., 3};

    for (unsigned int i = 0; i < cutVar_TPCnSigmaEl_low.size(); i++) {
      if (nameStr == Form("electronPID_TPCnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
          cut->AddCut(VarManager::kTPCnSigmaEl, cutVar_TPCnSigmaEl_low.at(i), cutVar_TPCnSigmaEl_up.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPi, cutVar_TPCnSigmaPi_low.at(i), cutVar_TPCnSigmaPi_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaKa, cutVar_TPCnSigmaKa_low.at(i), cutVar_TPCnSigmaKa_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPr, cutVar_TPCnSigmaPr_low.at(i), cutVar_TPCnSigmaPr_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
          cut->AddCut(VarManager::kTPCnSigmaEl_Corr, cutVar_TPCnSigmaEl_low.at(i), cutVar_TPCnSigmaEl_up.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPi_Corr, cutVar_TPCnSigmaPi_low.at(i), cutVar_TPCnSigmaPi_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaKa, cutVar_TPCnSigmaKa_low.at(i), cutVar_TPCnSigmaKa_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPr_Corr, cutVar_TPCnSigmaPr_low.at(i), cutVar_TPCnSigmaPr_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
          cut->AddCut(VarManager::kTPCnSigmaEl_Corr, cutVar_TPCnSigmaEl_low.at(i), cutVar_TPCnSigmaEl_up.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPi_Corr, cutVar_TPCnSigmaPi_low.at(i), cutVar_TPCnSigmaPi_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaKa_Corr, cutVar_TPCnSigmaKa_low.at(i), cutVar_TPCnSigmaKa_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPr_Corr, cutVar_TPCnSigmaPr_low.at(i), cutVar_TPCnSigmaPr_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
        }
        return cut;
      }
    }

    if (nameStr == Form("lmee_pp_502TeV_TPC%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_lowB_TPC%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.5, 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCloosenopkrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) {
        cut->AddCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) {
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCPbPbnopkrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TPCloose%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -2., 2., true, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }
  }

  if (nameStr == "electronPIDnsigmaOpen") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.0, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaVeryVeryLoose") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.0, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaVeryLoose") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.4, 1e+10, false);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaLoose") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaMedium") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.7, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.7, 3000.0);
    return cut;
  }
  if (nameStr == "electronPIDnsigmaMedium_withLargeTOFPID") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 2.7, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 2.7, 3000.0);
    cut->AddCut(VarManager::kTOFnSigmaEl, -5.0, 5.0);
    return cut;
  }
  if (nameStr == "electronPIDnsigmaSkewed") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -2.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaSkewed_2") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -0.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRej") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRej_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRejLoose") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDPrKaPiRejLoose_Corr") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "electronPIDnsigmaEMu") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -1.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.5, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaKa, 3.5, 3000.0);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma") {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "kaonRejNsigma") {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma2") {
    cut->AddCut(VarManager::kTPCnSigmaKa, -2.0, 2.0);
    return cut;
  }

  if (nameStr == "kaonPID_TPCnTOF") {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    cut->AddCut(VarManager::kTOFnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "kaonPIDnsigma700") {
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0);
    cut->AddCut(VarManager::kPin, 0.0, 0.7);
    return cut;
  }

  if (nameStr == "AssocKine") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (nameStr == "electronPIDworseRes") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0 * 0.8, 3000.0); // emulates a 20% degradation in PID resolution
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0 * 0.8, 3000.0); // proton and pion rejections are effectively relaxed by 20%
    return cut;
  }

  if (nameStr == "electronPIDshift") {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0 - 0.2, 3000.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0 - 0.2, 3000.0);
    return cut;
  }

  if (nameStr == "pionPIDnsigma") {
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "pionPID_TPCnTOF") {
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0);
    cut->AddCut(VarManager::kTOFnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "protonPID_TPCnTOF") {
    cut->AddCut(VarManager::kTPCnSigmaPr, -3.0, 3.0);
    cut->AddCut(VarManager::kTOFnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "protonPID_TPCnTOF2") {
    cut->AddCut(VarManager::kTPCnSigmaPr, -2.5, 2.5);
    return cut;
  }

  if (nameStr == "tpc_pion_rejection") {
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+[1]*x", 0, 10);
    f1maxPi->SetParameters(85, -50);
    cut->AddCut(VarManager::kTPCsignal, 70, f1maxPi, true, VarManager::kPin, 0.0, 0.4, false);
    return cut;
  }

  if (nameStr == "tpc_pion_band_rejection") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+[1]*log(x)", 0, 10);
    f1minPi->SetParameters(-115, -90);
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+[1]*log(x)", 0, 10);
    f1maxPi->SetParameters(-70, -90);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.05, 0.3, false);
    return cut;
  }

  if (nameStr == "tpc_pion_muon_band_rejection") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPi->SetParameters(37, -18, 5.5);
    std::shared_ptr<TF1> f1maxPi = std::make_shared<TF1>("f1maxPi", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPi->SetParameters(67, -17, 5.9);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, f1maxPi, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_pion_rejection_highp") {
    std::shared_ptr<TF1> f1minPi = std::make_shared<TF1>("f1minPi", "[0]+[1]*x", 0, 10);
    f1minPi->SetParameters(65, 4.);
    cut->AddCut(VarManager::kTPCsignal, f1minPi, 110., false, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_kaon_rejection") {
    std::shared_ptr<TF1> f1minKa = std::make_shared<TF1>("f1minKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1minKa->SetParameters(37, -4, 5.6);
    std::shared_ptr<TF1> f1maxKa = std::make_shared<TF1>("f1maxKa", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxKa->SetParameters(60, -4.1, 6.);
    cut->AddCut(VarManager::kTPCsignal, f1minKa, f1maxKa, true, VarManager::kPin, 0.0, 10.0, false);
    return cut;
  }

  if (nameStr == "tpc_proton_rejection") {
    std::shared_ptr<TF1> f1minPr = std::make_shared<TF1>("f1minPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1minPr->SetParameters(37, -2.6, 6.1);
    std::shared_ptr<TF1> f1maxPr = std::make_shared<TF1>("f1maxPr", "[0]+exp([1]*x+[2])", 0, 10);
    f1maxPr->SetParameters(60, -2.4, 6.2);
    cut->AddCut(VarManager::kTPCsignal, f1minPr, f1maxPr, true, VarManager::kPin, 0.0, 10, false);
    return cut;
  }

  if (nameStr == "tpc_electron") {
    cut->AddCut(VarManager::kTPCsignal, 60, 110, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron") {
    cut->AddCut(VarManager::kTOFbeta, 0.99, 1.01, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron_sigma") {
    cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (nameStr == "tof_electron_sigma_2") {
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3.);
    return cut;
  }

  if (nameStr == "tof_electron_loose") {
    cut->AddCut(VarManager::kTOFbeta, 0.95, 1.05, false, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  // loop to define TOF PID cuts with and without post calibration
  for (size_t icase = 0; icase < vecPIDcase.size(); icase++) {
    if (nameStr == Form("electronPID_TOFnsigma%s_loose", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TOFnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongHadRej", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    // List of nSigma values for lower and upper edge of TPC: El, Pi and TOF: El PID selections/rejection
    std::vector<double> cutVar_TPCnSigmaEl_low = {-4., -4., -4., -2., -3., -2., -3., -4., -2., -4., -3., -3., -3., -4., -3., -3., -4., -3., -4., -4., -3., -4., -3., -3., -2., -4., -4., -3., -4., -2};
    std::vector<double> cutVar_TPCnSigmaEl_up = {2., 3., 2., 2., 2., 2., 2., 4., 2., 2., 4., 3., 3., 2., 4., 2., 2., 4., 3., 4., 2., 4., 2., 3., 2., 3., 4., 2., 3., 2};
    std::vector<double> cutVar_TPCnSigmaPi_low = {-3., -2., -3., -4., -4., -3., -4., -2., -2., -2., -3., -3., -2., -2., -4., -3., -3., -2., -3., -2., -4., -2., -4., -4., -3., -3., -3., -2., -4., -4};
    std::vector<double> cutVar_TPCnSigmaPi_up = {3., 3., 4., 4., 3., 2., 4., 4., 3., 4., 4., 3., 4., 4., 3., 2., 4., 2., 4., 2., 3., 4., 2., 2., 3., 2., 3., 4., 2., 4};
    std::vector<double> cutVar_TOFnSigmaEl_low = {-4., -2., -4., -4., -3., -2., -4., -4., -4., -2., -2., -4., -3., -3., -4., -4., -4., -2., -4., -4., -2., -2., -3., -4., -4., -2., -4., -2., -3., -3};
    std::vector<double> cutVar_TOFnSigmaEl_up = {4., 2., 4., 2., 4., 3., 2., 3., 3., 3., 4., 3., 2., 3., 4., 3., 3., 3., 4., 4., 2., 2., 2., 3., 3., 3., 2., 3., 2., 4};

    for (unsigned int i = 0; i < cutVar_TOFnSigmaEl_low.size(); i++) {
      if (nameStr == Form("electronPID_TOFnsigma_cutVar%s%i", vecPIDcase.at(icase).Data(), i)) {
        if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
          cut->AddCut(VarManager::kTPCnSigmaEl, cutVar_TPCnSigmaEl_low.at(i), cutVar_TPCnSigmaEl_up.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPi, cutVar_TPCnSigmaPi_low.at(i), cutVar_TPCnSigmaPi_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTOFnSigmaEl, cutVar_TOFnSigmaEl_low.at(i), cutVar_TOFnSigmaEl_up.at(i), false, VarManager::kPin, 0.3, 1e+10, false);
        } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
          cut->AddCut(VarManager::kTPCnSigmaEl_Corr, cutVar_TPCnSigmaEl_low.at(i), cutVar_TPCnSigmaEl_up.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTPCnSigmaPi_Corr, cutVar_TPCnSigmaPi_low.at(i), cutVar_TPCnSigmaPi_up.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
          cut->AddCut(VarManager::kTOFnSigmaEl, cutVar_TOFnSigmaEl_low.at(i), cutVar_TOFnSigmaEl_up.at(i), false, VarManager::kPin, 0.3, 1e+10, false);
        }
        return cut;
      }
    }

    if (nameStr == Form("electronPID_TPC_TOFnsigma%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // previously known as electronPID_TOFnsigma_tight  // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TOFnsigma%s_strongNSigE", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_lowB_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.3, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 4., true, VarManager::kPin, 0.25, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0.0, 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigE_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFnsigma%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      cut->AddCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cut->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
      return cut;
    }

    if (nameStr == Form("electronPID_TOFreq%s_strongNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -1., 2., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -2., 2., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("electronPID_TOFreq%s_tightNSigEPbPb_rejBadTOF", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, 0., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 1., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOF%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_lowB_TOF%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaMu, -3., 3.5, true, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOFloose%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 3., true, VarManager::kPin, 0.0, 4.0, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 4.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 3., true, VarManager::kPin, 0.0, 4.0, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 4.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }

    if (nameStr == Form("lmee_pp_502TeV_TOFloose_pionrej%s", vecPIDcase.at(icase).Data())) {
      if (icase == 0) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 3.5, true, VarManager::kPin, 0.0, 2.0, false);
        cut->AddCut(VarManager::kTPCnSigmaPi, -99., 2.5, true, VarManager::kPin, 2.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      } else if (icase == 1 || icase == 2) { // o2-linter: disable=magic-number (number of cuts)
        cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4., 4., false, VarManager::kPin, 0.0, 1e+10, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 3.5, true, VarManager::kPin, 0.0, 2.0, false);
        cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -99., 2.5, true, VarManager::kPin, 2.0, 1e+10, false);
        cut->AddCut(VarManager::kTOFnSigmaEl, -4., 4., false, VarManager::kPin, 0.3, 1e+10, false);
      }
      return cut;
    }
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts
  if (nameStr == "GlobalMuonTrack") {
    cut->AddCut(VarManager::kMuonTrackType, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "MFTMCH") {
    cut->AddCut(VarManager::kMuonTrackType, 1.5, 2.5);
    return cut;
  }

  if (nameStr == "MCHMID") {
    cut->AddCut(VarManager::kMuonTrackType, 2.5, 3.5);
    return cut;
  }

  if (nameStr == "MCHStandalone") {
    cut->AddCut(VarManager::kMuonTrackType, 3.5, 4.5);
    return cut;
  }

  if (nameStr == "muonMinimalCuts") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    return cut;
  }

  if (nameStr == "muonMinimalCuts10SigmaPDCA") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 990.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 540.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    return cut;
  }

  if (nameStr == "muonQualityCuts") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonQualityCuts5SigmaPDCA_Run3") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 500.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 335.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonQualityCuts10SigmaPDCA") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 990.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 540.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "matchedQualityCuts") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "matchedQualityCutsMFTeta") {
    cut->AddCut(VarManager::kEta, -3.6, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "muonQualityCutsMatchingOnly") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonLowPt") {
    cut->AddCut(VarManager::kPt, 0.5, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt2") {
    cut->AddCut(VarManager::kPt, 0.7, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt3") {
    cut->AddCut(VarManager::kPt, 0.8, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt4") {
    cut->AddCut(VarManager::kPt, 0.9, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt5") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonLowPt6") {
    cut->AddCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt") {
    cut->AddCut(VarManager::kPt, 3.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt2") {
    cut->AddCut(VarManager::kPt, 4.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt3") {
    cut->AddCut(VarManager::kPt, 6.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt4") {
    cut->AddCut(VarManager::kPt, 8.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt5") {
    cut->AddCut(VarManager::kPt, 10.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonHighPt6") {
    cut->AddCut(VarManager::kPt, 20.0, 1000.0);
    return cut;
  }

  if (nameStr == "muonTightQualityCutsForTests") {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 20.0, 60.0);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    return cut;
  }

  if (nameStr == "mchTrack") {
    cut->AddCut(VarManager::kMuonTrackType, 3.5, 4.5);
    return cut;
  }

  if (nameStr == "matchedMchMid") {
    cut->AddCut(VarManager::kMuonTrackType, 2.5, 3.5);
    return cut;
  }

  if (nameStr == "matchedFwd") {
    cut->AddCut(VarManager::kMuonTrackType, 1.5, 2.5);
    return cut;
  }

  if (nameStr == "matchedGlobal") {
    cut->AddCut(VarManager::kMuonTrackType, -0.5, 0.5);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut1") {
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0, 30);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut2") {
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0, 40);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut3") {
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0, 50);
    return cut;
  }

  if (nameStr == "Chi2MCHMFTCut4") {
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0, 60);
    return cut;
  }

  // -----------------------------------------------------------------------------------------------
  // Pair cuts
  if (nameStr == "pairDalitz1") {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    std::shared_ptr<TF1> fcutHigh = std::make_shared<TF1>("f1", "[0] - [0]/[1]*x", -1.5, 1.5);
    fcutHigh->SetParameters(0.6, 0.12);
    std::shared_ptr<TF1> fcutLow = std::make_shared<TF1>("f2", "-[0] + [0]/[1]*x", -1.5, 1.5);
    fcutLow->SetParameters(0.6, 0.12);
    cut->AddCut(VarManager::kPsiPair, fcutLow, fcutHigh, true, VarManager::kDeltaPhiPair, 0, 0.12);
    return cut;
  }

  if (nameStr == "pairDalitz1Strong") {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    cut->AddCut(VarManager::kDeltaPhiPair, -1., 0.);
    return cut;
  }

  if (nameStr == "pairDalitz2") {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    return cut;
  }

  if (nameStr == "pairDalitz3") {
    cut->AddCut(VarManager::kMass, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "paira_prefilter1") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    return cut;
  }

  if (nameStr == "paira_prefilter2") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "paira_prefilter3") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "paira_prefilter4") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "paira_prefilter5") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "paira_prefilter6") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "paira_prefilter7") {
    cut->AddCut(VarManager::kMass, 0.0, 0.06);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairb_prefilter1") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairb_prefilter2") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "pairb_prefilter3") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairb_prefilter4") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairb_prefilter5") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairb_prefilter6") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairb_prefilter7") {
    cut->AddCut(VarManager::kMass, 0.0, 0.05);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairc_prefilter1") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    return cut;
  }

  if (nameStr == "pairc_prefilter2") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "pairc_prefilter3") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairc_prefilter4") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairc_prefilter5") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairc_prefilter6") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairc_prefilter7") {
    cut->AddCut(VarManager::kMass, 0.0, 0.04);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "paird_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.025);
    return cut;
  }

  if (nameStr == "paire_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.05);
    return cut;
  }

  if (nameStr == "pairf_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.075);
    return cut;
  }

  if (nameStr == "pairg_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.1);
    return cut;
  }

  if (nameStr == "pairh_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.125);
    return cut;
  }

  if (nameStr == "pairi_prefilter1") {
    cut->AddCut(VarManager::kOpeningAngle, 0.0, 0.15);
    return cut;
  }

  if (nameStr == "pairNoCut") {
    cut->AddCut(VarManager::kMass, 0.0, 1000.0);
    return cut;
  }

  if (nameStr == "DipionMassCut1") {
    cut->AddCut(VarManager::kMass, 0.5, 1.0);
    return cut;
  }

  if (nameStr == "DipionMassCut2") {
    cut->AddCut(VarManager::kMass, 0.0, 1.0);
    return cut;
  }

  if (nameStr == "pairMassLow1") {
    cut->AddCut(VarManager::kMass, 1.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow2") {
    cut->AddCut(VarManager::kMass, 1.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow3") {
    cut->AddCut(VarManager::kMass, 1.6, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow4") {
    cut->AddCut(VarManager::kMass, 1.7, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow5") {
    cut->AddCut(VarManager::kMass, 1.8, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow6") {
    cut->AddCut(VarManager::kMass, 1.85, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow7") {
    cut->AddCut(VarManager::kMass, 1.9, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow8") {
    cut->AddCut(VarManager::kMass, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow9") {
    cut->AddCut(VarManager::kMass, 2.2, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow10") {
    cut->AddCut(VarManager::kMass, 2.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow11") {
    cut->AddCut(VarManager::kMass, 3.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairMassLow12") {
    cut->AddCut(VarManager::kMass, 3.5, 1000.0);
    return cut;
  }

  if (nameStr == "pairMass1to2") {
    cut->AddCut(VarManager::kMass, 1., 2.);
    return cut;
  }

  if (nameStr == "pairMassIMR") {
    cut->AddCut(VarManager::kMass, 1.1, 2.7);
    return cut;
  }

  if (nameStr == "pairMass1_5to2_7") {
    cut->AddCut(VarManager::kMass, 1.5, 2.7);
    return cut;
  }

  if (nameStr == "pairMass1_3to3_5") {
    cut->AddCut(VarManager::kMass, 1.3, 3.5);
    return cut;
  }

  if (nameStr == "pairMass1_3") {
    cut->AddCut(VarManager::kMass, 1.3, 1000.0);
    return cut;
  }

  if (nameStr == "pairMass1_5to3_5") {
    cut->AddCut(VarManager::kMass, 1.5, 3.5);
    return cut;
  }

  if (nameStr == "pairD0") {
    cut->AddCut(VarManager::kMass, 1.814, 1.914);
    return cut;
  }

  if (nameStr == "pairJpsi") {
    cut->AddCut(VarManager::kMass, 2.8, 3.3);
    return cut;
  }

  if (nameStr == "pairJpsi2") {
    cut->AddCut(VarManager::kMass, 2.72, 3.2);
    return cut;
  }

  if (nameStr == "pairJpsi3") {
    cut->AddCut(VarManager::kMass, 2.92, 3.14);
    return cut;
  }

  if (nameStr == "pairPsi2S") {
    cut->AddCut(VarManager::kMass, 3.4, 3.9);
    return cut;
  }

  if (nameStr == "pairUpsilon") {
    cut->AddCut(VarManager::kMass, 8.0, 11.0);
    return cut;
  }

  if (nameStr == "pairX3872") {
    cut->AddCut(VarManager::kQ, 0.0, 0.3);
    return cut;
  }

  if (nameStr == "pairX3872_2") {
    cut->AddCut(VarManager::kQuadDefaultDileptonMass, 3.0, 5.0);
    cut->AddCut(VarManager::kQ, 0.0, 0.5);
    cut->AddCut(VarManager::kDeltaR, 0.0, 5.0);
    cut->AddCut(VarManager::kQuadPt, 5.0, 40.0);
    return cut;
  }

  if (nameStr == "pairX3872_3") {
    cut->AddCut(VarManager::kQuadDefaultDileptonMass, 3.0, 5.0);
    cut->AddCut(VarManager::kQ, 0.0, 0.5);
    cut->AddCut(VarManager::kDeltaR, 0.0, 5.0);
    cut->AddCut(VarManager::kQuadPt, 0.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairPtLow1") {
    cut->AddCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairPtLow2") {
    cut->AddCut(VarManager::kPt, 5.0, 1000.0);
    return cut;
  }

  if (nameStr == "pairPtLow3") {
    cut->AddCut(VarManager::kPt, 0, 0.15);
    return cut;
  }

  if (nameStr == "pairPtLow4") {
    cut->AddCut(VarManager::kPt, 0, 10.0);
    return cut;
  }

  if (nameStr == "pairPtLow5") {
    cut->AddCut(VarManager::kPt, 0.8, 1000.0);
    return cut;
  }

  if (nameStr == "pairRapidityForward") {
    cut->AddCut(VarManager::kRap, 2.5, 4.0);
    return cut;
  }

  if (nameStr == "pairDCA") {
    cut->AddCut(VarManager::kQuadDCAabsXY, .0, .50);
    return cut;
  }

  if (nameStr == "singleDCA") {
    cut->AddCut(VarManager::kTrackDCAsigXY, 0.0, 5.);
    return cut;
  }

  if (nameStr == "pairPhiV") {
    cut->AddCut(VarManager::kPairPhiv, 2., 3.2);
    return cut;
  }

  if (nameStr == "excludePairPhiV") {
    cut->AddCut(VarManager::kPairPhiv, 2., 3.2, true);
    return cut;
  }

  if (nameStr == "pairLowMass") {
    cut->AddCut(VarManager::kMass, 0., 0.1);
    return cut;
  }

  if (nameStr == "excludePairLowMass") {
    cut->AddCut(VarManager::kMass, 0., 0.1, true);
    return cut;
  }

  if (nameStr == "pairLxyzProjected3sigma") {
    cut->AddCut(VarManager::kVertexingLxyzProjected, 0.015, 10.);
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1") {
    cut->AddCut(VarManager::kVertexingTauxyzProjected, 0.0005, 10.);
    return cut;
  }

  if (nameStr == "pairTauxyzProjected1sigma") {
    cut->AddCut(VarManager::kVertexingTauxyzProjected, 0.003, 10.);
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaLambdacCand") {
    std::shared_ptr<TF1> f1minLxyProjected = std::make_shared<TF1>("f1minLxyProjected", "[0]+[1]*x", 0., 20.);
    f1minLxyProjected->SetParameters(0.0065, -0.00023);
    cut->AddCut(VarManager::kVertexingLxyProjected, f1minLxyProjected, 1., false, VarManager::kPt, 0., 20.);
    return cut;
  }

  if (nameStr == "pairLxyProjected3sigmaDplusCand") {
    cut->AddCut(VarManager::kVertexingLxyProjected, 0.009, 10.);
    return cut;
  }

  if (nameStr == "pairCosPointingPos") {
    cut->AddCut(VarManager::kCosPointingAngle, 0.9, 1000.);
    return cut;
  }

  if (nameStr == "pairCosPointingNeg90") {
    cut->AddCut(VarManager::kCosPointingAngle, -1000., -0.9);
    return cut;
  }

  if (nameStr == "pairCosPointingNeg85") {
    cut->AddCut(VarManager::kCosPointingAngle, -1000., -0.85);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and pair selection in order or optimize the trigger
  // trigger selection cuts

  if (nameStr == "electronStandardQualityTriggerTest") {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 50.0, 161.);
    return cut;
  }

  if (nameStr == "muonLooseTriggerTestCuts") {
    cut->AddCut(VarManager::kEta, -4.5, -2.0);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (nameStr == "muonMatchingMFTMCHTriggerTestCuts") {
    cut->AddCut(VarManager::kEta, -4.5, -2.0);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest1") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest2") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest3") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (nameStr == "jpsi_TPCPID_TriggerTest4") {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  //---------------------------------------------------------------
  // ALICE 3 Analysis Cuts

  if (nameStr == "alice3KineSkim") {
    cut->AddCut(VarManager::kPt, 0.05, 1000.0);
    cut->AddCut(VarManager::kEta, -4.0, 4.0);
    return cut;
  }

  if (nameStr == "alice3StandardKine") {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -2.5, 2.5); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineTOFAcceptance") {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -2.0, 2.0); // TOF acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3KineRICHAcceptance") {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8); // RICH acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3JpsiKine") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -2.5, 2.5); // Total tracker acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3JpsiKineTOFAcceptance") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -2.0, 2.0); // TOF acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3JpsiKineRICHAcceptance") {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8); // RICH acceptance in v3b geomety
    return cut;
  }

  if (nameStr == "alice3TrackQuality") {
    cut->AddCut(VarManager::kIsReconstructed, 0.5, 1.5);
    cut->AddCut(VarManager::kNSiliconHits, 6.0, 12.0);
    cut->AddCut(VarManager::kTrackDCAxy, -3.0, 3.0);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3TrackQualityTightDCA") {
    cut->AddCut(VarManager::kIsReconstructed, 0.5, 1.5);
    cut->AddCut(VarManager::kNSiliconHits, 6.0, 12.0);
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->AddCut(VarManager::kTrackDCAz, -1.0, 1.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDEl") {
    cut->AddCut(VarManager::kOTnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDPi") {
    cut->AddCut(VarManager::kOTnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDKa") {
    cut->AddCut(VarManager::kOTnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3OTPIDPr") {
    cut->AddCut(VarManager::kOTnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDEl") {
    cut->AddCut(VarManager::kInnerTOFnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPi") {
    cut->AddCut(VarManager::kInnerTOFnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDKa") {
    cut->AddCut(VarManager::kInnerTOFnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3iTOFPIDPr") {
    cut->AddCut(VarManager::kInnerTOFnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDEl") {
    cut->AddCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDEl") {
    cut->AddCut(VarManager::kOuterTOFnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPi") {
    cut->AddCut(VarManager::kOuterTOFnSigmaPi, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDKa") {
    cut->AddCut(VarManager::kOuterTOFnSigmaKa, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3oTOFPIDPr") {
    cut->AddCut(VarManager::kOuterTOFnSigmaPr, -3.0, 3.0);
    return cut;
  }

  if (nameStr == "alice3RICHPIDEl") {
    cut->AddCut(VarManager::kRICHnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kHasRICHSigEl, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDPi") {
    cut->AddCut(VarManager::kRICHnSigmaPi, -3.0, 3.0);
    cut->AddCut(VarManager::kHasRICHSigPi, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDKa") {
    cut->AddCut(VarManager::kRICHnSigmaKa, -3.0, 3.0);
    cut->AddCut(VarManager::kHasRICHSigKa, 0.5, 1.5);
    return cut;
  }

  if (nameStr == "alice3RICHPIDPr") {
    cut->AddCut(VarManager::kRICHnSigmaPr, -3.0, 3.0);
    cut->AddCut(VarManager::kHasRICHSigPr, 0.5, 1.5);
    return cut;
  }

  delete cut;
  LOGF(fatal, Form("Did not find cut %s", cutName));
  return nullptr;
}

//________________________________________________________________________________________________
std::vector<AnalysisCut*> o2::aod::dqcuts::GetCutsFromJSON(const char* json)
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
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")";
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
      AnalysisCut* anaCut = ParseJSONAnalysisCut(&cut, cutName);
      if (anaCut != nullptr) {
        cuts.push_back(anaCut);
      } else {
        LOG(fatal) << "Something went wrong in creating the AnalysisCut " << cutName;
        return cuts;
      }
    }
    if (isAnalysisCompositeCut) {
      AnalysisCompositeCut* anaCut = ParseJSONAnalysisCompositeCut(&cut, cutName);
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
bool o2::aod::dqcuts::ValidateJSONAnalysisCut(T cut)
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
bool o2::aod::dqcuts::ValidateJSONAddCut(T addcut, bool isSimple)
{
  //
  // Validate AddCut definition in JSON file
  //

  // Check if this AddCut is adding an analysis cut (if the mother is a composite cut)
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
    LOG(fatal) << "One cannot call AddCut with an AnalysisCut in this case";
    return false;
  }
  if (!isSimple && isAnalysisCut) {
    return true;
  }

  // check that the variable to cut on is a valid one (implemented in the VarManager)
  const char* var = addcut->FindMember("var")->value.GetString();
  if (VarManager::fgVarNamesMap.find(var) == VarManager::fgVarNamesMap.end()) {
    LOG(fatal) << "Bad cut variable (" << var << ") specified for this AddCut";
    return false;
  }

  // check whether the specified cut low and high are numbers or functions
  bool cutLow_isNumber = addcut->FindMember("cutLow")->value.IsNumber();
  bool cutHigh_isNumber = addcut->FindMember("cutHigh")->value.IsNumber();
  if (!cutLow_isNumber) {
    auto& cutLow = addcut->FindMember("cutLow")->value;
    if (!cutLow.HasMember("funcName") || !cutLow.HasMember("funcBody") ||
        !cutLow.HasMember("xLow") || !cutLow.HasMember("xHigh")) {
      LOG(fatal) << "Missing fields for the cutLow TF1 definition";
      return false;
    }
  }
  if (!cutHigh_isNumber) {
    auto& cutHigh = addcut->FindMember("cutHigh")->value;
    if (!cutHigh.HasMember("funcName") || !cutHigh.HasMember("funcBody") ||
        !cutHigh.HasMember("xLow") || !cutHigh.HasMember("xHigh")) {
      LOG(fatal) << "Missing fields for the cutLow TF1 definition";
      return false;
    }
  }
  if (!cutHigh_isNumber || !cutLow_isNumber) {
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
  if (addcut->HasMember("dependentVar1") && cutLow_isNumber && cutHigh_isNumber) {
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
AnalysisCut* o2::aod::dqcuts::ParseJSONAnalysisCut(T cut, const char* cutName)
{

  // Parse the json object and build an AnalysisCut
  if (!ValidateJSONAnalysisCut(cut)) {
    LOG(fatal) << "AnalysisCut not properly defined in the JSON file. Skipping it";
    return nullptr;
  }

  // If the analysis cut has the field "library", its just loaded from the preexisting cuts in the library and return
  if (cut->HasMember("library")) {
    return GetAnalysisCut(cut->FindMember("library")->value.GetString());
  }

  // construct the AnalysisCut object and add the AddCuts
  auto* retCut = new AnalysisCut(cutName, cut->HasMember("library") ? cut->FindMember("title")->value.GetString() : "");

  // loop over all the members for this cut and configure the AddCut objects
  for (rapidjson::Value::ConstMemberIterator it = cut->MemberBegin(); it != cut->MemberEnd(); it++) {

    TString itName = it->name.GetString();

    if (itName.Contains("AddCut")) {

      LOG(debug) << "Parsing " << itName.Data();
      const auto& addcut = it->value;
      if (!ValidateJSONAddCut(&addcut, true)) {
        LOG(fatal) << "AddCut statement not properly defined";
        return nullptr;
      }

      const char* var = addcut.FindMember("var")->value.GetString();
      LOG(info) << "var " << var;
      bool cutLow_isNumber = addcut.FindMember("cutLow")->value.IsNumber();
      LOG(info) << "cutLow_isNumber " << cutLow_isNumber;
      bool cutHigh_isNumber = addcut.FindMember("cutHigh")->value.IsNumber();
      LOG(info) << "cutHigh_isNumber " << cutHigh_isNumber;

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
      if (cutLow_isNumber) {
        cutLowNumber = addcut.FindMember("cutLow")->value.GetDouble();
        LOG(info) << "cutLowNumber " << cutLowNumber;
      } else {
        auto& cutLow = addcut.FindMember("cutLow")->value;
        cutLowFunc = std::make_shared<TF1>(cutLow.FindMember("funcName")->value.GetString(), cutLow.FindMember("funcBody")->value.GetString(),
                                           cutLow.FindMember("xLow")->value.GetDouble(), cutLow.FindMember("xHigh")->value.GetDouble());
        LOG(info) << "cutLowFunc " << cutLow.FindMember("funcName")->value.GetString() << ", " << cutLow.FindMember("funcBody")->value.GetString()
                  << ", " << cutLow.FindMember("xLow")->value.GetDouble() << ", " << cutLow.FindMember("xHigh")->value.GetDouble();
      }
      if (cutHigh_isNumber) {
        cutHighNumber = addcut.FindMember("cutHigh")->value.GetDouble();
        LOG(info) << "cutHighNumber " << cutHighNumber;
      } else {
        auto& cutHigh = addcut.FindMember("cutHigh")->value;
        cutHighFunc = std::make_shared<TF1>(cutHigh.FindMember("funcName")->value.GetString(), cutHigh.FindMember("funcBody")->value.GetString(),
                                            cutHigh.FindMember("xLow")->value.GetDouble(), cutHigh.FindMember("xHigh")->value.GetDouble());
        LOG(info) << "cutHighFunc " << cutHigh.FindMember("funcName")->value.GetString() << ", " << cutHigh.FindMember("funcBody")->value.GetString()
                  << ", " << cutHigh.FindMember("xLow")->value.GetDouble() << ", " << cutHigh.FindMember("xHigh")->value.GetDouble();
      }
      if (cutLow_isNumber) {
        if (cutHigh_isNumber) {
          retCut->AddCut(VarManager::fgVarNamesMap[var], cutLowNumber, cutHighNumber, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        } else {
          retCut->AddCut(VarManager::fgVarNamesMap[var], cutLowNumber, cutHighFunc, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        }
      } else {
        if (cutHigh_isNumber) {
          retCut->AddCut(VarManager::fgVarNamesMap[var], cutLowFunc, cutHighNumber, exclude,
                         VarManager::fgVarNamesMap[dependentVar], depCutLow, depCutHigh, depCutExclude,
                         VarManager::fgVarNamesMap[dependentVar2], depCut2Low, depCut2High, depCut2Exclude);
        } else {
          retCut->AddCut(VarManager::fgVarNamesMap[var], cutLowFunc, cutHighFunc, exclude,
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
bool o2::aod::dqcuts::ValidateJSONAnalysisCompositeCut(T cut)
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
AnalysisCompositeCut* o2::aod::dqcuts::ParseJSONAnalysisCompositeCut(T cut, const char* cutName)
{

  // Configure an AnalysisCompositeCut
  if (!ValidateJSONAnalysisCompositeCut(cut)) {
    LOG(fatal) << "Composite Cut not properly defined in the JSON file. Skipping it";
    return nullptr;
  }

  if (cut->HasMember("library")) {
    return GetCompositeCut(cut->FindMember("library")->value.GetString());
  }

  auto* retCut = new AnalysisCompositeCut(cutName, cut->HasMember("library") ? cut->FindMember("title")->value.GetString() : "", cut->FindMember("useAND")->value.GetBool());

  // Loop to find AddCut objects
  for (rapidjson::Value::ConstMemberIterator it = cut->MemberBegin(); it != cut->MemberEnd(); it++) {

    TString itName = it->name.GetString();

    if (itName.Contains("AddCut")) {

      LOG(debug) << "Parsing " << itName.Data();
      const auto& addcut = it->value;
      if (!ValidateJSONAddCut(&addcut, false)) {
        LOG(fatal) << "AddCut statement not properly defined";
        return nullptr;
      }

      // For an AnalysisCompositeCut, one can call AddCut with either an AnalysisCut or another AnalysisCompositeCut
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
        AnalysisCut* cutMember = ParseJSONAnalysisCut(&addcut, itName.Data());
        if (cutMember != nullptr) {
          retCut->AddCut(cutMember);
        } else {
          return nullptr;
        }
      }
      // Add an AnalysisCompositeCut
      if (isAnalysisCompositeCut) {
        AnalysisCompositeCut* cutMember = ParseJSONAnalysisCompositeCut(&addcut, itName.Data());
        if (cutMember != nullptr) {
          retCut->AddCut(cutMember);
        } else {
          return nullptr;
        }
      }
    }
  }

  return retCut;
}

//________________________________________________________________________________________________
o2::aod::dqmlcuts::BdtScoreConfig o2::aod::dqmlcuts::GetBdtScoreCutsAndConfigFromJSON(const char* json)
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
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")";
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

        if (!cutDirsFilled) {
          cutDirsFilled = true;
        }

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

      // MultiClass
    } else if (typeStr == "MultiClass") {
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
