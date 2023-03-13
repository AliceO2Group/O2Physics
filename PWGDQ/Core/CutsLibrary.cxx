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
#include "PWGDQ/Core/CutsLibrary.h"

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

  if (!nameStr.compare("jpsiO2MCdebugCuts2")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts2_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts2_prefiltered1")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    cut->AddCut(GetAnalysisCut("notDalitzLeg1"));
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

  if (!nameStr.compare("jpsiO2MCdebugCuts4_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
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

  if (!nameStr.compare("jpsiO2MCdebugCuts7_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts7_Corr_2")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine2"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts7_Corr_3")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine3"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));

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

  if (!nameStr.compare("jpsiO2MCdebugCuts10_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts10_Corr_Amb")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    cut->AddCut(GetAnalysisCut("ambiguousTrack")); // IsAmbiguous
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts11_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug3")); // cut on 1 ITS cluster
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts12")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaVeryLoose"));     // with 3 sigma El TOF
    return cut;
  }

  if (!nameStr.compare("jpsiO2MCdebugCuts13_Corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly")); // no cut on ITS clusters
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrackDCA")); // with DCA cut
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

  if (!nameStr.compare("lowMultTrackCut")) {
    cut->AddCut(GetAnalysisCut("lowMultTrackCut"));
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

  if (!nameStr.compare("kaonPID2")) {
    cut->AddCut(GetAnalysisCut("PIDStandardKine")); // standard kine cuts usually are applied via Filter in the task
    cut->AddCut(GetAnalysisCut("electronStandardQualityForO2MCdebug"));
    cut->AddCut(GetAnalysisCut("kaonPIDnsigma2"));
    return cut;
  }
  // NOTE Below there are several TPC pid cuts used for studies of the Run3 TPC post PID calib.
  if (!nameStr.compare("Jpsi_TPCPost_calib_debug1")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (!nameStr.compare("Jpsi_TPCPost_calib_debug2")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (!nameStr.compare("Jpsi_TPCPost_calib_noITSCuts_debug2")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }
  if (!nameStr.compare("Jpsi_TPCPost_calib_debug3")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug3"));
    return cut;
  }
  if (!nameStr.compare("Jpsi_TPCPost_calib_debug4")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }
  if (!nameStr.compare("Jpsi_TPCPost_calib_noITSCuts_debug4")) {
    cut->AddCut(GetAnalysisCut("jpsi_trackCut_noITSCuts_debug"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug4"));
    return cut;
  }

  if (!nameStr.compare("LMee_TPCPost_calib_debug1")) {
    cut->AddCut(GetAnalysisCut("lmee_trackCut_debug"));
    cut->AddCut(GetAnalysisCut("lmee_TPCPID_debug1"));
    return cut;
  }

  //---------------------------------------------------------------
  // Cuts for the selection of legs from dalitz decay
  //
  if (!nameStr.compare("DalitzCut1")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2_Corr")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut3")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut3_Corr")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut1SPDfirst")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut1SPDfirst_Corr")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDOnly_Corr"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2SPDfirst")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut2SPDfirst_Corr")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej_Corr"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut3SPDfirst")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose"));
    return cut;
  }

  if (!nameStr.compare("DalitzCut3SPDfirst_Corr")) {
    cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
    cut->AddCut(GetAnalysisCut("SPDfirst"));
    cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRejLoose_Corr"));
    return cut;
  }

  for (int i = 1; i <= 8; i++) {
    if (!nameStr.compare(Form("DalitzCut1_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("electronPIDOnly"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("DalitzCut2_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("DalitzCut2_Corr_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("DalitzCut1SPDfirst_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("SPDfirst"));
      cut->AddCut(GetAnalysisCut("electronPIDOnly"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("DalitzCut2SPDfirst_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("SPDfirst"));
      cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("DalitzCut2SPDfirst_Corr_dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut("dalitzStandardKine"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityITSOnly"));
      cut->AddCut(GetAnalysisCut("electronStandardQualityTPCOnly"));
      cut->AddCut(GetAnalysisCut("SPDfirst"));
      cut->AddCut(GetAnalysisCut("electronPIDPrKaPiRej"));
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }

    if (!nameStr.compare(Form("dalitzSelected%d", i))) {
      cut->AddCut(GetAnalysisCut(Form("dalitzLeg%d", i)));
      return cut;
    }
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
  //
  // LMee cuts
  // List of cuts used for low mass dielectron analyses
  //
  // Skimming cuts:
  if (!nameStr.compare("lmee_skimming")) {
    cut->AddCut(GetAnalysisCut("lmee_skimming_cuts"));
    return cut;
  }

  // LMee Run2 PID cuts

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

  if (!nameStr.compare("lmee_eNSigmaRun2")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
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

  if (!nameStr.compare("lmee_GlobalTrackRun2")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun2_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_TPCTrackRun2_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_TPCTrackRun2")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrack"));
    cut->AddCut(GetAnalysisCut("standardPrimaryTrack"));
    return cut;
  }

  // LMee Run3 PID cuts

  if (!nameStr.compare("lmeePID_TPChadrejTOFrecRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

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

  if (!nameStr.compare("lmee_eNSigmaRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

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

  if (!nameStr.compare("lmee_eNSigmaRun3_corr")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma_corr", "pid_TPCnSigma_corr", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_corr"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma_corr", "e_NSigma_corr", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_eNSigmaRun3_tight")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_eNSigmaRun3_strongPID_looseDCAxy")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCAxy"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_eNSigmaRun3_tight_pt04")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine_pt04"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  // 4 cuts to separate pos & neg tracks in pos & neg eta range
  if (!nameStr.compare("lmee_posTrack_posEta_selection")) {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("posEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_negTrack_posEta_selection")) {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("posEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_posTrack_negEta_selection")) {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("negEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_negTrack_negEta_selection")) {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("negEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  // 4 cuts to separate pos & neg tracks in pos & neg eta range applying electron PID
  if (!nameStr.compare("lmee_posNSigmaRun3_posEta_tight")) {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("posEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_negNSigmaRun3_posEta_tight")) {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("posEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_posNSigmaRun3_negEta_tight")) {
    cut->AddCut(GetAnalysisCut("posTrack"));
    cut->AddCut(GetAnalysisCut("negEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_negNSigmaRun3_negEta_tight")) {
    cut->AddCut(GetAnalysisCut("negTrack"));
    cut->AddCut(GetAnalysisCut("negEtaSel"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut("pid_TPCnSigma", "pid_TPCnSigma", kTRUE);
    cut_tpc_nSigma->AddCut(GetAnalysisCut("electronPID_TPCnsigma_tight"));

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut("pid_TOFnSigma", "pid_TOFnSigma", kTRUE);
    cut_tof_nSigma->AddCut(GetAnalysisCut("electronPID_TOFnsigma_tight"));

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut("e_NSigma", "e_NSigma", kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    cut->AddCut(cut_pid_OR);
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_GlobalTrackRun3_looseDCAxy")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightGlobalTrackRun3"));
    cut->AddCut(GetAnalysisCut("PrimaryTrack_looseDCAxy"));
    return cut;
  }

  if (!nameStr.compare("lmee_TPCTrackRun3_lowPt")) {
    cut->AddCut(GetAnalysisCut("lmeeLowBKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  if (!nameStr.compare("lmee_TPCTrackRun3")) {
    cut->AddCut(GetAnalysisCut("lmeeStandardKine"));
    cut->AddCut(GetAnalysisCut("TightTPCTrackRun3"));
    cut->AddCut(GetAnalysisCut("tightPrimaryTrack"));
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // lmee pair cuts

  if (!nameStr.compare("pairPhiV")) {
    AnalysisCompositeCut* cut_pairPhiV = new AnalysisCompositeCut("cut_pairPhiV", "cut_pairPhiV", kTRUE);
    cut_pairPhiV->AddCut(GetAnalysisCut("pairLowMass"));
    cut_pairPhiV->AddCut(GetAnalysisCut("pairPhiV"));
    cut->AddCut(cut_pairPhiV);
    return cut;
  }

  if (!nameStr.compare("excludePairPhiV")) {
    AnalysisCompositeCut* cut_pairlowPhiV = new AnalysisCompositeCut("cut_pairlowPhiV", "cut_pairlowPhiV", kFALSE);
    cut_pairlowPhiV->AddCut(GetAnalysisCut("excludePairLowMass"));
    cut_pairlowPhiV->AddCut(GetAnalysisCut("excludePairPhiV"));
    cut->AddCut(cut_pairlowPhiV);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts

  if (!nameStr.compare("muonQualityCutsMatchingOnly")) {
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonQualityCuts")) {
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("matchedQualityCuts")) {
    cut->AddCut(GetAnalysisCut("matchedQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt2")) {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt3")) {
    cut->AddCut(GetAnalysisCut("muonLowPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt4")) {
    cut->AddCut(GetAnalysisCut("muonLowPt4"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt5")) {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPt6")) {
    cut->AddCut(GetAnalysisCut("muonLowPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly2")) {
    cut->AddCut(GetAnalysisCut("muonLowPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly3")) {
    cut->AddCut(GetAnalysisCut("muonLowPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly4")) {
    cut->AddCut(GetAnalysisCut("muonLowPt4"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly5")) {
    cut->AddCut(GetAnalysisCut("muonLowPt5"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonLowPtMatchingOnly6")) {
    cut->AddCut(GetAnalysisCut("muonLowPt6"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonHighPt")) {
    cut->AddCut(GetAnalysisCut("muonHighPt"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonHighPt2")) {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonHighPt3")) {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCuts"));
    return cut;
  }

  if (!nameStr.compare("muonHighPtMatchingOnly2")) {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonHighPtMatchingOnly3")) {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
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

  // -----------------------------------------------------------
  // Pair cuts
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

  if (!nameStr.compare("pairMassLow4")) {
    cut->AddCut(GetAnalysisCut("pairMassLow4"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow5")) {
    cut->AddCut(GetAnalysisCut("pairMassLow5"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow6")) {
    cut->AddCut(GetAnalysisCut("pairMassLow6"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow7")) {
    cut->AddCut(GetAnalysisCut("pairMassLow7"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow8")) {
    cut->AddCut(GetAnalysisCut("pairMassLow8"));
    return cut;
  }

  if (!nameStr.compare("pairMassLow9")) {
    cut->AddCut(GetAnalysisCut("pairMassLow9"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz1")) {
    cut->AddCut(GetAnalysisCut("pairDalitz1"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz1Neg")) {
    cut->AddCut(GetAnalysisCut("pairDalitz1Neg"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz1Strong")) {
    cut->AddCut(GetAnalysisCut("pairDalitz1Strong"));
    return cut;
  }

  if (!nameStr.compare("pairDalitz1StrongNeg")) {
    cut->AddCut(GetAnalysisCut("pairDalitz1StrongNeg"));
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

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and in order or optimize the trigger
  // trigger selection cuts

  if (!nameStr.compare("jpsiO2TriggerTestCuts_LooseNsigma")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_LooseNsigma_corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug5"));
    return cut;
  }
  if (!nameStr.compare("jpsiO2TriggerTestCuts_MediumNsigma")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigmaOpen"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_MediumNsigma_corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug1"));
    return cut;
  }
  if (!nameStr.compare("jpsiO2TriggerTestCuts_TightNsigma")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("electronPIDnsigma"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_TightNsigma_corr")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_debug2"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_TPCPID1")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest1"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_TPCPID2")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest2"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_TPCPID3")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest3"));
    return cut;
  }

  if (!nameStr.compare("jpsiO2TriggerTestCuts_TPCPID4")) {
    cut->AddCut(GetAnalysisCut("jpsiStandardKine"));
    cut->AddCut(GetAnalysisCut("electronStandardQualityTriggerTest"));
    cut->AddCut(GetAnalysisCut("jpsi_TPCPID_TriggerTest4"));
    return cut;
  }

  if (!nameStr.compare("muonLooseTriggerTestCuts")) {
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLooseTriggerTestCuts_LowPt")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLooseTriggerTestCuts_HighPt2")) {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (!nameStr.compare("muonLooseTriggerTestCuts_HighPt3")) {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonLooseTriggerTestCuts"));
    return cut;
  }

  if (!nameStr.compare("muonHighPtMatchingOnly2")) {
    cut->AddCut(GetAnalysisCut("muonHighPt2"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonHighPtMatchingOnly3")) {
    cut->AddCut(GetAnalysisCut("muonHighPt3"));
    cut->AddCut(GetAnalysisCut("muonQualityCutsMatchingOnly"));
    return cut;
  }

  if (!nameStr.compare("muonMatchingMFTMCHTriggerTestCuts")) {
    cut->AddCut(GetAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  if (!nameStr.compare("muonMatchingMFTMCHTriggerTestCuts_LowPt")) {
    cut->AddCut(GetAnalysisCut("muonLowPt"));
    cut->AddCut(GetAnalysisCut("muonMatchingMFTMCHTriggerTestCuts"));
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

AnalysisCut* o2::aod::dqcuts::GetAnalysisCut(const char* cutName)
{
  //
  // define here cuts which are likely to be used often
  //
  AnalysisCut* cut = new AnalysisCut(cutName, cutName);
  std::string nameStr = cutName;
  // ---------------------------------------------------------------
  // Event cuts
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

  if (!nameStr.compare("eventTPCMultLow")) {
    cut->AddCut(VarManager::kMultTPC, 2, 50);
    return cut;
  }

  if (!nameStr.compare("int7vtxZ5")) {
    cut->AddCut(VarManager::kVtxZ, -5.0, 5.0);
    cut->AddCut(VarManager::kIsINT7, 0.5, 1.5);
    return cut;
  }

  // ---------------------------------------------------
  // Barrel track kine cuts
  if (!nameStr.compare("negTrack")) {
    cut->AddCut(VarManager::kCharge, -99., 0.);
    return cut;
  }

  if (!nameStr.compare("posTrack")) {
    cut->AddCut(VarManager::kCharge, 0., 99.);
    return cut;
  }

  if (!nameStr.compare("posEtaSel")) {
    cut->AddCut(VarManager::kEta, 0., 0.8);
    return cut;
  }

  if (!nameStr.compare("negEtaSel")) {
    cut->AddCut(VarManager::kEta, -0.8, 0.);
    return cut;
  }

  if (!nameStr.compare("etaSel")) {
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("lowMultTrackCut")) {
    cut->AddCut(VarManager::kPt, 0.15, 1000.0);
    cut->AddCut(VarManager::kTPCncls, 50.0, 159);
    return cut;
  }

  if (!nameStr.compare("jpsiStandardKine")) {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("jpsiStandardKine2")) {
    cut->AddCut(VarManager::kPt, 0.9, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("jpsiStandardKine3")) {
    cut->AddCut(VarManager::kPin, 1.0, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("lmeeStandardKine")) {
    cut->AddCut(VarManager::kPt, 0.2, 10.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("lmeeStandardKine_pt04")) {
    cut->AddCut(VarManager::kPt, 0.4, 10.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("lmeeLowBKine")) {
    cut->AddCut(VarManager::kPt, 0.075, 10.0);
    cut->AddCut(VarManager::kEta, -0.8, 0.8);
    return cut;
  }

  if (!nameStr.compare("dalitzStandardKine")) {
    cut->AddCut(VarManager::kPt, 0.15, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  if (!nameStr.compare("PIDStandardKine")) {
    cut->AddCut(VarManager::kPt, 0.1, 1000.0);
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    return cut;
  }

  // -----------------------------------------------
  // Barrel track quality cuts

  // Run 2 only

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

  if (!nameStr.compare("TightGlobalTrack")) {
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

  // Run 2 or run 3

  if (!nameStr.compare("jpsi_trackCut_debug")) {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    cut->AddCut(VarManager::kITSncls, 2.5, 7.5);
    return cut;
  }

  if (!nameStr.compare("jpsi_trackCut_noITSCuts_debug")) {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 90., 159);
    return cut;
  }

  if (!nameStr.compare("lmee_trackCut_debug")) {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 80., 159);
    cut->AddCut(VarManager::kITSncls, 0.0, 7.5);
    return cut;
  }

  if (!nameStr.compare("lmee_skimming_cuts")) {
    cut->AddCut(VarManager::kEta, -0.9, 0.9);
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 60.0, 161.);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
    return cut;
  }

  if (!nameStr.compare("TightGlobalTrackRun3")) {
    cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 4.5, 7.5);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kTPCncls, 90.0, 170.);
    return cut;
  }

  if (!nameStr.compare("TightTPCTrackRun3")) {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    return cut;
  }

  if (!nameStr.compare("TightTPCTrack")) {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCnclsCR, 80.0, 161.);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
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

  if (!nameStr.compare("electronStandardQualityForO2MCdebug3")) {
    cut->AddCut(VarManager::kITSncls, 0.5, 10);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (!nameStr.compare("electronStandardQualityITSOnly")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kITSchi2, 0.0, 5.0);
    cut->AddCut(VarManager::kITSncls, 3.5, 7.5);
    return cut;
  }

  if (!nameStr.compare("electronStandardQualityTPCOnly")) {
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kTPCncls, 70, 161.);
    return cut;
  }

  if (!nameStr.compare("pidbasic")) {
    cut->AddCut(VarManager::kTPCnclsCR, 70, 161);
    cut->AddCut(VarManager::kTPCchi2, 0, 4);
    return cut;
  }

  if (!nameStr.compare("standardPrimaryTrack")) {
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);
    cut->AddCut(VarManager::kTrackDCAz, -3.0, 3.0);
    return cut;
  }
  if (!nameStr.compare("standardPrimaryTrackDCA")) {
    cut->AddCut(VarManager::kTrackDCAxy, -0.1, 0.1);
    cut->AddCut(VarManager::kTrackDCAz, -0.15, 0.15);
    return cut;
  }

  if (!nameStr.compare("PrimaryTrack_looseDCAxy")) {
    cut->AddCut(VarManager::kTrackDCAxy, -1.0, 1.0);   // in cm
    cut->AddCut(VarManager::kTrackDCAsigZ, -3.0, 3.0); // in sigma
    return cut;
  }

  if (!nameStr.compare("tightPrimaryTrack")) {
    cut->AddCut(VarManager::kTrackDCAsigXY, -3.0, 3.0);
    cut->AddCut(VarManager::kTrackDCAsigZ, -3.0, 3.0);
    return cut;
  }

  // -----------------------------------------------------
  // V0 and Dalitz legs selections

  for (int i = 1; i <= 8; i++) {
    if (!nameStr.compare(Form("dalitzLeg%d", i))) {
      cut->AddCut(VarManager::kIsDalitzLeg + i - 1, 0.5, 1.5);
      return cut;
    }

    if (!nameStr.compare(Form("notDalitzLeg%d", i))) {
      cut->AddCut(VarManager::kIsDalitzLeg + i - 1, -0.5, 0.5);
      return cut;
    }
  }

  if (!nameStr.compare("pidcalib_ele")) {
    cut->AddCut(VarManager::kIsLegFromGamma, 0.5, 1.5, false);
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

  // ------------------------------------------------
  // Barrel PID cuts
  if (!nameStr.compare("jpsi_TPCPID_debug1")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_debug2")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_debug3")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.5, 999);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_debug4")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPi, 3.0, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr, 3.0, 999);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_debug5")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -4.0, 4.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999);
    return cut;
  }

  if (!nameStr.compare("lmee_TPCPID_debug1")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -5.0, 5.0);
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

  if (!nameStr.compare("electronPID_TPCnsigma_corr")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("electronPID_TPCnsigma_tight")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPr, -3., 4., true, VarManager::kPin, 0.0, 1e+10, false);
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
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.4, 1e+10, false);
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

  if (!nameStr.compare("electronPIDPrKaPiRej_Corr")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (!nameStr.compare("electronPIDPrKaPiRejLoose")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (!nameStr.compare("electronPIDPrKaPiRejLoose_Corr")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, -2.0, 2.0, true);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 3.0, true, VarManager::kPin, 0.0, 1.0, false);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, -3.0, 2.0, true, VarManager::kPin, 0.0, 1.0, true);
    cut->AddCut(VarManager::kTPCnSigmaKa, -3.0, 3.0, true);
    return cut;
  }

  if (!nameStr.compare("electronPIDOnly")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3.0, 3.0);
    return cut;
  }

  if (!nameStr.compare("electronPIDOnly_Corr")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 3.0);
    return cut;
  }

  if (!nameStr.compare("kaonPIDnsigma2")) {
    cut->AddCut(VarManager::kTPCnSigmaKa, -2.0, 2.0);
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

  if (!nameStr.compare("electronPID_TOFnsigma")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3., 3., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.4, 1e+10, false);
    return cut;
  }

  if (!nameStr.compare("electronPID_TOFnsigma_tight")) {
    cut->AddCut(VarManager::kTPCnSigmaEl, -3., 2., false, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTPCnSigmaPi, -3., 3., true, VarManager::kPin, 0.0, 1e+10, false);
    cut->AddCut(VarManager::kTOFnSigmaEl, -3., 3., false, VarManager::kPin, 0.3, 1e+10, false);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  // Muon cuts
  if (!nameStr.compare("muonQualityCuts")) {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (!nameStr.compare("matchedQualityCuts")) {
    cut->AddCut(VarManager::kEta, -4.0, -2.5);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
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
    cut->AddCut(VarManager::kPt, 0.7, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonLowPt3")) {
    cut->AddCut(VarManager::kPt, 0.8, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonLowPt4")) {
    cut->AddCut(VarManager::kPt, 0.9, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonLowPt5")) {
    cut->AddCut(VarManager::kPt, 1.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonLowPt6")) {
    cut->AddCut(VarManager::kPt, 2.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonHighPt")) {
    cut->AddCut(VarManager::kPt, 3.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonHighPt2")) {
    cut->AddCut(VarManager::kPt, 4.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("muonHighPt3")) {
    cut->AddCut(VarManager::kPt, 6.0, 1000.0);
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

  // -----------------------------------------------------------------------------------------------
  // Pair cuts
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

  if (!nameStr.compare("pairDalitz1Neg")) { // To be used when negative field
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    TF1* fcutHigh = new TF1("f1", "[0] + [0]/[1]*x", -1.5, 1.5);
    fcutHigh->SetParameters(0.6, 0.12);
    TF1* fcutLow = new TF1("f2", "-[0] - [0]/[1]*x", -1.5, 1.5);
    fcutLow->SetParameters(0.6, 0.12);
    cut->AddCut(VarManager::kPsiPair, fcutLow, fcutHigh, true, VarManager::kDeltaPhiPair, -0.12, 0.);
    return cut;
  }

  if (!nameStr.compare("pairDalitz1Strong")) {
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    cut->AddCut(VarManager::kDeltaPhiPair, -1., 0.);
    return cut;
  }

  if (!nameStr.compare("pairDalitz1StrongNeg")) { // To be used when negative field
    cut->AddCut(VarManager::kMass, 0.0, 0.015, false, VarManager::kPt, 0., 1.);
    cut->AddCut(VarManager::kMass, 0.0, 0.035, false, VarManager::kPt, 0., 1., true);
    cut->AddCut(VarManager::kDeltaPhiPair, 0., 1.);
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
    cut->AddCut(VarManager::kMass, 1.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow2")) {
    cut->AddCut(VarManager::kMass, 1.5, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow3")) {
    cut->AddCut(VarManager::kMass, 1.8, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow5")) {
    cut->AddCut(VarManager::kMass, 1.85, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow6")) {
    cut->AddCut(VarManager::kMass, 1.9, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow7")) {
    cut->AddCut(VarManager::kMass, 2.0, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow8")) {
    cut->AddCut(VarManager::kMass, 2.2, 1000.0);
    return cut;
  }

  if (!nameStr.compare("pairMassLow9")) {
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

  if (!nameStr.compare("pairPhiV")) {
    cut->AddCut(VarManager::kPairPhiv, 2., 3.2);
    return cut;
  }

  if (!nameStr.compare("excludePairPhiV")) {
    cut->AddCut(VarManager::kPairPhiv, 2., 3.2, true);
    return cut;
  }

  if (!nameStr.compare("pairLowMass")) {
    cut->AddCut(VarManager::kMass, 0., 0.1);
    return cut;
  }

  if (!nameStr.compare("excludePairLowMass")) {
    cut->AddCut(VarManager::kMass, 0., 0.1, true);
    return cut;
  }

  // -------------------------------------------------------------------------------------------------
  //
  // Below are a list of single electron single muon and pair selection in order or optimize the trigger
  // trigger selection cuts

  if (!nameStr.compare("electronStandardQualityTriggerTest")) {
    cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    cut->AddCut(VarManager::kIsITSrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kIsTPCrefit, 0.5, 1.5);
    cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    cut->AddCut(VarManager::kITSchi2, 0.1, 36.0);
    cut->AddCut(VarManager::kTPCncls, 50.0, 161.);
    return cut;
  }

  if (!nameStr.compare("muonLooseTriggerTestCuts")) {
    cut->AddCut(VarManager::kEta, -4.5, -2.0);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    return cut;
  }

  if (!nameStr.compare("muonMatchingMFTMCHTriggerTestCuts")) {
    cut->AddCut(VarManager::kEta, -4.5, -2.0);
    cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 10, 100);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 1500, false, VarManager::kMuonRAtAbsorberEnd, 10, 26.5);
    cut->AddCut(VarManager::kMuonPDca, 0.0, 800, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 100);
    cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_TriggerTest1")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_TriggerTest2")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 2.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_TriggerTest3")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 3.0, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  if (!nameStr.compare("jpsi_TPCPID_TriggerTest4")) {
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -3.0, 4.0, false, VarManager::kPin, 0, 3.0);
    cut->AddCut(VarManager::kTPCnSigmaEl_Corr, -2.0, 4.0, false, VarManager::kPin, 3.0, 9999.0);
    cut->AddCut(VarManager::kTPCnSigmaPi_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    cut->AddCut(VarManager::kTPCnSigmaPr_Corr, 2.5, 999, false, VarManager::kPin, 0, 2.0);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}
