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
#include <vector>
#include <algorithm>
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "VarManager.h"
#include "CommonConstants/MathConstants.h"

void o2::aod::dqhistograms::DefineHistograms(HistogramManager* hm, const char* histClass, const char* groupName, const char* subGroupName)
{
  //
  // Add a predefined group of histograms to the HistogramManager hm and histogram class histClass
  // NOTE: The subGroupName argument may contain several keywords, but the user should take care of ambiguities. TODO: fix it!
  // NOTE: All of the histograms which match any of the group or subgroup names will be added to the same histogram class !!
  //            So one has to make sure not to mix e.g. event-wise with track-wise histograms
  // NOTE: The subgroup name can be empty. In this case just a minimal set of histograms corresponding to the group name will be defined
  //
  TString groupStr = groupName;
  groupStr.ToLower();
  TString subGroupStr = subGroupName;
  subGroupStr.ToLower();
  if (!groupStr.CompareTo("event")) {
    if (!subGroupStr.Contains("generator")) {
      hm->AddHistogram(histClass, "VtxZ", "Vtx Z", false, 60, -15.0, 15.0, VarManager::kVtxZ);
      hm->AddHistogram(histClass, "VtxZ_Run", "Vtx Z", true, 1, -0.5, 0.5, VarManager::kRunNo, 60, -15.0, 15.0, VarManager::kVtxZ, 1, 0, 1, VarManager::kNothing, "", "", "", VarManager::kNothing, VarManager::kNothing, false, true);
      hm->AddHistogram(histClass, "BC", "Event per BC", false, 3564, 0.0, 3564.0, VarManager::kBCOrbit);
    }
    if (subGroupStr.Contains("trigger")) {
      hm->AddHistogram(histClass, "IsINT7", "Is INT7", false, 2, -0.5, 1.5, VarManager::kIsINT7);
      if (subGroupStr.Contains("muon") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsINT7inMUON", "INT7inMUON", false, 2, -0.5, 1.5, VarManager::kIsINT7inMUON);
        hm->AddHistogram(histClass, "IsMuonSingleLowPt7", "Is MuonSingleLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonSingleLowPt7);
        hm->AddHistogram(histClass, "IsMuonSingleHighPt7", "Is MuonSingleHighPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonSingleHighPt7);
        hm->AddHistogram(histClass, "IsMuonUnlikeLowPt7", "Is MuonUnlikeLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonUnlikeLowPt7);
        hm->AddHistogram(histClass, "IsMuonLikeLowPt7", "Is MuonLikeLowPt7", false, 2, -0.5, 1.5, VarManager::kIsMuonLikeLowPt7);
      }
      if (subGroupStr.Contains("up") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsCUP8", "CUP8", false, 2, -0.5, 1.5, VarManager::kIsCUP8);
        hm->AddHistogram(histClass, "IsCUP9", "CUP9", false, 2, -0.5, 1.5, VarManager::kIsCUP9);
        hm->AddHistogram(histClass, "IsMUP10", "MUP10", false, 2, -0.5, 1.5, VarManager::kIsMUP10);
        hm->AddHistogram(histClass, "IsMUP11", "MUP11", false, 2, -0.5, 1.5, VarManager::kIsMUP11);
      }
      if (subGroupStr.Contains("emc") || subGroupStr.Contains("all")) {
        hm->AddHistogram(histClass, "IsEMC7", "EMC7", false, 2, -0.5, 1.5, VarManager::kIsEMC7);
      }
    }
    if (subGroupStr.Contains("time")) {
      hm->AddHistogram(histClass, "CollTime", "Coll. time wrt BC", false, 100, 0.0, 100.0, VarManager::kCollisionTime);
      hm->AddHistogram(histClass, "CollTimeRes", "Coll. time resolution", false, 100, 0.0, 200.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTime_VtxZ", "Coll. time wrt BC vs vtx-z", false, 50, -15.0, 15., VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kCollisionTime);
      hm->AddHistogram(histClass, "CollTimeRes_VtxZ", "Coll. time resolution ", false, 50, -15.0, 15., VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTimeRes_MultTPC", "Coll. time resolution ", false, 50, 0.0, 50000., VarManager::kMultTPC, 100, 0.0, 200.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "CollTimeRes_MultPV", "Coll. time resolution ", false, 100, 0.0, 4000., VarManager::kVtxNcontribReal, 100, 0.0, 200.0, VarManager::kCollisionTimeRes);
      hm->AddHistogram(histClass, "TimeFromSOR", "Time since SOR", false, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR);
    }
    if (subGroupStr.Contains("vtx")) {
      hm->AddHistogram(histClass, "VtxX", "Vtx X", false, 200, -0.1, 0.1, VarManager::kVtxX);
      hm->AddHistogram(histClass, "VtxY", "Vtx Y", false, 200, -0.1, 0.1, VarManager::kVtxY);
      hm->AddHistogram(histClass, "VtxYVtxX", "Vtx Y vs Vtx X", false, 200, -0.06, 0.0, VarManager::kVtxX, 200, -0.03, 0.03, VarManager::kVtxY);
    }
    if (subGroupStr.Contains("vtxpp")) {
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("vtxPbPb")) {
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 20000.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("cent")) {
      hm->AddHistogram(histClass, "CentFT0C", "CentFT0C", false, 100, 0., 100., VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "CentFT0C_vtxZ", "CentFT0C vs Vtx Z", false, 60, -15.0, 15.0, VarManager::kVtxZ, 20, 0., 100., VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "CentFT0C_MultTPC", "CentFT0C vs MultTPC", false, 100, 0., 100., VarManager::kCentFT0C, 100, 0., 50000., VarManager::kMultTPC);
      hm->AddHistogram(histClass, "CentFT0C_Run", "Cent FT0C", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, 0., 100., VarManager::kCentFT0C, 1, 0, 1, VarManager::kNothing, "", "", "", VarManager::kNothing, VarManager::kNothing, false, true);
    }
    if (subGroupStr.Contains("mult")) {
      if (subGroupStr.Contains("pp")) {
        hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 250, 0.0, 500.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "MultFV0A", "MultFV0A", false, 250, 0.0, 500.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 300, 0.0, 300.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 300, 0.0, 300.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFDDA", "MultFDDA", false, 300, 0.0, 300.0, VarManager::kMultFDDA);
        hm->AddHistogram(histClass, "MultFDDC", "MultFDDC", false, 50, 0.0, 50.0, VarManager::kMultFDDC);
        hm->AddHistogram(histClass, "MultTracklets", "MultTracklets", false, 250, 0.0, 250.0, VarManager::kMultTracklets);
        hm->AddHistogram(histClass, "VtxNContribReal", "Vtx n contributors", false, 200, 0.0, 200.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MultNTracksPVeta1", "MultNTracksPVeta1", false, 200, 0, 200.0, VarManager::kMultNTracksPVeta1);
        hm->AddHistogram(histClass, "MultNTracksPVetaHalf", "MultNTracksPVetaHalf", false, 200, 0, 200.0, VarManager::kMultNTracksPVetaHalf);
        hm->AddHistogram(histClass, "MultTPC_MultFV0A", "MultTPC vs MultFV0A", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 500.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0A", "MultTPC vs MultFT0A", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 200.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0C", "MultTPC vs MultFT0C", false, 100, 0, 500.0, VarManager::kMultTPC, 100, 0, 300.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFT0A_MultFT0C", "MultFT0A vs MultFT0C", false, 100, 0, 200.0, VarManager::kMultFT0A, 100, 0, 300.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultITSWithPV", "MultITSWithPV", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasITS);
        hm->AddHistogram(histClass, "MultTPCWithPV", "MultTPCWithPV", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasTPC);
        hm->AddHistogram(histClass, "MultITSTPCWithPV", "MultITSTPCWithPV", false, 150, 0.0, 150.0, VarManager::kMultNTracksITSTPC);
        hm->AddHistogram(histClass, "MultITSOnly", "MultITSOnly", false, 150, 0.0, 150.0, VarManager::kMultNTracksITSOnly);
        hm->AddHistogram(histClass, "MultITSWithPV_MultTPCWithPV", "MultITSWithPV_MultTPCWithPV", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasITS, 150, 0.0, 150.0, VarManager::kMultNTracksHasTPC);
        hm->AddHistogram(histClass, "MultITSWithPV_MultITSTPCWithPV", "MultITSWithPV_MultTPCWithPV", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasITS, 150, 0.0, 150.0, VarManager::kMultNTracksITSTPC);
        hm->AddHistogram(histClass, "MultITSWithPV_MultFT0C", "MultITSWithPV_MultFT0C", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasITS, 250, 0.0, 2500.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultITSWithPV_MultFT0A", "MultITSWithPV_MultFT0A", false, 150, 0.0, 150.0, VarManager::kMultNTracksHasITS, 250, 0.0, 2500.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultITSTPCWithPV_MultFT0C", "MultITSTPCWithPV_MultFT0C", false, 150, 0.0, 150.0, VarManager::kMultNTracksITSTPC, 250, 0.0, 2500.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultITSTPCWithPV_MultFT0A", "MultITSTPCWithPV_MultFT0A", false, 150, 0.0, 150.0, VarManager::kMultNTracksITSTPC, 250, 0.0, 2500.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "VtxZ_MultITSWithPV", "VtxZ vs MultITSWithPV", false, 240, -12.0, 12.0, VarManager::kVtxZ, 400, 0, 400.0, VarManager::kMultNTracksHasITS);
        hm->AddHistogram(histClass, "VtxZ_MultTPCWithPV", "VtxZ vs MultTPCWithPV", false, 240, -12.0, 12.0, VarManager::kVtxZ, 400, 0, 400.0, VarManager::kMultNTracksHasTPC);
        hm->AddHistogram(histClass, "VtxZ_MultITSTPCWithPV", "VtxZ vs MultITSTPCWithPV", false, 240, -12.0, 12.0, VarManager::kVtxZ, 400, 0, 400.0, VarManager::kMultNTracksITSTPC);
        hm->AddHistogram(histClass, "VtxZ_MultITSOnly", "VtxZ vs MultITSOnly", false, 240, -12.0, 12.0, VarManager::kVtxZ, 400, 0, 400.0, VarManager::kMultNTracksITSOnly);
        hm->AddHistogram(histClass, "VtxZ_VtxNcontribReal", "VtxZ vs VtxNcontribReal", false, 240, -12.0, 12.0, VarManager::kVtxZ, 200, 0, 200.0, VarManager::kVtxNcontribReal);

      } else {
        hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 200, 0.0, 50000.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "MultTPC_vsTimeSOR", "MultTPC vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 50000.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "MultFV0A", "MultFV0A", false, 200, 0.0, 300000.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 200, 0.0, 300000.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 200, 0.0, 100000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFDDA", "MultFDDA", false, 100, 0.0, 100000.0, VarManager::kMultFDDA);
        hm->AddHistogram(histClass, "MultFDDC", "MultFDDC", false, 100, 0.0, 100000.0, VarManager::kMultFDDC);
        hm->AddHistogram(histClass, "MultZNA", "MultZNA", false, 400, 0.0, 400.0, VarManager::kMultZNA);
        hm->AddHistogram(histClass, "MultZNC", "MultZNC", false, 400, 0.0, 400.0, VarManager::kMultZNC);
        hm->AddHistogram(histClass, "MultZNA_ZNC", "MultZNA vs ZNC", false, 400, 0.0, 400.0, VarManager::kMultZNA, 400, 0.0, 400.0, VarManager::kMultZNC);
        hm->AddHistogram(histClass, "MultTracklets", "MultTracklets", false, 100, 0.0, 25000.0, VarManager::kMultTracklets);
        hm->AddHistogram(histClass, "VtxNContribReal", "Vtx n contributors (real)", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "VtxNContribReal_vsTimeSOR", "VtxNContribReal vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 5000.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 20000.0, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MultTPC_MultFV0A", "MultTPC vs MultFV0A", false, 100, 0, 50000.0, VarManager::kMultTPC, 100, 0, 300000.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0A", "MultTPC vs MultFT0A", false, 100, 0, 50000.0, VarManager::kMultTPC, 100, 0, 300000.0, VarManager::kMultFT0A);
        hm->AddHistogram(histClass, "MultTPC_MultFT0C", "MultTPC vs MultFT0C", false, 100, 0, 50000.0, VarManager::kMultTPC, 100, 0, 100000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultFT0A_MultFT0C", "MultFT0A vs MultFT0C", false, 100, 0, 100000.0, VarManager::kMultFT0A, 100, 0, 300000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "VtxNContribReal_MultTPC", "Vtx n contributors (real) vs mult TPC", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 200, 0.0, 50000.0, VarManager::kMultTPC);
        hm->AddHistogram(histClass, "VtxNContribReal_ZNA", "Vtx n contributors (real) vs ZNA", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 200, 0.0, 400.0, VarManager::kMultZNA);
        hm->AddHistogram(histClass, "VtxNContribReal_ZNC", "Vtx n contributors (real) vs ZNC", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 200, 0.0, 400.0, VarManager::kMultZNC);
        hm->AddHistogram(histClass, "MultZNA_FT0C", "MultZNA vs FT0C", false, 400, 0.0, 400.0, VarManager::kMultZNA, 200, 0.0, 100000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "MultZNC_FT0C", "MultZNC vs FT0C", false, 400, 0.0, 400.0, VarManager::kMultZNC, 200, 0.0, 100000.0, VarManager::kMultFT0C);
        hm->AddHistogram(histClass, "TPCpileupZA", "TPC pileup Z, A-side", false, 200, -50.0, 50.0, VarManager::kNTPCpileupZA);
        hm->AddHistogram(histClass, "TPCpileupZC", "TPC pileup Z, C-side", false, 200, -50.0, 50.0, VarManager::kNTPCpileupZC);
        hm->AddHistogram(histClass, "TPCpileupNcontribA", "TPC pileup n-contributors, A-side", false, 300, 0.0, 3000.0, VarManager::kNTPCpileupContribA);
        hm->AddHistogram(histClass, "TPCpileupNcontribC", "TPC pileup n-contributors, C-side", false, 300, 0.0, 3000.0, VarManager::kNTPCpileupContribC);
        hm->AddHistogram(histClass, "TPCoccupContribLongA", "TPC occupancy from pileup, n-contrib, A-side, long time range", false, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongA);
        hm->AddHistogram(histClass, "TPCoccupContribLongAvsTime", "TPC occupancy from pileup, n-contrib, A-side, long time range", true, 1000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 10000.0, VarManager::kNTPCcontribLongA);
        hm->AddHistogram(histClass, "TPCoccupContribLongAvsContribPV", "TPC occupancy from pileup, n-contrib, A-side, long time range, vs n.contrib", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongA);
        hm->AddHistogram(histClass, "TPCoccupContribLongC", "TPC occupancy from pileup, n-contrib, C-side, long time range", false, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongC);
        hm->AddHistogram(histClass, "TPCoccupContribLongCvsTime", "TPC occupancy from pileup, n-contrib, C-side, long time range", true, 1000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 10000.0, VarManager::kNTPCcontribLongC);
        hm->AddHistogram(histClass, "TPCoccupContribLongCvsContribPV", "TPC occupancy from pileup, n-contrib, C-side, long time range, vs n.contrib", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongC);
        hm->AddHistogram(histClass, "TPCoccupContribLongAvsC", "TPC occupancy from pileup, n-contrib, A-side vs C-side, long time range", false, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongA, 100, 0.0, 10000.0, VarManager::kNTPCcontribLongC);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeLongA", "TPC occupancy from pileup, mean time, A-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmeanTimeLongA);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeLongC", "TPC occupancy from pileup, mean time, C-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmeanTimeLongC);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeLongAvsC", "TPC occupancy from pileup, mean time, A-side vs C-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmeanTimeLongA, 100, -100.0, 100.0, VarManager::kNTPCmeanTimeLongC);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeLongA", "TPC occupancy from pileup, median time, A-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmedianTimeLongA);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeLongC", "TPC occupancy from pileup, median time, C-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmedianTimeLongC);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeLongAvsC", "TPC occupancy from pileup, median time, A-side vs C-side, long time range", false, 100, -100.0, 100.0, VarManager::kNTPCmedianTimeLongA, 100, -100.0, 100.0, VarManager::kNTPCmedianTimeLongC);
        hm->AddHistogram(histClass, "TPCoccupContribShortA", "TPC occupancy from pileup, n-contrib, A-side, short time range", false, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortA);
        hm->AddHistogram(histClass, "TPCoccupContribShortAvsTime", "TPC occupancy from pileup, n-contrib, A-side, short time range", true, 1000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 10000.0, VarManager::kNTPCcontribShortA);
        hm->AddHistogram(histClass, "TPCoccupContribShortAvsContribPV", "TPC occupancy from pileup, n-contrib, A-side, short time range, vs n.contrib", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortA);
        hm->AddHistogram(histClass, "TPCoccupContribShortC", "TPC occupancy from pileup, n-contrib, C-side, short time range", false, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortC);
        hm->AddHistogram(histClass, "TPCoccupContribShortCvsTime", "TPC occupancy from pileup, n-contrib, C-side, short time range", true, 1000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, 0.0, 10000.0, VarManager::kNTPCcontribShortC);
        hm->AddHistogram(histClass, "TPCoccupContribShortAvsC", "TPC occupancy from pileup, n-contrib, A-side vs C-side, short time range", false, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortA, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortC);
        hm->AddHistogram(histClass, "TPCoccupContribShortCvsContribPV", "TPC occupancy from pileup, n-contrib, C-side, short time range, vs n.contrib", false, 100, 0.0, 5000.0, VarManager::kVtxNcontribReal, 100, 0.0, 7000.0, VarManager::kNTPCcontribShortC);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeShortA", "TPC occupancy from pileup, mean time, A-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmeanTimeShortA);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeShortC", "TPC occupancy from pileup, mean time, C-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmeanTimeShortC);
        hm->AddHistogram(histClass, "TPCoccupMeanTimeShortAvsC", "TPC occupancy from pileup, mean time, A-side vs C-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmeanTimeShortA, 100, -20.0, 20.0, VarManager::kNTPCmeanTimeShortC);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeShortA", "TPC occupancy from pileup, median time, A-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmedianTimeShortA);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeShortC", "TPC occupancy from pileup, median time, C-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmedianTimeShortC);
        hm->AddHistogram(histClass, "TPCoccupMedianTimeShortAvsC", "TPC occupancy from pileup, median time, A-side vs C-side, short time range", false, 100, -20.0, 20.0, VarManager::kNTPCmedianTimeShortA, 100, -20.0, 20.0, VarManager::kNTPCmedianTimeShortC);
        hm->AddHistogram(histClass, "NcontribReal_centT0C", "Ncontrib vs Cent", false, 100, 0, 100, VarManager::kCentFT0C, 4000, 0, 4000, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "globalTracks_centT0C", "globalTracks vs Cent", false, 100, 0, 100, VarManager::kCentFT0C, 4000, 0, 4000, VarManager::kMultA);
        hm->AddHistogram(histClass, "ITSTPCTracks_centT0C", "ITSTPCTracks vs Cent", false, 100, 0, 100, VarManager::kCentFT0C, 4000, 0, 4000, VarManager::kMultAllTracksITSTPC);
      }
    }
    if (subGroupStr.Contains("ftmulpbpb")) {
      hm->AddHistogram(histClass, "MultTPC", "MultTPC", false, 100, 0.0, 50000.0, VarManager::kMultTPC);
      hm->AddHistogram(histClass, "MultFT0C", "MultFT0C", false, 100, 0.0, 60000.0, VarManager::kMultFT0C);
      hm->AddHistogram(histClass, "MultFT0A", "MultFT0A", false, 100, 0.0, 180000.0, VarManager::kMultFT0A);
      hm->AddHistogram(histClass, "VtxNContribReal", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontribReal);
      hm->AddHistogram(histClass, "VtxNContrib", "Vtx n contributors", false, 100, 0.0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultTPC_VtxNContrib", "MultTPC vs VtxNContrib", false, 100, 0, 50000.0, VarManager::kMultTPC, 100, 0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultFT0C_VtxNContrib", "MultFT0C vs VtxNContrib", false, 100, 0, 60000.0, VarManager::kMultFT0C, 100, 0, 10000.0, VarManager::kVtxNcontrib);
      hm->AddHistogram(histClass, "MultFT0A_VtxNContrib", "MultFT0A vs VtxNContrib", false, 100, 0, 180000.0, VarManager::kMultFT0A, 100, 0, 10000.0, VarManager::kVtxNcontrib);
    }
    if (subGroupStr.Contains("occupancy")) {
      hm->AddHistogram(histClass, "ITStrackOccupancy", "ITStrackOccupancy", false, 200, 0.0, 20000.0, VarManager::kTrackOccupancyInTimeRange);
      hm->AddHistogram(histClass, "Ft0cOccupancy", "Ft0cOccupancy", false, 200, 0.0, 20000.0, VarManager::kFT0COccupancyInTimeRange);
    }
    if (subGroupStr.Contains("mc")) {
      hm->AddHistogram(histClass, "MCVtxX_VtxX", "Vtx X (MC vs rec)", false, 100, -0.5, 0.5, VarManager::kVtxX, 100, -0.5, 0.5, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxY_VtxY", "Vtx Y (MC vs rec)", false, 100, -0.5, 0.5, VarManager::kVtxY, 100, -0.5, 0.5, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxZ_VtxZ", "Vtx Z (MC vs rec)", false, 75, -15.0, 15.0, VarManager::kVtxZ, 75, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCVtxZ", "Vtx Z (MC)", false, 75, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCImpPar_CentVZERO", "MC impact param vs CentVZERO", false, 50, 0.0, 100.0, VarManager::kCentVZERO, 20, 0.0, 20.0, VarManager::kMCEventImpParam);
    }
    if (subGroupStr.Contains("generator")) {
      hm->AddHistogram(histClass, "MCVtxX", "Vtx X", false, 1000, -0.5, 0.5, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxY", "Vtx Y", false, 1000, -0.5, 0.5, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxX_VtxY", "Vtx X vs Vtx Y", false, 200, -0.2, 0.2, VarManager::kMCVtxX, 200, -0.2, 0.2, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCVtxZ", "Vtx Z", false, 60, -15.0, 15.0, VarManager::kMCVtxZ);
      hm->AddHistogram(histClass, "MCVtxZ_VtxX", "Vtx X vs Vtx Z", false, 60, -15.0, 15.0, VarManager::kMCVtxZ, 200, -0.2, 0.2, VarManager::kMCVtxX);
      hm->AddHistogram(histClass, "MCVtxX_VtxY", "Vtx X vs Vtx Y", false, 200, 15.0, 15.0, VarManager::kMCVtxZ, 200, -0.2, 0.2, VarManager::kMCVtxY);
      hm->AddHistogram(histClass, "MCImpPar", "MC impact param", false, 20, 0.0, 20.0, VarManager::kMCEventImpParam);
    }
    if (subGroupStr.Contains("subgen")) {
      hm->AddHistogram(histClass, "SubGenID", "SubGenerator ID", false, 11, -0.5, 10.5, VarManager::kMCEventSubGeneratorId);
    }
    if (subGroupStr.Contains("qvector")) {
      int varZNA[3] = {VarManager::kQ1ZNAX, VarManager::kQ1ZNAY, VarManager::kCentFT0C};
      int varZNC[3] = {VarManager::kQ1ZNCX, VarManager::kQ1ZNCY, VarManager::kCentFT0C};

      int bins[3] = {500, 500, 18};
      double minBins[3] = {-10, -10, 0};
      double maxBins[3] = {10, 10, 90};
      hm->AddHistogram(histClass, "Q1ZNAX_Q1ZNAY_CentFT0C", "", 3, varZNA, bins, minBins, maxBins, 0, -1, kTRUE);
      hm->AddHistogram(histClass, "Q1ZNCX_Q1ZNCY_CentFT0C", "", 3, varZNC, bins, minBins, maxBins, 0, -1, kTRUE);

      hm->AddHistogram(histClass, "IntercalibZNA_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -50.0, 50.0, VarManager::KIntercalibZNA);
      hm->AddHistogram(histClass, "IntercalibZNC_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -50.0, 50.0, VarManager::KIntercalibZNC);

      hm->AddHistogram(histClass, "EnergyCommonZNA", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyCommonZNA);
      hm->AddHistogram(histClass, "EnergyCommonZNC", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyCommonZNC);

      hm->AddHistogram(histClass, "EnergyZNA1", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNA1);
      hm->AddHistogram(histClass, "EnergyZNA2", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNA2);
      hm->AddHistogram(histClass, "EnergyZNA3", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNA3);
      hm->AddHistogram(histClass, "EnergyZNA4", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNA4);

      hm->AddHistogram(histClass, "EnergyZNC1", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNC1);
      hm->AddHistogram(histClass, "EnergyZNC2", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNC2);
      hm->AddHistogram(histClass, "EnergyZNC3", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNC3);
      hm->AddHistogram(histClass, "EnergyZNC4", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 2000, 0, 2000, VarManager::kEnergyZNC4);

      hm->AddHistogram(histClass, "Q2X0A", "", false, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q2X0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q2X0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q3X0A", "", false, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B", "", false, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C", "", false, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q3X0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C_VtxZ", "", true, 60, -15.0, 15.0, VarManager::kVtxZ, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q3X0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Q3X0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0B);
      hm->AddHistogram(histClass, "Q3Y0B_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0B);
      hm->AddHistogram(histClass, "Q3X0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3X0C);
      hm->AddHistogram(histClass, "Q3Y0C_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kQ3Y0C);
      hm->AddHistogram(histClass, "Q2X0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0A);
      hm->AddHistogram(histClass, "Q2Y0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0A);
      hm->AddHistogram(histClass, "Q2X0B_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0B);
      hm->AddHistogram(histClass, "Q2Y0B_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0B);
      hm->AddHistogram(histClass, "Q2X0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2X0C);
      hm->AddHistogram(histClass, "Q2Y0C_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ2Y0C);
      hm->AddHistogram(histClass, "Q3X0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ3X0A);
      hm->AddHistogram(histClass, "Q3Y0A_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kQ3Y0A);
      hm->AddHistogram(histClass, "Psi2A", "", false, 100, -2.0, 2.0, VarManager::kPsi2A);
      hm->AddHistogram(histClass, "Psi2B", "", false, 100, -2.0, 2.0, VarManager::kPsi2B);
      hm->AddHistogram(histClass, "Psi2C", "", false, 100, -2.0, 2.0, VarManager::kPsi2C);
      hm->AddHistogram(histClass, "Psi2A_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2A);
      hm->AddHistogram(histClass, "Psi2B_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2B);
      hm->AddHistogram(histClass, "Psi2C_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 100, -2.0, 2.0, VarManager::kPsi2C);
      hm->AddHistogram(histClass, "centrFT0C_Corr2REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR2REF, VarManager::kM11REF);
      hm->AddHistogram(histClass, "centrFT0C_Corr2REFetagap_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR2REFetagap, VarManager::kM11REFetagap);
      hm->AddHistogram(histClass, "centrFT0C_Corr4REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR4REF, VarManager::kM1111REF);
      hm->AddHistogram(histClass, "centrFT0C_Corr2Corr4REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR2CORR4REF, VarManager::kM11M1111REF);
      hm->AddHistogram(histClass, "Run2_centrFT0C_Corr2REF_ev", "", true, 9, std::array<double, 10>{0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}.data(), VarManager::kCentFT0C, 250, std::array<double, 2>{-1.0, 1.0}.data(), VarManager::kCORR2REF, 0, nullptr, -1, "", "", "", VarManager::kCORR2REF, VarManager::kM11REF);
      hm->AddHistogram(histClass, "Run2_centrFT0C_Corr2REFetagap_ev", "", true, 9, std::array<double, 10>{0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}.data(), VarManager::kCentFT0C, 250, std::array<double, 2>{-1.0, 1.0}.data(), VarManager::kCORR2REFetagap, 0, nullptr, -1, "", "", "", VarManager::kCORR2REFetagap, VarManager::kM11REFetagap);
      hm->AddHistogram(histClass, "Run2_centrFT0C_Corr4REF_ev", "", true, 9, std::array<double, 10>{0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}.data(), VarManager::kCentFT0C, 250, std::array<double, 2>{-1.0, 1.0}.data(), VarManager::kCORR4REF, 0, nullptr, -1, "", "", "", VarManager::kCORR4REF, VarManager::kM1111REF);
      hm->AddHistogram(histClass, "Run2_centrFT0C_Corr2Corr4REF_ev", "", true, 9, std::array<double, 10>{0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}.data(), VarManager::kCentFT0C, 250, std::array<double, 2>{-1.0, 1.0}.data(), VarManager::kCORR2CORR4REF, 0, nullptr, -1, "", "", "", VarManager::kCORR2CORR4REF, VarManager::kM11M1111REF);
      hm->AddHistogram(histClass, "centrFT0C_M11REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 1000000.0, VarManager::kM11REF);
      hm->AddHistogram(histClass, "centrFT0C_M11etagap_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 10000000.0, VarManager::kM11REFetagap);
      hm->AddHistogram(histClass, "centrFT0C_M1111REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 100000000000000.0, VarManager::kM1111REF);
      hm->AddHistogram(histClass, "centrFT0C_M11M1111REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 10000000000000000.0, VarManager::kM11M1111REF);
      if (subGroupStr.Contains("cross")) {
        hm->AddHistogram(histClass, "Q1ZNACXX_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACXX);
        hm->AddHistogram(histClass, "Q1ZNACYY_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACYY);
        hm->AddHistogram(histClass, "Q1ZNACYX_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACYX);
        hm->AddHistogram(histClass, "Q1ZNACXY_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 4000, -2, 2, VarManager::kQ1ZNACXY);
        hm->AddHistogram(histClass, "Q2X0A_Q2Y0A", "", false, 500, -10.0, 10.0, VarManager::kQ2X0A, 500, -10.0, 10.0, VarManager::kQ2Y0A);
        hm->AddHistogram(histClass, "Q2X0B_Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B, 500, -10.0, 10.0, VarManager::kQ2Y0B);
        hm->AddHistogram(histClass, "Q2X0C_Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C, 500, -10.0, 10.0, VarManager::kQ2Y0C);
        hm->AddHistogram(histClass, "Q2X0B_Q2Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ2X0B, 500, -10.0, 10.0, VarManager::kQ2Y0C);
        hm->AddHistogram(histClass, "Q2X0C_Q2Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ2X0C, 500, -10.0, 10.0, VarManager::kQ2Y0B);
        hm->AddHistogram(histClass, "Q3X0B_Q3Y0C", "", false, 500, -10.0, 10.0, VarManager::kQ3X0B, 500, -10.0, 10.0, VarManager::kQ3Y0C);
        hm->AddHistogram(histClass, "Q3X0C_Q3Y0B", "", false, 500, -10.0, 10.0, VarManager::kQ3X0C, 500, -10.0, 10.0, VarManager::kQ3Y0B);
        hm->AddHistogram(histClass, "Q2YYAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYAB);
        hm->AddHistogram(histClass, "Q2XXAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXAB);
        hm->AddHistogram(histClass, "Q2XYAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYAB);
        hm->AddHistogram(histClass, "Q2YXAB_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXAB);
        hm->AddHistogram(histClass, "Q2YYAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYAC);
        hm->AddHistogram(histClass, "Q2XXAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXAC);
        hm->AddHistogram(histClass, "Q2XYAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYAC);
        hm->AddHistogram(histClass, "Q2YXAC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXAC);
        hm->AddHistogram(histClass, "Q2YYBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YYBC);
        hm->AddHistogram(histClass, "Q2XXBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XXBC);
        hm->AddHistogram(histClass, "Q2XYBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2XYBC);
        hm->AddHistogram(histClass, "Q2YXBC_Cent", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -1.0, 1.0, VarManager::kQ2YXBC);
        hm->AddHistogram(histClass, "Q2YYAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYAB);
        hm->AddHistogram(histClass, "Q2XXAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXAB);
        hm->AddHistogram(histClass, "Q2XYAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYAB);
        hm->AddHistogram(histClass, "Q2YXAB_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXAB);
        hm->AddHistogram(histClass, "Q2YYAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYAC);
        hm->AddHistogram(histClass, "Q2XXAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXAC);
        hm->AddHistogram(histClass, "Q2XYAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYAC);
        hm->AddHistogram(histClass, "Q2YXAC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXAC);
        hm->AddHistogram(histClass, "Q2YYBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YYBC);
        hm->AddHistogram(histClass, "Q2XXBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XXBC);
        hm->AddHistogram(histClass, "Q2XYBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2XYBC);
        hm->AddHistogram(histClass, "Q2YXBC_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -1.0, 1.0, VarManager::kQ2YXBC);
      }
    }
    if (subGroupStr.Contains("res")) {
      hm->AddHistogram(histClass, "R2SP_TPCFT0A_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_AB);
      hm->AddHistogram(histClass, "R2SP_TPCFT0C_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_AC);
      hm->AddHistogram(histClass, "R2SP_FT0AFT0C_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2SP_BC);
      hm->AddHistogram(histClass, "R3SP_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR3SP);
      hm->AddHistogram(histClass, "R3EP_CentV0M", "", true, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR3EP);
      hm->AddHistogram(histClass, "R2EP_TPCFT0A_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_AB);
      hm->AddHistogram(histClass, "R2EP_TPCFT0C_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_AC);
      hm->AddHistogram(histClass, "R2EP_FT0CFT0A_CentV0M", "", false, 18, 0.0, 90.0, VarManager::kCentVZERO, 500, -10.0, 10.0, VarManager::kR2EP_BC);
      hm->AddHistogram(histClass, "R2SP_TPCFT0A_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_AB);
      hm->AddHistogram(histClass, "R2SP_TPCFT0C_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_AC);
      hm->AddHistogram(histClass, "R2SP_FT0AFT0C_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_BC);
      hm->AddHistogram(histClass, "R2SP_FT0CTPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0CTPCPOS);
      hm->AddHistogram(histClass, "R2SP_FT0CTPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0CTPCNEG);
      hm->AddHistogram(histClass, "R2SP_FT0ATPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0ATPCPOS);
      hm->AddHistogram(histClass, "R2SP_FT0ATPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2SP_FT0ATPCNEG);
      hm->AddHistogram(histClass, "R3SP_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR3SP);
      hm->AddHistogram(histClass, "R2EP_TPCFT0A_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_AB);
      hm->AddHistogram(histClass, "R2EP_TPCFT0C_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_AC);
      hm->AddHistogram(histClass, "R2EP_FT0CFT0A_CentFT0C", "", false, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_BC);
      hm->AddHistogram(histClass, "R2EP_FT0CTPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0CTPCPOS);
      hm->AddHistogram(histClass, "R2EP_FT0CTPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0CTPCNEG);
      hm->AddHistogram(histClass, "R2EP_FT0ATPCPOS_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0ATPCPOS);
      hm->AddHistogram(histClass, "R2EP_FT0ATPCNEG_CentFT0C", "", false, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR2EP_FT0ATPCNEG);
      hm->AddHistogram(histClass, "R3EP_CentFT0C", "", true, 18, 0.0, 90.0, VarManager::kCentFT0C, 500, -10.0, 10.0, VarManager::kR3EP);
    }
    if (subGroupStr.Contains("reso-profile")) {
      hm->AddHistogram(histClass, "Profile_R2SP_TPCFT0A_CentFT0C", "Profile_R2SP_TPCFT0A_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_AB, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_AB);
      hm->AddHistogram(histClass, "Profile_R2SP_TPCFT0C_CentFT0C", "Profile_R2SP_TPCFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_AC, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_AC);
      hm->AddHistogram(histClass, "Profile_R2SP_FT0AFT0C_CentFT0C", "Profile_R2SP_FT0AFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_BC, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_BC);
      hm->AddHistogram(histClass, "Profile_R2EP_TPCFT0A_CentFT0C", "Profile_R2EP_TPCFT0A_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_AB, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_AB);
      hm->AddHistogram(histClass, "Profile_R2EP_TPCFT0C_CentFT0C", "Profile_R2EP_TPCFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_AC, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_AC);
      hm->AddHistogram(histClass, "Profile_R2EP_FT0AFT0C_CentFT0C", "Profile_R2EP_FT0AFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_BC, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_BC);
      hm->AddHistogram(histClass, "Profile_R2SP_Im_TPCFT0A_CentFT0C", "Profile_R2SP_Im_TPCFT0A_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_AB_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_AB_Im);
      hm->AddHistogram(histClass, "Profile_R2SP_Im_TPCFT0C_CentFT0C", "Profile_R2SP_Im_TPCFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_AC_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_AC_Im);
      hm->AddHistogram(histClass, "Profile_R2SP_Im_FT0AFT0C_CentFT0C", "Profile_R2SP_Im_FT0AFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2SP_BC_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2SP_BC_Im);
      hm->AddHistogram(histClass, "Profile_R2EP_Im_TPCFT0A_CentFT0C", "Profile_R2EP_Im_TPCFT0A_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_AB_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_AB_Im);
      hm->AddHistogram(histClass, "Profile_R2EP_Im_TPCFT0C_CentFT0C", "Profile_R2EP_Im_TPCFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_AC_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_AC_Im);
      hm->AddHistogram(histClass, "Profile_R2EP_Im_FT0AFT0C_CentFT0C", "Profile_R2EP_Im_FT0AFT0C_CentFT0C", true, 90, 0.0, 90.0, VarManager::kCentFT0C, 500, 0.0, 5.0, VarManager::kR2EP_BC_Im, 0, 0, 0, -1, "", "", "", -1, VarManager::kWR2EP_BC_Im);
    }
    if (subGroupStr.Contains("filter")) {
      hm->AddHistogram(histClass, "IsDoubleGap", "Is double gap", false, 2, -0.5, 1.5, VarManager::kIsDoubleGap);
      hm->AddHistogram(histClass, "IsSingleGapA", "Is single gap on side A", false, 2, -0.5, 1.5, VarManager::kIsSingleGapA);
      hm->AddHistogram(histClass, "IsSingleGapC", "Is single gap on side C", false, 2, -0.5, 1.5, VarManager::kIsSingleGapC);
      hm->AddHistogram(histClass, "IsITSUPCMode", "UPC settings used", false, 2, -0.5, 1.5, VarManager::kIsITSUPCMode);
      hm->AddHistogram(histClass, "IsITSUPCMode_IsSingleGap", "UPC settings used vs Is single gap", false, 2, -0.5, 1.5, VarManager::kIsITSUPCMode, 2, -0.5, 1.5, VarManager::kIsSingleGap);
    }
    if (subGroupStr.Contains("zdc")) {
      hm->AddHistogram(histClass, "energyCommonZNA_energyCommonZNC", "Common ZNA energy vs common ZNC energy", false, 1050, -10.0, 200.0, VarManager::kEnergyCommonZNA, 1050, -10.0, 200.0, VarManager::kEnergyCommonZNC);
      hm->AddHistogram(histClass, "energyCommonZNA_energyCommonZNC_lowRange", "Common ZNA energy vs common ZNC energy", false, 220, -2.0, 20.0, VarManager::kEnergyCommonZNA, 220, -2.0, 20.0, VarManager::kEnergyCommonZNC);
    }
  } // end "event"

  if (!groupStr.CompareTo("two-collisions")) {
    hm->AddHistogram(histClass, "DeltaZ", "z_{1} - z_{2}", false, 400, -20., 20., VarManager::kTwoEvDeltaZ);
    hm->AddHistogram(histClass, "DeltaZ_Z1", "z_{1} - z_{2} vs z_{1}", false, 24, -12., 12., VarManager::kTwoEvPosZ1, 300, -15., 15., VarManager::kTwoEvDeltaZ);
    hm->AddHistogram(histClass, "DeltaR", "r_{1} - r_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaR);
    hm->AddHistogram(histClass, "DeltaX", "x_{1} - x_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaX);
    hm->AddHistogram(histClass, "DeltaY", "y_{1} - y_{2}", false, 200, -0.1, 0.1, VarManager::kTwoEvDeltaY);
    hm->AddHistogram(histClass, "Z1vsZ2", "z_{1} vs z_{2}", false, 120, -12.0, 12.0, VarManager::kTwoEvPosZ1, 120, -12.0, 12.0, VarManager::kTwoEvPosZ2);
    hm->AddHistogram(histClass, "R1vsR2", "r_{1} vs r_{2}", false, 100, -0.1, 0.1, VarManager::kTwoEvPosR1, 100, -0.1, 0.1, VarManager::kTwoEvPosR2);
    hm->AddHistogram(histClass, "NContrib1vs2", "n.contrib 1 vs 2", false, 100, 0.0, 100.0, VarManager::kTwoEvPVcontrib1, 100, 0.0, 100.0, VarManager::kTwoEvPVcontrib2);
  }

  if (!groupStr.CompareTo("track")) {
    hm->AddHistogram(histClass, "Pt", "p_{T} distribution", false, 2000, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta);
    hm->AddHistogram(histClass, "Phi", "#varphi distribution", false, 500, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
    hm->AddHistogram(histClass, "Phi_Pt", "#varphi distribution", false, 50, 0.0, 10.0, VarManager::kPt, 720, 0.0, o2::constants::math::TwoPI, VarManager::kPhi);
    hm->AddHistogram(histClass, "IsPVcontrib_pt", "is PV contributor vs pt", false, 50, 0.0, 50.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kPVContributor);
    hm->AddHistogram(histClass, "IsPVcontrib_pt_prof", "is PV contributor vs pt", true, 50, 0.0, 50.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kPVContributor);
    if (subGroupStr.Contains("ambiguity") && !subGroupStr.Contains("muon")) {
      hm->AddHistogram(histClass, "AmbiguityInBunch", "in bunch collision ambiguity", false, 10, 0., 10., VarManager::kBarrelNAssocsInBunch);
      hm->AddHistogram(histClass, "AmbiguityInBunch_pt", "in bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kBarrelNAssocsInBunch);
      hm->AddHistogram(histClass, "AmbiguityOutOfBunch", "out of bunch collision ambiguity", false, 10, 0., 10., VarManager::kBarrelNAssocsOutOfBunch);
      hm->AddHistogram(histClass, "AmbiguityOutOfBunch_pt", "out of bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kBarrelNAssocsOutOfBunch);
    }
    if (subGroupStr.Contains("cent")) {
      hm->AddHistogram(histClass, "Pt_CentFT0C", "p_{T} distribution", false, 2000, 0.0, 20.0, VarManager::kPt, 20, 0.0, 100.0, VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "Eta_CentFT0C", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta, 20, 0.0, 100.0, VarManager::kCentFT0C);
      hm->AddHistogram(histClass, "Phi_CentFT0C", "#varphi distribution", false, 500, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi, 20, 0.0, 100.0, VarManager::kCentFT0C);
    }
    if (subGroupStr.Contains("kine")) {
      hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kEta, 200, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
      hm->AddHistogram(histClass, "Eta_Pt", "", false, 20, -1.0, 1.0, VarManager::kEta, 100, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Eta_VtxZ", "", false, 100, -1.0, 1.0, VarManager::kEta, 300, -15.0, 15.0, VarManager::kVtxZ);
      hm->AddHistogram(histClass, "Px", "p_{x} distribution", false, 200, 0.0, 20.0, VarManager::kPx);
      hm->AddHistogram(histClass, "Py", "p_{y} distribution", false, 200, 0.0, 20.0, VarManager::kPy);
      hm->AddHistogram(histClass, "Pz", "p_{z} distribution", false, 200, 0.0, 20.0, VarManager::kPz);
      hm->AddHistogram(histClass, "InvPt", "1/p_{T} distribution", false, 1500, 0.0, 15.0, VarManager::kInvPt);
    }
    if (subGroupStr.Contains("time")) {
      hm->AddHistogram(histClass, "TrackTime", "", false, 800, -10000.0, 10000.0, VarManager::kTrackTime);
      hm->AddHistogram(histClass, "TrackTimeRes", "", false, 400, 0.0, 1000.0, VarManager::kTrackTimeRes);
      hm->AddHistogram(histClass, "TrackTimeResRelative", "", false, 100, 0.0, 100.0, VarManager::kTrackTimeResRelative);
      hm->AddHistogram(histClass, "TOFExpMom", "", false, 100, 0.0, 10.0, VarManager::kTOFExpMom);
    }
    if (subGroupStr.Contains("map")) {
      hm->AddHistogram(histClass, "HasITSandTPC", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTPC);
      hm->AddHistogram(histClass, "HasITSandTOF", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTOF);
      hm->AddHistogram(histClass, "HasITSandTRD", "", false, 2, -0.5, 1.5, VarManager::kHasITS, 2, -0.5, 1.5, VarManager::kHasTRD);
      hm->AddHistogram(histClass, "HasTPCandTOF", "", false, 2, -0.5, 1.5, VarManager::kHasTPC, 2, -0.5, 1.5, VarManager::kHasTOF);
    }
    if (subGroupStr.Contains("its")) {
      hm->AddHistogram(histClass, "ITSncls", "Number of cluster in ITS", false, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSncls_vs_BC", "Average Number of cluster in ITS vs BC", true, 3564, 0.0, 3564, VarManager::kBCOrbit, 2, 0.0, 7.0, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2", "ITS chi2", false, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit", "", false, 2, -0.5, 1.5, VarManager::kIsITSrefit);
      hm->AddHistogram(histClass, "IsSPDany", "", false, 2, -0.5, 1.5, VarManager::kIsSPDany);
      hm->AddHistogram(histClass, "IsSPDfirst", "", false, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
      hm->AddHistogram(histClass, "ITSncls_vsTimeFromSOR", "Number of cluster in ITS vs time from SOR", true, 10000, 0.0, 1000, VarManager::kTimeFromSOR, 8, -0.5, 7.5, VarManager::kITSncls);
      if (subGroupStr.Contains("cluster")) {
        hm->AddHistogram(histClass, "ITSClusterMap", "", false, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_pin", "ITSClustermap vs pin", false, 200, 0.0, 20.0, VarManager::kPin, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_SignedPin", "ITSClustermap vs SignedPin", false, 400, -20.0, 20.0, VarManager::kSignedPin, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_pt", "ITSClustermap vs pt", false, 200, 0.0, 20.0, VarManager::kPt, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_eta", "ITSClustermap vs eta", false, 100, -1.0, 1.0, VarManager::kEta, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "ITSClustermap_vs_phi", "ITSClustermap vs phi", false, 315, 0.0, 6.3, VarManager::kPhi, 128, -0.5, 127.5, VarManager::kITSClusterMap);
        hm->AddHistogram(histClass, "SignedPin_P_ITSMap", "SignedPin vs P vs ITSMap", false, 400, -20.0, 20.0, VarManager::kSignedPin, 200, 0.0, 20.0, VarManager::kP, 2, -0.5, 1.5, VarManager::kHasITS);
      }
      if (subGroupStr.Contains("cent")) {
        hm->AddHistogram(histClass, "ITSncls_CentFT0C", "Number of cluster in ITS", false, 8, -0.5, 7.5, VarManager::kITSncls, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "ITSchi2_CentFT0C", "ITS chi2", false, 100, 0.0, 50.0, VarManager::kITSchi2, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "ITSClusterMap_CentFT0C", "", false, 128, -0.5, 127.5, VarManager::kITSClusterMap, 20, 0.0, 100.0, VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("clssize")) {
        hm->AddHistogram(histClass, "ITSclssize_p", "Mean ITS cluster size vs P", false, 200, 0.0, 10.0, VarManager::kP, 150, 0.0, 15., VarManager::kITSmeanClsSize);
      }
    }
    if (subGroupStr.Contains("itsvspt")) {
      hm->AddHistogram(histClass, "ITSncls_Pt", "Number of cluster in ITS vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2_Pt", "ITS chi2 vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsITSrefit);
      hm->AddHistogram(histClass, "IsSPDany_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsSPDany);
      hm->AddHistogram(histClass, "IsSPDfirst_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
      hm->AddHistogram(histClass, "ITSncls_signedPt", "Number of cluster in ITS vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 8, -0.5, 7.5, VarManager::kITSncls);
      hm->AddHistogram(histClass, "ITSchi2_signedPt", "ITS chi2 vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 100, 0.0, 50.0, VarManager::kITSchi2);
      hm->AddHistogram(histClass, "IsITSrefit_signedPt", "", false, 400, -10.0, 10.0, VarManager::kSignedPt, 2, -0.5, 1.5, VarManager::kIsITSrefit);
    }
    if (subGroupStr.Contains("tpc")) {
      hm->AddHistogram(histClass, "TPCncls", "Number of cluster in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_vsTimeFromSOR", "Number of cluster in TPC vs time from SOR", true, 10000, 0.0, 1000., VarManager::kTimeFromSOR, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_Phi", "Number of cluster in TPC vs #varphi", true, 720, 0.0, o2::constants::math::TwoPI, VarManager::kPhi, 10, 0.0, 200.0, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCncls_PhiPt", "Number of cluster in TPC vs p_{T} and #varphi", true, 20, 0.0, 10.0, VarManager::kPt, 720, 0.0, o2::constants::math::TwoPI, VarManager::kPhi, 10, 0.0, 200.0, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCnclsCR", "Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "TPCncls_TPCnclsCR", "Number of TPC cluster vs Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "IsTPCrefit", "", false, 2, -0.5, 1.5, VarManager::kIsTPCrefit);
      hm->AddHistogram(histClass, "IsGoldenChi2", "", false, 2, -0.5, 1.5, VarManager::kIsGoldenChi2);
      hm->AddHistogram(histClass, "TPCchi2", "TPC chi2", false, 100, 0.0, 10.0, VarManager::kTPCchi2);
      hm->AddHistogram(histClass, "pin_vs_p", "", false, 400, -20.0, 20.0, VarManager::kSignedPin, 200, 0.0, 20, VarManager::kP);
      if (subGroupStr.Contains("cent")) {
        hm->AddHistogram(histClass, "TPCncls_CentFT0C", "Number of cluster in TPC", false, 160, -0.5, 159.5, VarManager::kTPCncls, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "TPCnclsCR_CentFT0C", "Number of crossed rows in TPC", false, 160, -0.5, 159.5, VarManager::kTPCnclsCR, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "TPCchi2_CentFT0C", "TPC chi2", false, 100, 0.0, 10.0, VarManager::kTPCchi2, 20, 0.0, 100.0, VarManager::kCentFT0C);
      }
    }
    if (subGroupStr.Contains("tpcvspt")) {
      hm->AddHistogram(histClass, "TPCncls_Pt", "Number of cluster in TPC vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCnclsCR_Pt", "Number of crossed rows in TPC vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "IsTPCrefit_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsTPCrefit);
      hm->AddHistogram(histClass, "IsGoldenChi2_Pt", "", false, 200, 0.0, 10.0, VarManager::kPt, 2, -0.5, 1.5, VarManager::kIsGoldenChi2);
      hm->AddHistogram(histClass, "TPCchi2_Pt", "TPC chi2 vs Pt", false, 200, 0.0, 10.0, VarManager::kPt, 100, 0.0, 10.0, VarManager::kTPCchi2);
      hm->AddHistogram(histClass, "TPCncls_signedPt", "Number of cluster in TPC vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 160, -0.5, 159.5, VarManager::kTPCncls);
      hm->AddHistogram(histClass, "TPCnclsCR_signedPt", "Number of crossed rows in TPC vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 160, -0.5, 159.5, VarManager::kTPCnclsCR);
      hm->AddHistogram(histClass, "TPCchi2_signedPt", "TPC chi2 vs signed Pt", false, 400, -10.0, 10.0, VarManager::kSignedPt, 100, 0.0, 10.0, VarManager::kTPCchi2);
    }
    if (subGroupStr.Contains("tpcpid")) {
      if (subGroupStr.Contains("tpcpid_fine")) {
        // fine binning for pIN: steps in 10 MeV/c from 0 to 1 GeV/c and 100 MeV/c up to 10 GeV/c
        double pIN_bins[281];
        for (int i = 0; i <= 200; i++)
          pIN_bins[i] = 0.01 * i;
        for (int i = 1; i <= 80; i++)
          pIN_bins[200 + i] = 2. + 0.1 * i;
        int nbins_pIN = sizeof(pIN_bins) / sizeof(*pIN_bins) - 1;

        double TPCdEdx_bins[201];
        for (int i = 0; i <= 200; i++)
          TPCdEdx_bins[i] = i;
        int nbins_TPCdEdx = sizeof(TPCdEdx_bins) / sizeof(*TPCdEdx_bins) - 1;

        double nSigma_bins[101];
        for (int i = 0; i <= 100; i++)
          nSigma_bins[i] = -5. + 0.1 * i;
        int nbins_nSigma = sizeof(nSigma_bins) / sizeof(*nSigma_bins) - 1;

        hm->AddHistogram(histClass, "TPCdedx_pIN", "TPC dE/dx vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_TPCdEdx, TPCdEdx_bins, VarManager::kTPCsignal);
        hm->AddHistogram(histClass, "TPCnSigEle_pIN", "TPC n-#sigma(e) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigPi_pIN", "TPC n-#sigma(#pi) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigKa_pIN", "TPC n-#sigma(K) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaKa);
        hm->AddHistogram(histClass, "TPCnSigPr_pIN", "TPC n-#sigma(p) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPr);
        if (subGroupStr.Contains("tpcpid_fine_Corr")) {
          hm->AddHistogram(histClass, "TPCnSigEl_Corr_pIN", "TPC n-#sigma(e) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaEl_Corr);
          hm->AddHistogram(histClass, "TPCnSigPi_Corr_pIN", "TPC n-#sigma(#pi) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPi_Corr);
          hm->AddHistogram(histClass, "TPCnSigPr_Corr_pIN", "TPC n-#sigma(p) Corr. vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTPCnSigmaPr_Corr);
        }
      } else {
        hm->AddHistogram(histClass, "TPCdedx_pIN", "TPC dE/dx vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 200, 0.0, 200., VarManager::kTPCsignal);
        hm->AddHistogram(histClass, "TPCnSigEle_pIN", "TPC n-#sigma(e) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupancy", "TPC n-#sigma(e) vs occupancy", false, 200, 0., 20000., VarManager::kTrackOccupancyInTimeRange, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_timeFromSOR", "TPC n-#sigma(e) vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCcontribLongA", "TPC n-#sigma(e) vs pileup n-contrib, long time range A-side", false, 20, 0.0, 10000.0, VarManager::kNTPCcontribLongA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmeanTimeLongA", "TPC n-#sigma(e) vs pileup mean time, long time range, A-side", false, 20, -100.0, 100.0, VarManager::kNTPCmeanTimeLongA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmedianTimeLongA", "TPC n-#sigma(e) vs pileup mean time, long time range, A-side", false, 20, -100.0, 100.0, VarManager::kNTPCmedianTimeLongA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCcontribShortA", "TPC n-#sigma(e) vs pileup n-contrib, short time range A-side", false, 50, 0.0, 10000.0, VarManager::kNTPCcontribShortA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmeanTimeShortA", "TPC n-#sigma(e) vs pileup mean time, short time range, A-side", false, 20, -15.0, 15.0, VarManager::kNTPCmeanTimeShortA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmedianTimeShortA", "TPC n-#sigma(e) vs pileup mean time, short time range, A-side", false, 20, -15.0, 15.0, VarManager::kNTPCmedianTimeShortA, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCcontribLongC", "TPC n-#sigma(e) vs pileup n-contrib, long time range C-side", false, 20, 0.0, 10000.0, VarManager::kNTPCcontribLongC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmeanTimeLongC", "TPC n-#sigma(e) vs pileup mean time, long time range, C-side", false, 20, -100.0, 100.0, VarManager::kNTPCmeanTimeLongC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmedianTimeLongC", "TPC n-#sigma(e) vs pileup mean time, long time range, C-side", false, 20, -100.0, 100.0, VarManager::kNTPCmedianTimeLongC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCcontribShortC", "TPC n-#sigma(e) vs pileup n-contrib, short time range C-side", false, 50, 0.0, 10000.0, VarManager::kNTPCcontribShortC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmeanTimeShortC", "TPC n-#sigma(e) vs pileup mean time, short time range, C-side", false, 20, -15.0, 15.0, VarManager::kNTPCmeanTimeShortC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigEle_occupTPCmedianTimeShortC", "TPC n-#sigma(e) vs pileup mean time, short time range, C-side", false, 20, -15.0, 15.0, VarManager::kNTPCmedianTimeShortC, 200, -5.0, 5.0, VarManager::kTPCnSigmaEl);
        hm->AddHistogram(histClass, "TPCnSigPi_pIN", "TPC n-#sigma(#pi) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_timeFromSOR", "TPC n-#sigma(#pi) vs time from SOR", true, 1000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_eta", "TPC n-#sigma(#pi) vs #eta", false, 20, -1.0, 1.0, VarManager::kEta, 200, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_etaPin_prof", "<TPC n-#sigma(#pi)> vs (#eta,p_{IN}), --s--", true, 20, -1.0, 1.0, VarManager::kEta, 20, 0.0, 10.0, VarManager::kPin, 10, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_etaCent_prof", "<TPC n-#sigma(#pi)> vs (#eta,cent), --s--", true, 20, -1.0, 1.0, VarManager::kEta, 20, 0.0, 100.0, VarManager::kCentFT0C, 10, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_etaContrib_prof", "<TPC n-#sigma(#pi)> vs (#eta,n-contrib real), --s--", true, 20, -1.0, 1.0, VarManager::kEta, 20, 0.0, 4000.0, VarManager::kVtxNcontribReal, 10, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_centFT0C", "TPC n-#sigma(#pi) vs centrality", false, 20, 0.0, 100.0, VarManager::kCentFT0C, 200, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_vtxContrib", "TPC n-#sigma(#pi) vs vtx. contrib real", false, 50, 0.0, 4000.0, VarManager::kVtxNcontribReal, 200, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigPi_occupancy", "TPC n-#sigma(#pi) vs occupancy", false, 200, 0., 20000., VarManager::kTrackOccupancyInTimeRange, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi);
        hm->AddHistogram(histClass, "TPCnSigKa_pIN", "TPC n-#sigma(K) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaKa);
        hm->AddHistogram(histClass, "TPCnSigPr_pIN", "TPC n-#sigma(p) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr);
        hm->AddHistogram(histClass, "TPCnSigPr_timeFromSOR", "TPC n-#sigma(p) vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, -5.0, 5.0, VarManager::kTPCnSigmaPr);
        hm->AddHistogram(histClass, "TPCnSigPr_occupancy", "TPC n-#sigma(p) vs. occupancy", false, 200, 0., 20000., VarManager::kTrackOccupancyInTimeRange, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr);
        if (subGroupStr.Contains("tpcpid_Corr")) {
          hm->AddHistogram(histClass, "TPCnSigEl_Corr_pIN", "TPC n-#sigma(e) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl_Corr);
          hm->AddHistogram(histClass, "TPCnSigPi_Corr_pIN", "TPC n-#sigma(#pi) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPi_Corr);
          hm->AddHistogram(histClass, "TPCnSigKa_Corr_pIN", "TPC n-#sigma(K) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaKa_Corr);
          hm->AddHistogram(histClass, "TPCnSigPr_Corr_pIN", "TPC n-#sigma(p) Corr. vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTPCnSigmaPr_Corr);
        }
      }
    }
    if (subGroupStr.Contains("postcalib")) {
      const int kNvarsPID = 4;
      const int kTPCnsigmaNbins = 70;
      double tpcNsigmaBinLims[kTPCnsigmaNbins + 1];
      for (int i = 0; i <= kTPCnsigmaNbins; ++i)
        tpcNsigmaBinLims[i] = -7.0 + 0.2 * i;

      const int kPinEleNbins = 20;
      double pinEleBinLims[kPinEleNbins + 1] = {0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0};

      const int kEtaNbins = 9;
      double etaBinLimsI[kEtaNbins + 1] = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};

      const int kTPCnClusterbins = 16;
      double tpcNclusterBinLims[kTPCnClusterbins + 1];
      for (int i = 0; i <= kTPCnClusterbins; ++i)
        tpcNclusterBinLims[i] = 10 * i;

      TArrayD nSigBinLimits[kNvarsPID];
      nSigBinLimits[0] = TArrayD(kTPCnsigmaNbins + 1, tpcNsigmaBinLims);
      nSigBinLimits[1] = TArrayD(kTPCnClusterbins + 1, tpcNclusterBinLims);
      nSigBinLimits[2] = TArrayD(kPinEleNbins + 1, pinEleBinLims);
      nSigBinLimits[3] = TArrayD(kEtaNbins + 1, etaBinLimsI);

      const int kPinKaNbins = 15;
      double pinKaBinLims[kPinKaNbins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0};

      TArrayD nSigBinLimitsKa[kNvarsPID];
      nSigBinLimitsKa[0] = TArrayD(kTPCnsigmaNbins + 1, tpcNsigmaBinLims);
      nSigBinLimitsKa[1] = TArrayD(kTPCnClusterbins + 1, tpcNclusterBinLims);
      nSigBinLimitsKa[2] = TArrayD(kPinKaNbins + 1, pinKaBinLims);
      nSigBinLimitsKa[3] = TArrayD(kEtaNbins + 1, etaBinLimsI);

      if (subGroupStr.Contains("electron")) {
        int varsPIDnSigEle[kNvarsPID] = {VarManager::kTPCnSigmaEl, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigEle_Corr[kNvarsPID] = {VarManager::kTPCnSigmaEl_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCelectron", "TPC n_{#sigma}(e) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigEle, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCelectron_Corr", "TPC n_{#sigma}^{Corr}(e) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigEle_Corr, nSigBinLimits);
      }
      if (subGroupStr.Contains("pion")) {
        int varsPIDnSigPion[kNvarsPID] = {VarManager::kTPCnSigmaPi, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigPion_Corr[kNvarsPID] = {VarManager::kTPCnSigmaPi_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCpion", "TPC n_{#sigma}(pion) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigPion, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCpion_Corr", "TPC n_{#sigma}^{Corr}(pion) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigPion_Corr, nSigBinLimits);
      }
      if (subGroupStr.Contains("kaon")) {
        int varsPIDnSigKaon[kNvarsPID] = {VarManager::kTPCnSigmaKa, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigKaon_Corr[kNvarsPID] = {VarManager::kTPCnSigmaKa_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCkaon", "TPC n_{#sigma}(kaon) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigKaon, nSigBinLimitsKa);
        hm->AddHistogram(histClass, "nSigmaTPCkaon_Corr", "TPC n_{#sigma}^{Corr}(kaon) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigKaon_Corr, nSigBinLimitsKa);
      }
      if (subGroupStr.Contains("proton")) {
        int varsPIDnSigProton[kNvarsPID] = {VarManager::kTPCnSigmaPr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        int varsPIDnSigProton_Corr[kNvarsPID] = {VarManager::kTPCnSigmaPr_Corr, VarManager::kTPCncls, VarManager::kPin, VarManager::kEta};
        hm->AddHistogram(histClass, "nSigmaTPCproton", "TPC n_{#sigma}(proton) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigProton, nSigBinLimits);
        hm->AddHistogram(histClass, "nSigmaTPCproton_Corr", "TPC n_{#sigma}^{Corr}(proton) Vs normNcluster Vs Pin Vs Eta", kNvarsPID, varsPIDnSigProton_Corr, nSigBinLimits);
      }
    }
    if (subGroupStr.Contains("tofpid")) {
      if (subGroupStr.Contains("tofpid_fine")) {
        // fine binning for pIN: steps in 10 MeV/c from 0 to 1 GeV/c and 100 MeV/c up to 10 GeV/c
        double pIN_bins[281];
        for (int i = 0; i <= 200; i++)
          pIN_bins[i] = 0.01 * i;
        for (int i = 1; i <= 80; i++)
          pIN_bins[200 + i] = 2. + 0.1 * i;
        int nbins_pIN = sizeof(pIN_bins) / sizeof(*pIN_bins) - 1;

        double TOFbeta_bins[241];
        for (int i = 0; i <= 240; i++)
          TOFbeta_bins[i] = 0.005 * i;
        int nbins_TOFbeta = sizeof(TOFbeta_bins) / sizeof(*TOFbeta_bins) - 1;

        double nSigma_bins[101];
        for (int i = 0; i <= 100; i++)
          nSigma_bins[i] = -5. + 0.1 * i;
        int nbins_nSigma = sizeof(nSigma_bins) / sizeof(*nSigma_bins) - 1;

        hm->AddHistogram(histClass, "TOFbeta_pIN", "TOF #beta vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_TOFbeta, TOFbeta_bins, VarManager::kTOFbeta);
        hm->AddHistogram(histClass, "TOFnSigEle_pIN", "TOF n-#sigma(e) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaEl);
        hm->AddHistogram(histClass, "TOFnSigPi_pIN", "TOF n-#sigma(#pi) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaPi);
        hm->AddHistogram(histClass, "TOFnSigKa_pIN", "TOF n-#sigma(K) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaKa);
        hm->AddHistogram(histClass, "TOFnSigPr_pIN", "TOF n-#sigma(p) vs pIN", false, nbins_pIN, pIN_bins, VarManager::kPin, nbins_nSigma, nSigma_bins, VarManager::kTOFnSigmaPr);
      } else {
        hm->AddHistogram(histClass, "TOFbeta_pIN", "TOF #beta vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 240, 0.0, 1.2, VarManager::kTOFbeta);
        hm->AddHistogram(histClass, "TOFnSigEle_pIN", "TOF n-#sigma(e) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaEl);
        hm->AddHistogram(histClass, "TOFnSigPi_pIN", "TOF n-#sigma(#pi) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaPi);
        hm->AddHistogram(histClass, "TOFnSigKa_pIN", "TOF n-#sigma(K) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaKa);
        hm->AddHistogram(histClass, "TOFnSigPr_pIN", "TOF n-#sigma(p) vs pIN", false, 100, 0.0, 10.0, VarManager::kPin, 100, -5.0, 5.0, VarManager::kTOFnSigmaPr);
      }
    }
    if (subGroupStr.Contains("pidcorre")) {
      const int kNvarsPID = 3;
      const int kNbins_pIN = 169;
      double pIN_bins[kNbins_pIN + 1];
      for (int i = 0; i <= 140; i++)
        pIN_bins[i] = 0.01 * i + 0.1;
      for (int i = 1; i <= 15; i++)
        pIN_bins[140 + i] = 1.5 + 0.1 * i;
      for (int i = 1; i <= 14; i++)
        pIN_bins[155 + i] = 3. + 0.5 * i;

      const int kNbins_pINmore = 135;
      double pIN_binsmore[kNbins_pINmore + 1];
      for (int i = 0; i <= 120; i++)
        pIN_binsmore[i] = 0.01 * i + 0.3;
      for (int i = 1; i <= 10; i++)
        pIN_binsmore[120 + i] = 1.5 + 0.2 * i;
      pIN_binsmore[131] = 4.;
      pIN_binsmore[132] = 5.;
      pIN_binsmore[133] = 6.;
      pIN_binsmore[134] = 8.;
      pIN_binsmore[135] = 10.;

      const int kNbins_nSigma = 100;
      double nSigma_bins[kNbins_nSigma + 1];
      for (int i = 0; i <= kNbins_nSigma; i++)
        nSigma_bins[i] = -5. + 0.1 * i;

      const int kNbins_nSigmamore = 50;
      double nSigma_binsmore[kNbins_nSigmamore + 1];
      for (int i = 0; i <= kNbins_nSigmamore; i++)
        nSigma_binsmore[i] = -5. + 0.2 * i;

      const int kNbins_nSigmagrob = 24;
      double nSigma_binsgrob[kNbins_nSigmagrob + 1];
      for (int i = 0; i <= kNbins_nSigmagrob; i++)
        nSigma_binsgrob[i] = -6. + 0.5 * i;

      const int kNbins_TOFbeta = 120;
      double TOFbeta_bins[kNbins_TOFbeta + 1];
      for (int i = 0; i <= kNbins_TOFbeta; i++)
        TOFbeta_bins[i] = 0.01 * i;

      const int kNbins_TPCdEdx = 140;
      double TPCdEdx_bins[kNbins_TPCdEdx + 1];
      for (int i = 0; i <= kNbins_TPCdEdx; i++)
        TPCdEdx_bins[i] = i + 20;

      TArrayD nSigmaBinLimits[kNvarsPID];
      nSigmaBinLimits[0] = TArrayD(kNbins_pIN + 1, pIN_bins);
      nSigmaBinLimits[1] = TArrayD(kNbins_nSigma + 1, nSigma_bins);
      nSigmaBinLimits[2] = TArrayD(kNbins_nSigma + 1, nSigma_bins);

      TArrayD nSignalBinLimits[kNvarsPID];
      nSignalBinLimits[0] = TArrayD(kNbins_pIN + 1, pIN_bins);
      nSignalBinLimits[1] = TArrayD(kNbins_TPCdEdx + 1, TPCdEdx_bins);
      nSignalBinLimits[2] = TArrayD(kNbins_TOFbeta + 1, TOFbeta_bins);

      int varsPIDnSignal[kNvarsPID] = {VarManager::kPin, VarManager::kTPCsignal, VarManager::kTOFbeta};
      hm->AddHistogram(histClass, "nSignalTPCTOF", "", kNvarsPID, varsPIDnSignal, nSignalBinLimits);

      if (subGroupStr.Contains("more")) {
        const int kNvarsPIDmore = 4;
        TArrayD nSigmaBinLimitsmore[kNvarsPIDmore];
        nSigmaBinLimitsmore[0] = TArrayD(kNbins_pINmore + 1, pIN_binsmore);
        nSigmaBinLimitsmore[1] = TArrayD(kNbins_nSigmamore + 1, nSigma_binsmore);
        nSigmaBinLimitsmore[2] = TArrayD(kNbins_nSigmagrob + 1, nSigma_binsgrob);
        nSigmaBinLimitsmore[3] = TArrayD(kNbins_nSigmamore + 1, nSigma_binsmore);
        int varsPIDnSigmamore[kNvarsPIDmore] = {VarManager::kPin, VarManager::kTPCnSigmaEl, VarManager::kTPCnSigmaPi, VarManager::kTOFnSigmaEl};
        hm->AddHistogram(histClass, "nSigmaTPCTOF", "", kNvarsPIDmore, varsPIDnSigmamore, nSigmaBinLimitsmore);
      } else {
        int varsPIDnSigma[kNvarsPID] = {VarManager::kPin, VarManager::kTPCnSigmaEl, VarManager::kTOFnSigmaEl};
        hm->AddHistogram(histClass, "nSigmaTPCTOF", "", kNvarsPID, varsPIDnSigma, nSigmaBinLimits);
      }
    }
    if (subGroupStr.Contains("runbyrun")) {
      hm->AddHistogram(histClass, "TPCncls_Run", "Number of cluster in TPC vs RunNumber", true, 1, -0.5, 0.5, VarManager::kRunNo, 160, -0.5, 159.5, VarManager::kTPCncls, 1, 0., 1., VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "TPCdEdx_Run", "TPCdEdx vs RunNumber", true, 1, -0.5, 0.5, VarManager::kRunNo, 300, 0., 300., VarManager::kTPCsignal, 1, 0., 1., VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "TPCchi2_Run", "TPCchi2 vs RunNumber", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, 0., 10., VarManager::kTPCchi2, 1, 0., 1., VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "Pt_Run", "p_{T} distribution", true, 1, -0.5, 0.5, VarManager::kRunNo, 2000, 0.0, 20.0, VarManager::kPt, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "ITSncls_Run", "Number of cluster in ITS", true, 1, -0.5, 0.5, VarManager::kRunNo, 8, -0.5, 7.5, VarManager::kITSncls, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "ITSchi2_Run", "ITS chi2", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, 0.0, 50.0, VarManager::kITSchi2, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "TPCnSigEle_Run", "TPC n-#sigma(e)", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, -5.0, 5.0, VarManager::kTPCnSigmaEl, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "DCAxy_Run", "DCA_{xy}", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, -1.0, 1.0, VarManager::kTrackDCAxy, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
      hm->AddHistogram(histClass, "DCAz_Run", "DCA_{z}", true, 1, -0.5, 0.5, VarManager::kRunNo, 100, -1.0, 1.0, VarManager::kTrackDCAz, 1, 0, 1, VarManager::kNothing, "", "", "", -1, -1, false, true);
    }
    if (subGroupStr.Contains("dca")) {
      hm->AddHistogram(histClass, "DCAxy", "DCA_{xy}", false, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "DCAxy_vsTimeFromSOR", "DCA_{xy} vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, -2.0, 2.0, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "DCAz", "DCA_{z}", false, 800, -4.0, 4.0, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "DCAz_vsTimeFromSOR", "DCA_{z} vs time from SOR", true, 10000, 0.0, 1000.0, VarManager::kTimeFromSOR, 10, -2.0, 2.0, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "DCAsigXY", "DCA_{XY} [#sigma]", false, 200, -20.0, 20.0, VarManager::kTrackDCAsigXY);
      hm->AddHistogram(histClass, "DCAsigZ", "DCA_{Z} [#sigma]", false, 200, -20.0, 20.0, VarManager::kTrackDCAsigZ);
      hm->AddHistogram(histClass, "DCAxy_DCAz", "DCA_{xy}  vs DCA_{z}", false, 200, -4.0, 4.0, VarManager::kTrackDCAxy, 200, -4.0, 4.0, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "DCAsigXY_DCAsigZ", "DCA_{XY} [#sigma] vs DCA_{Z} [#sigma]", false, 200, -20.0, 20.0, VarManager::kTrackDCAsigXY, 200, -20.0, 20.0, VarManager::kTrackDCAsigZ);
      if (subGroupStr.Contains("pt")) {
        hm->AddHistogram(histClass, "Pt_DCAxy", "p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "Pt_DCAz", "p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 800, -4.0, 4.0, VarManager::kTrackDCAz);
        hm->AddHistogram(histClass, "Pt_DCAsigXY", "p_{T} vs DCA_{XY} [#sigma]", false, 200, 0.0, 20.0, VarManager::kPt, 200, -20.0, 20.0, VarManager::kTrackDCAsigXY); // JJ:edit
        hm->AddHistogram(histClass, "Pt_DCAsigZ", "p_{T} vs DCA_{Z} [#sigma]", false, 200, 0.0, 20.0, VarManager::kPt, 200, -20.0, 20.0, VarManager::kTrackDCAsigZ);
        hm->AddHistogram(histClass, "Pt_DCAresXY", "p_{T} vs #DeltaDCA_{XY}", false, 200, 0.0, 10.0, VarManager::kPt, 100, -0.03, 0.03, VarManager::kTrackDCAresXY);
        hm->AddHistogram(histClass, "Pt_DCAresZ", "p_{T} vs #DeltaDCA_{Z}", false, 200, 0.0, 10.0, VarManager::kPt, 100, -0.03, 0.03, VarManager::kTrackDCAresZ);
      }
      if (subGroupStr.Contains("eta")) {
        hm->AddHistogram(histClass, "Eta_DCAxy", "#eta vs DCA_{xy}", false, 20, -1.0, 1.0, VarManager::kEta, 400, -2.0, 2.0, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "Eta_DCAz", "#eta vs DCA_{z}", false, 20, -1.0, 1.0, VarManager::kEta, 800, -4.0, 4.0, VarManager::kTrackDCAz);
      }
      if (subGroupStr.Contains("dca_fine")) { // Fine binning
        hm->AddHistogram(histClass, "DCAxy_fine", "DCA_{xy}", false, 1000, -0.2, 0.2, VarManager::kTrackDCAxy);
        hm->AddHistogram(histClass, "DCAz_fine", "DCA_{z}", false, 1000, -0.2, 0.2, VarManager::kTrackDCAz);
        hm->AddHistogram(histClass, "IsSPDfirst_Pt_DCAxy_fine", "IsSPDfirst vs p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAxy, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
        hm->AddHistogram(histClass, "IsSPDfirst_Pt_DCAz_fine", "IsSPDfirst vs p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAz, 2, -0.5, 1.5, VarManager::kIsSPDfirst);
        if (subGroupStr.Contains("pt")) {
          hm->AddHistogram(histClass, "PtLow_DCAxy_fine", "p_{T} vs DCA_{xy}", false, 100, 0.0, 2.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAxy);
          hm->AddHistogram(histClass, "PtLow_DCAz_fine", "p_{T} vs DCA_{z}", false, 100, 0.0, 2.0, VarManager::kPt, 1000, -0.2, 0.2, VarManager::kTrackDCAz);
          hm->AddHistogram(histClass, "PtHigh_DCAxy_fine", "p_{T} vs DCA_{xy}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.05, 0.05, VarManager::kTrackDCAxy);
          hm->AddHistogram(histClass, "PtHigh_DCAz_fine", "p_{T} vs DCA_{z}", false, 200, 0.0, 20.0, VarManager::kPt, 1000, -0.05, 0.05, VarManager::kTrackDCAz);
        }
      }
    }
    if (subGroupStr.Contains("muon")) {
      if (!subGroupStr.Contains("ambiguity")) {
        hm->AddHistogram(histClass, "MuonNClusters", "", false, 100, 0.0, 10.0, VarManager::kMuonNClusters);
        hm->AddHistogram(histClass, "pdca", "", false, 200, 0.0, 1000., VarManager::kMuonPDca);
        hm->AddHistogram(histClass, "RAtAbsorberEnd", "", false, 100, 0.0, 200., VarManager::kMuonRAtAbsorberEnd);
        hm->AddHistogram(histClass, "Chi2", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2);
        hm->AddHistogram(histClass, "Chi2MCHMID", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2MatchMCHMID);
        hm->AddHistogram(histClass, "Chi2MCHMFT", "", false, 100, 0.0, 200.0, VarManager::kMuonChi2MatchMCHMFT);
        hm->AddHistogram(histClass, "Chi2MatchScoreMCHMFT", "", false, 100, 0.0, 200.0, VarManager::kMuonMatchScoreMCHMFT);
        hm->AddHistogram(histClass, "MuonCXX", "", false, 100, -1.0, 1.0, VarManager::kMuonCXX);
        hm->AddHistogram(histClass, "MuonCYY", "", false, 100, -1.0, 1.0, VarManager::kMuonCYY);
        hm->AddHistogram(histClass, "MuonCPhiPhi", "", false, 100, -1.0, 1.0, VarManager::kMuonCPhiPhi);
        hm->AddHistogram(histClass, "MuonCTglTgl", "", false, 100, -1.0, 1.0, VarManager::kMuonCTglTgl);
        hm->AddHistogram(histClass, "MuonC1Pt21Pt2", "", false, 100, -1.0, 1.0, VarManager::kMuonC1Pt21Pt2);
        hm->AddHistogram(histClass, "MCHBitMap_vs_pt", "MCH vs pt", false, 1025, 0.0, 1025.0, VarManager::kMCHBitMap, 400, 0, 100, VarManager::kPt);
        hm->AddHistogram(histClass, "MuonTime", "", false, 100, -1.0, 1.0, VarManager::kMuonTime);
        hm->AddHistogram(histClass, "MuonTimeRes", "", false, 100, -1.0, 1.0, VarManager::kMuonTimeRes);
        hm->AddHistogram(histClass, "MuonDcaX_vs_phi", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAx, 200, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
        hm->AddHistogram(histClass, "MuonDcaY_vs_phi", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAy, 200, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
        hm->AddHistogram(histClass, "MuonDcaX_vs_eta", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAx, 500, -5.0, 5.0, VarManager::kEta);
        hm->AddHistogram(histClass, "MuonDcaY_vs_eta", "", false, 2000, -20.0, 20.0, VarManager::kMuonDCAy, 500, -5.0, 5.0, VarManager::kEta);
      } else {
        hm->AddHistogram(histClass, "Pt", "p_{T} distribution", false, 2000, 0.0, 20.0, VarManager::kPt);
        hm->AddHistogram(histClass, "Eta", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta);
        hm->AddHistogram(histClass, "Phi", "#varphi distribution", false, 500, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
        hm->AddHistogram(histClass, "AmbiguityInBunch", "", false, 10, 0.0, 10., VarManager::kMuonNAssocsInBunch);
        hm->AddHistogram(histClass, "AmbiguityOutOfBunch", "", false, 10, 0.0, 10., VarManager::kMuonNAssocsOutOfBunch);
        hm->AddHistogram(histClass, "AmbiguityInBunch_pt", "in bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kMuonNAssocsInBunch);
        hm->AddHistogram(histClass, "AmbiguityOutOfBunch_pt", "out of bunch collision ambiguity vs p_{T}", false, 50, 0.0, 10.0, VarManager::kPt, 10, 0., 10., VarManager::kMuonNAssocsOutOfBunch);
      }
    }

    if (subGroupStr.Contains("muon-pdca")) {
      hm->AddHistogram(histClass, "p", "p", false, 200, 0.0, 20.0, VarManager::kP);
      hm->AddHistogram(histClass, "pdca_vs_p", "pDCA vs p", false, 2000, 0.0, 20.0, VarManager::kP, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_pt", "pDCA vs pt", false, 2000, 0.0, 20.0, VarManager::kPt, 200, 0.0, 1000., VarManager::kMuonPDca);
      hm->AddHistogram(histClass, "pdca_vs_Rabs", "pDCA vs R_{abs}", false, 100, 0., 200., VarManager::kMuonRAtAbsorberEnd, 200, 0.0, 1000., VarManager::kMuonPDca);
    }
    if (subGroupStr.Contains("mft-pid")) {
      hm->AddHistogram(histClass, "hMftTrackEtaVsPt", "", false, 100, -5.f, -2.f, VarManager::kEta, 100, 0.f, 20.f, VarManager::kPt);
      hm->AddHistogram(histClass, "hMftNClusters", "", false, 16, -0.5f, 15.5f, VarManager::kMftNClusters);
      hm->AddHistogram(histClass, "hMftClusterSize", "", false, 101, -0.5f, 100.5f, VarManager::kMftClusterSize);
      hm->AddHistogram(histClass, "hMftMeanClusterSize", "", false, 200, 0.f, 20.f, VarManager::kMftMeanClusterSize);
    }
    if (subGroupStr.Contains("mc")) {
      hm->AddHistogram(histClass, "Pt_vs_PtMC", "pT vs MC pT", false, 200, 0.0, 20.0, VarManager::kPt, 200, 0.0, 20.0, VarManager::kMCPt);
      hm->AddHistogram(histClass, "Eta_vs_EtaMC", "#eta vs MC #eta", false, 50, -1.0, 1.0, VarManager::kEta, 50, -1.0, 1.0, VarManager::kMCEta);
      hm->AddHistogram(histClass, "Phi_vs_PhiMC", "#varphi vs MC #varphi", false, 50, 0.0, 2. * o2::constants::math::PI, VarManager::kPhi, 50, 0.0, 2. * o2::constants::math::PI, VarManager::kMCPhi);
      hm->AddHistogram(histClass, "TrackPDGcode", "PDG code of track", false, 10001, -5000, 5000, VarManager::kMCPdgCode);
    }
    if (subGroupStr.Contains("mcmother")) {
      hm->AddHistogram(histClass, "MotherPDGcode", "PDG code of mother", false, 10001, -5000, 5000, VarManager::kMCMotherPdgCode);
    }
    if (subGroupStr.Contains("dmeson")) {
      hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
    }
    if (subGroupStr.Contains("tpcpidvstofpid")) {
      hm->AddHistogram(histClass, "tpcNSigmaKa_tofNSigmaKa", "", false, 200, -10., 10., VarManager::kTPCnSigmaKa, 200, -10., 10., VarManager::kTOFnSigmaKa);
      hm->AddHistogram(histClass, "tpcNSigmaPi_tofNSigmaPi", "", false, 200, -10., 10., VarManager::kTPCnSigmaPi, 200, -10., 10., VarManager::kTOFnSigmaPi);
    }
  }

  if (!groupStr.CompareTo("mctruth_triple")) {
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kPairEta, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt1);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt2);
    hm->AddHistogram(histClass, "Eta_Pt_Photon", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kPairEta, 200, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPairPhi);
    hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 4500, 0.0, 4.5, VarManager::kPairMassDau);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Photon", "", false, 500, 0.0, 5.0, VarManager::kPairMassDau, 500, 0.0, 5.0, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 2000, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Pt_Photon", "", false, 500, 0.0, 5.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Mass_DileptonPhoton", "", false, 4500, 0.0, 4.5, VarManager::kPairMass);
    hm->AddHistogram(histClass, "Pt_DileptonPhoton", "", false, 2000, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Mass_Pt_DileptonPhoton", "", false, 500, 0.0, 5.0, VarManager::kPairMass, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass", "", false, 1500, 0.0, 1.5, VarManager::kDeltaMass);
    hm->AddHistogram(histClass, "DeltaMass_ptdileptonphoton", "", false, 1000, 0.0, 1.0, VarManager::kDeltaMass, 3000, 0.0, 30.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass_Jpsi", "", false, 1500, 3, 4.5, (VarManager::kDeltaMass_jpsi));
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
  }
  if (!groupStr.CompareTo("mctruth_pair")) {
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Pt", "", false, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 200, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt1);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt2);
    hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 40, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi_Eta", "#phi vs #eta distribution", false, 200, -5.0, 5.0, VarManager::kEta, 200, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
  }
  if (!groupStr.CompareTo("mctruth_quad")) {
    hm->AddHistogram(histClass, "hMass_defaultDileptonMass", "", false, 1000, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass);
    hm->AddHistogram(histClass, "hPt", "", false, 150, 0.0, 15.0, VarManager::kQuadPt);
    hm->AddHistogram(histClass, "hMass_defaultDileptonMass_Pt", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 150, 0.0, 15.0, VarManager::kQuadPt);
    hm->AddHistogram(histClass, "hQ", "", false, 150, 0.0, 3.0, VarManager::kQ);
    hm->AddHistogram(histClass, "hDeltaR1", "", false, 100, 0.0, 10.0, VarManager::kDeltaR1);
    hm->AddHistogram(histClass, "hDeltaR2", "", false, 100, 0.0, 10.0, VarManager::kDeltaR2);
    hm->AddHistogram(histClass, "hDeltaR", "", false, 100, 0.0, 10.0, VarManager::kDeltaR);
    hm->AddHistogram(histClass, "hDiTrackMass", "", false, 300, 0.0, 3.0, VarManager::kDitrackMass);
    hm->AddHistogram(histClass, "hMCPt_MCRap", "", false, 200, 0.0, 20.0, VarManager::kMCPt, 100, -2.0, 2.0, VarManager::kMCY);
    hm->AddHistogram(histClass, "hMCPhi", "", false, 100, -TMath::Pi(), TMath::Pi(), VarManager::kMCPhi);
  }
  if (!groupStr.CompareTo("mctruth_track")) {
    hm->AddHistogram(histClass, "PtMC", "MC pT", false, 200, 0.0, 20.0, VarManager::kMCPt);
    hm->AddHistogram(histClass, "EtaMC", "MC #eta", false, 50, -5.0, 5.0, VarManager::kMCEta);
    hm->AddHistogram(histClass, "PhiMC", "MC #phi", false, 50, -6.3, 6.3, VarManager::kMCPhi);
    hm->AddHistogram(histClass, "YMC", "MC y", false, 50, -5.0, 5.0, VarManager::kMCY);
    hm->AddHistogram(histClass, "CentFT0CMC", "MC Cent. FT0C", false, 18, 0., 90., VarManager::kCentFT0C);
    hm->AddHistogram(histClass, "PtMC_YMC", "MC pT vs MC y", false, 120, 0.0, 30.0, VarManager::kMCPt, 1000, -5.0, 5.0, VarManager::kMCY);
    hm->AddHistogram(histClass, "EtaMC_PtMC", "", false, 40, -2.0, 2.0, VarManager::kMCEta, 200, 0.0, 20.0, VarManager::kMCPt);
    hm->AddHistogram(histClass, "VzMC", "MC vz", false, 100, -15.0, 15.0, VarManager::kMCVz);
    hm->AddHistogram(histClass, "VzMC_VtxZMC", "MC vz vs MC vtxZ", false, 50, -15.0, 15.0, VarManager::kMCVz, 50, -15.0, 15.0, VarManager::kMCVtxZ);
    hm->AddHistogram(histClass, "Weight", "", false, 50, 0.0, 5.0, VarManager::kMCParticleWeight);
  }

  if (!groupStr.CompareTo("pair")) {
    if (subGroupStr.Contains("barrel")) {
      hm->AddHistogram(histClass, "Mass", "", false, 500, 0.0, 5.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Mass_HighRange", "", false, 375, 0.0, 15.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Pt", "", false, 2000, 0.0, 20., VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 125, 0.0, 5.0, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
      double massBins[76];
      for (int i = 0; i < 76; i++) {
        massBins[i] = 1.5 + i * 0.04;
      }
      double ptBins[70];
      for (int i = 0; i <= 50; i++) {
        ptBins[i] = i * 0.01;
      }
      for (int i = 1; i <= 19; i++) {
        ptBins[50 + i] = 0.5 + i * 0.5;
      }
      hm->AddHistogram(histClass, "Mass_PtFine", "", false, 75, massBins, VarManager::kMass, 69, ptBins, VarManager::kPt);
      hm->AddHistogram(histClass, "Eta_Pt", "", false, 40, -2.0, 2.0, VarManager::kEta, 40, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Y_Pt", "", false, 40, -2.0, 2.0, VarManager::kRap, 40, 0.0, 20.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 500, 0.0, 5.0, VarManager::kMass);
      if (subGroupStr.Contains("pbpb")) {
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 125, 0.0, 5.0, VarManager::kMass, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Pt_CentFT0C", "", false, 100, 0.0, 10.0, VarManager::kPt, 20, 0.0, 100.0, VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C", "", false, 75, 1.5, 4.5, VarManager::kMass, 20, 0.0, 10.0, VarManager::kPt, 10, 0.0, 100.0, VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("mult")) {
        hm->AddHistogram(histClass, "Mass_Pt_MultFV0A", "", false, 200, 0.0, 5.0, VarManager::kMass, 40, 0.0, 40.0, VarManager::kPt, 100, 0.0, 25000.0, VarManager::kMultFV0A);
        hm->AddHistogram(histClass, "Mass_VtxNcontribReal", "Mass vs VtxNcontribReal", false, 200, 0.0, 5.0, VarManager::kMass, 200, 0, 200.0, VarManager::kVtxNcontribReal);
        hm->AddHistogram(histClass, "MultITSWithPV", "MultITSWithPV", false, 200, 0, 200.0, VarManager::kMultNTracksHasITS);
        hm->AddHistogram(histClass, "MultTPCWithPV", "MultTPCWithPV", false, 200, 0, 200.0, VarManager::kMultNTracksHasTPC);
        hm->AddHistogram(histClass, "MultITSTPCWithPV", "MultITSTPCWithPV", false, 200, 0, 200.0, VarManager::kMultNTracksITSTPC);
        hm->AddHistogram(histClass, "Mass_MultITSWithPV", "Mass vs MultITSWithPV", false, 200, 0.0, 5.0, VarManager::kMass, 200, 0, 200.0, VarManager::kMultNTracksHasITS);
        hm->AddHistogram(histClass, "Mass_MultTPCWithPV", "Mass vs MultTPCWithPV", false, 200, 0.0, 5.0, VarManager::kMass, 200, 0, 200.0, VarManager::kMultNTracksHasTPC);
        hm->AddHistogram(histClass, "Mass_MultITSTPCWithPV", "Mass vs MultITSTPCWithPV", false, 200, 0.0, 5.0, VarManager::kMass, 200, 0, 200.0, VarManager::kMultNTracksITSTPC);
      }
      if (subGroupStr.Contains("polarization")) {
        hm->AddHistogram(histClass, "cosThetaHE", "", false, 100, -1., 1., VarManager::kCosThetaHE);
        hm->AddHistogram(histClass, "cosThetaCS", "", false, 100, -1., 1., VarManager::kCosThetaCS);
        hm->AddHistogram(histClass, "PhiHE", "", false, 100, -o2::constants::math::PI, o2::constants::math::PI, VarManager::kPhiHE);
        hm->AddHistogram(histClass, "PhiCS", "", false, 100, -o2::constants::math::PI, o2::constants::math::PI, VarManager::kPhiCS);
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaHE", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -1., 1., VarManager::kCosThetaHE);
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaCS", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -1., 1., VarManager::kCosThetaCS);
        hm->AddHistogram(histClass, "Mass_Pt_PhiHE", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -o2::constants::math::PI, o2::constants::math::PI, VarManager::kPhiHE);
        hm->AddHistogram(histClass, "Mass_Pt_PhiCS", "", false, 100, 1.0, 5.0, VarManager::kMass, 250, 0.0, 25.0, VarManager::kPt, 40, -o2::constants::math::PI, o2::constants::math::PI, VarManager::kPhiCS);
      }
      if (subGroupStr.Contains("upsilon")) {
        hm->AddHistogram(histClass, "MassUpsilon_Pt", "", false, 500, 7.0, 12.0, VarManager::kMass, 400, 0.0, 40.0, VarManager::kPt);
      }
      if (subGroupStr.Contains("dalitz")) {
        hm->AddHistogram(histClass, "MassLow", "", false, 500, 0.0, 0.05, VarManager::kMass);
        hm->AddHistogram(histClass, "PsiPair", "", false, 200, -1.5, 1.5, VarManager::kPsiPair);
        hm->AddHistogram(histClass, "PsiPair_DeltaPhi", "", false, 100, -0.5, 0.5, VarManager::kDeltaPhiPair, 100, -1.5, 1.5, VarManager::kPsiPair);
      }
      if (subGroupStr.Contains("vertexing")) {
        hm->AddHistogram(histClass, "UsedKF", "", false, 2, -0.5, 1.5, VarManager::kUsedKF);
        hm->AddHistogram(histClass, "Lz", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLz);
        hm->AddHistogram(histClass, "Lxy", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxy);
        hm->AddHistogram(histClass, "Lxyz", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyz);
        hm->AddHistogram(histClass, "Tauz", "", false, 1000, -0.02, 0.02, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "Tauxy", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "Tauxy_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "Tauz_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "LzProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLzProjected);
        hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -1.0, 1.0, VarManager::kVertexingLxyzProjected);
        hm->AddHistogram(histClass, "TauzProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauzProjected);
        hm->AddHistogram(histClass, "TauxyProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "TauxyProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "TauzProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -0.03, 0.03, VarManager::kVertexingTauzProjected);
        hm->AddHistogram(histClass, "TauxyzProj", "", false, 1000, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
        hm->AddHistogram(histClass, "LxyProj_Pt", "", false, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LzProj_Mass_Pt", "", false, 50, 2.0, 4.0, VarManager::kMass, 10, 0.0, 20.0, VarManager::kPt, 1000, -1.0, 1.0, VarManager::kVertexingLzProjected);
        hm->AddHistogram(histClass, "CosPointingAngle", "", false, 200, -1.0, 1.0, VarManager::kCosPointingAngle);
      }

      if (subGroupStr.Contains("kalman-filter")) {
        hm->AddHistogram(histClass, "LxyErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyErr);
        hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
        hm->AddHistogram(histClass, "TauxyErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauxyErr);
        hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
        hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
        hm->AddHistogram(histClass, "KFMass", "", false, 500, 0.0, 5.0, VarManager::kKFMass);
        hm->AddHistogram(histClass, "LxyOverDLxy", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyOverErr);
        hm->AddHistogram(histClass, "LzOverDLz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzOverErr);
        hm->AddHistogram(histClass, "LxyzOverDLxyz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzOverErr);
        hm->AddHistogram(histClass, "KFTauxy", "", false, 1000, -0.2, 0.2, VarManager::kVertexingTauxy);
        hm->AddHistogram(histClass, "KFTrack0DCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFTrack0DCAxyz);
        hm->AddHistogram(histClass, "KFTrack1DCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFTrack1DCAxyz);
        hm->AddHistogram(histClass, "KFTracksDCAxyzMax", "", false, 400, -2.0, 2.0, VarManager::kKFTracksDCAxyzMax);
        hm->AddHistogram(histClass, "KFDCAxyzBetweenProngs", "", false, 400, -2.0, 2.0, VarManager::kKFDCAxyzBetweenProngs);
        hm->AddHistogram(histClass, "KFTrack0DCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFTrack0DCAxy);
        hm->AddHistogram(histClass, "KFTrack1DCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFTrack1DCAxy);
        hm->AddHistogram(histClass, "KFTracksDCAxyMax", "", false, 400, -2.0, 2.0, VarManager::kKFTracksDCAxyMax);
        hm->AddHistogram(histClass, "KFDCAxyBetweenProngs", "", false, 400, -2.0, 2.0, VarManager::kKFDCAxyBetweenProngs);
        hm->AddHistogram(histClass, "KFChi2OverNDFGeo", "", false, 150, -5, 10, VarManager::kKFChi2OverNDFGeo);
        hm->AddHistogram(histClass, "KFTrack0DeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack0DeviationFromPV);
        hm->AddHistogram(histClass, "KFTrack1DeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack1DeviationFromPV);
        hm->AddHistogram(histClass, "KFTrack0DeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack0DeviationxyFromPV);
        hm->AddHistogram(histClass, "KFTrack1DeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFTrack1DeviationxyFromPV);
        hm->AddHistogram(histClass, "KFPairDCAxyz", "", false, 400, -2.0, 2.0, VarManager::kKFJpsiDCAxyz);
        hm->AddHistogram(histClass, "KFPairDCAxy", "", false, 400, -2.0, 2.0, VarManager::kKFJpsiDCAxy);
        hm->AddHistogram(histClass, "KFPairDeviationFromPV", "", false, 150, 0, 15e+6, VarManager::kKFPairDeviationFromPV);
        hm->AddHistogram(histClass, "KFPairDeviationxyFromPV", "", false, 150, 0, 15e+6, VarManager::kKFPairDeviationxyFromPV);
        hm->AddHistogram(histClass, "KFCosPA", "", false, 300, -1.5, 1.5, VarManager::kKFCosPA);
        hm->AddHistogram(histClass, "KFMassGeoTop", "", false, 500, 0.0, 5.0, VarManager::kKFMassGeoTop);
        hm->AddHistogram(histClass, "KFChi2OverNDFGeoTop", "", false, 150, -5, 10, VarManager::kKFChi2OverNDFGeoTop);
        hm->AddHistogram(histClass, "KFNTrks2PV", "", false, 210, -10, 200, VarManager::kKFNContributorsPV);
        hm->AddHistogram(histClass, "Mass_DCAxyzTwoProngs", "", false, 500, 0.0, 5.0, VarManager::kMass, 400, -2.0, 2.0, VarManager::kKFDCAxyzBetweenProngs);
        hm->AddHistogram(histClass, "Mass_DCAxyTwoProngs", "", false, 500, 0.0, 5.0, VarManager::kMass, 400, -2.0, 2.0, VarManager::kKFDCAxyBetweenProngs);
        hm->AddHistogram(histClass, "Mass_KFChi2OverNDFGeo", "", false, 500, 0.0, 5.0, VarManager::kMass, 150, -5, 10, VarManager::kKFChi2OverNDFGeo);
      }
      if (subGroupStr.Contains("run2-vertexing-definitions")) {
        hm->AddHistogram(histClass, "LzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLzProjected);
        hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyProjected);
        hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyzProjected);
        hm->AddHistogram(histClass, "TauzProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauzProjected);
        hm->AddHistogram(histClass, "TauxyProjNs", "", false, 4000, -0.1, 0.1, VarManager::kVertexingTauxyProjectedNs);
        hm->AddHistogram(histClass, "TauxyProjCm", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauxyProjected);
      }
      if (subGroupStr.Contains("multidimentional-vertexing-histograms")) {
        hm->AddHistogram(histClass, "pT_TauxyProj", "", false, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected, 20, 0.0, 20., VarManager::kPt);
        hm->AddHistogram(histClass, "InvMass_TauxyProj", "", false, 500, 0.0, 5.0, VarManager::kMass, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "Eta_TauxyProj", "", false, 40, -2.0, 2.0, VarManager::kEta, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "Rap_TauxyProj", "", false, 200, -1.0, 1.0, VarManager::kRap, 1000, -0.2, 0.2, VarManager::kVertexingTauxyProjected);

        const int kNvarsPair = 4;
        const int kInvMassNbins = 3;
        double InvMassBinLims[kInvMassNbins + 1] = {2.2, 2.6, 3.4, 3.6};

        const int kPtNbins = 10;
        double PtBinLims[kPtNbins + 1] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 12., 20.};

        const int kTauNBins = 500;
        double TauBinLims[kTauNBins + 1];
        for (int i = 0; i <= kTauNBins; ++i)
          TauBinLims[i] = -0.3 + (0.0015 * i);

        TArrayD nCutsBinLimits[kNvarsPair];
        nCutsBinLimits[0] = TArrayD(kInvMassNbins + 1, InvMassBinLims);
        nCutsBinLimits[1] = TArrayD(kPtNbins + 1, PtBinLims);
        nCutsBinLimits[2] = TArrayD(kTauNBins + 1, TauBinLims);
        nCutsBinLimits[3] = TArrayD(kTauNBins + 1, TauBinLims);

        int varsPair[kNvarsPair] = {VarManager::kMass, VarManager::kPt, VarManager::kVertexingTauzProjected, VarManager::kVertexingTauxyProjected};
        hm->AddHistogram(histClass, "tau_MultiD", "Invariant mass vs. pT vs. eta vs. rapidity vs. Run2 tau", kNvarsPair, varsPair, nCutsBinLimits);
      }

      if (subGroupStr.Contains("flow")) {
        hm->AddHistogram(histClass, "Mass_u2q2", "u_{2}Q_{2}^{A} vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kU2Q2);
        hm->AddHistogram(histClass, "Mass_u3q3", "u_{3}Q_{3}^{A} vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kU3Q3);
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhi", "cos 2(#varphi-#Psi_{2}^{A}) vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhi);
        hm->AddHistogram(histClass, "Mass_cos3DeltaPhi", "cos 3(#varphi-#Psi_{3}^{A}) vs m", true, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos3DeltaPhi);
      }
    } else if (subGroupStr.Contains("dimuon")) {
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 750, 0.0, 15.0, VarManager::kMass, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Rapidity", "", false, 750, 0.0, 15.0, VarManager::kMass, 150, 2.5, 4.0, VarManager::kRap);
      hm->AddHistogram(histClass, "Mass_Phi", "", false, 750, 0.0, 15.0, VarManager::kMass, 180, constants::math::PI, 2 * constants::math::PI, VarManager::kPhi);
      if (subGroupStr.Contains("dimuon-multi-diff")) {
        int varsKine[3] = {VarManager::kMass, VarManager::kPt, VarManager::kRap};
        int binsKine[3] = {250, 120, 60};
        double xminKine[3] = {0.0, 0.0, 2.5};
        double xmaxKine[3] = {5.0, 30.0, 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Rapidity", "", 3, varsKine, binsKine, xminKine, xmaxKine, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-high-mass-multi-diff")) {
        int varsKine[3] = {VarManager::kMass, VarManager::kPt, VarManager::kRap};
        int binsKine[3] = {250, 120, 60};
        double xminKine[3] = {7.0, 0.0, 2.5};
        double xmaxKine[3] = {12.0, 30.0, 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Rapidity", "", 3, varsKine, binsKine, xminKine, xmaxKine, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-centr")) {
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("qc")) {
        hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 750, 0.0, 15.0, VarManager::kMass);
        hm->AddHistogram(histClass, "DeltaPtotTracks", "", false, 2000, -100., 100., VarManager::kDeltaPtotTracks);
        hm->AddHistogram(histClass, "Mass_DeltaPtotTracks", "", false, 150, 2.0, 5.0, VarManager::kMass, 200, -100., 100., VarManager::kDeltaPtotTracks);
      }
      if (subGroupStr.Contains("mixedevent")) {
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhiMu1", "cos 2(#varphi_{#mu1}-#phi_{#mu#mu}) vs m", false, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhiMu1);
        hm->AddHistogram(histClass, "Mass_cos2DeltaPhiMu2", "cos 2(#varphi_{#mu2}-#phi_{#mu#mu}) vs m", false, 125, 0.0, 5.0, VarManager::kMass, 100, -1.0, 1.0, VarManager::kCos2DeltaPhiMu2);
        hm->AddHistogram(histClass, "R2SP1_CentFT0C", "mass vs centrality vs. R2SP_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2SP1);
        hm->AddHistogram(histClass, "R2SP2_CentFT0C", "mass vs centrality vs. R2SP_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2SP2);
        hm->AddHistogram(histClass, "R2EP1_CentFT0C", "mass vs centrality vs. R2EP_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2EP1);
        hm->AddHistogram(histClass, "R2EP2_CentFT0C", "mass vs centrality vs. R2EP_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kTwoR2EP2);
        hm->AddHistogram(histClass, "U2Q2_CentFT0C_ev1", "mass vs. centrality vs. U2Q2_event1", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kU2Q2Ev1);
        hm->AddHistogram(histClass, "U2Q2_CentFT0C_ev2", "mass vs. centrality vs. U2Q2_event2", false, 125, 0.0, 5.0, VarManager::kMass, 9, 0.0, 90.0, VarManager::kCentFT0C, 100, -10.0, 10.0, VarManager::kU2Q2Ev2);
      }
      if (subGroupStr.Contains("metest")) {
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V2ME_SP", "Mass_Pt_CentFT0C_V2ME_SP", true, 250, 0.0, 5.0, VarManager::kMass, 200, 0.0, 20.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV2ME_SP, VarManager::kWV2ME_SP);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V2ME_EP", "Mass_Pt_CentFT0C_V2ME_EP", true, 250, 0.0, 5.0, VarManager::kMass, 200, 0.0, 20.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV2ME_EP, VarManager::kWV2ME_EP);
      }
      if (subGroupStr.Contains("cumulantme")) {
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M11REFoverMpME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM11REFoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M1111REFoverMpME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM1111REFoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M01POIoverMpME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM01POIoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M0111POIoverMpME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM0111POIoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2REFbydimuonsME, VarManager::kM11REFoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4REFME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR4REFbydimuonsME, VarManager::kM1111REFoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2POIME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2POIME, VarManager::kM01POIoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4POIME", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 9, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR4POIME, VarManager::kM0111POIoverMpME);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V22ME", "Mass_Pt_CentFT0C_V22ME", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV22ME, VarManager::kWV22ME);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V24ME", "Mass_Pt_CentFT0C_V24ME", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV24ME, VarManager::kWV24ME);
      }
      if (subGroupStr.Contains("dimuon-polarization-he")) {
        int varspTHE[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsrapHE[4] = {VarManager::kMass, VarManager::kRap, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binspT[4] = {100, 20, 20, 20};
        int binsy[4] = {100, 10, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {5., 20., 1., +3.14};
        double xminy[4] = {1., 2.5, -1., -3.14};
        double xmaxy[4] = {5., 4.0, 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaHE_phiHE", "", 4, varspTHE, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_y_cosThetaHE_phiHE", "", 4, varsrapHE, binsy, xminy, xmaxy, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-cs")) {
        int varspTCS[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int varsrapCS[4] = {VarManager::kMass, VarManager::kRap, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binspT[4] = {100, 20, 20, 20};
        int binsy[4] = {100, 10, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {5., 20., 1., +3.14};
        double xminy[4] = {1., 2.5, -1., -3.14};
        double xmaxy[4] = {5., 4.0, 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaCS_phiCS", "", 4, varspTCS, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_y_cosThetaCS_phiCS", "", 4, varsrapCS, binsy, xminy, xmaxy, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("upsilon-polarization-he")) {
        int varspTHE[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsrapHE[4] = {VarManager::kMass, VarManager::kRap, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binspT[4] = {100, 20, 20, 20};
        int binsy[4] = {100, 10, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {15., 20., 1., +3.14};
        double xminy[4] = {1., 2.5, -1., -3.14};
        double xmaxy[4] = {15., 4.0, 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaHE_phiHE", "", 4, varspTHE, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_y_cosThetaHE_phiHE", "", 4, varsrapHE, binsy, xminy, xmaxy, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("upsilon-polarization-cs")) {
        int varspTCS[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int varsrapCS[4] = {VarManager::kMass, VarManager::kRap, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binspT[4] = {100, 20, 20, 20};
        int binsy[4] = {100, 10, 20, 20};
        double xminpT[4] = {1., 0., -1., -3.14};
        double xmaxpT[4] = {15., 20., 1., +3.14};
        double xminy[4] = {1., 2.5, -1., -3.14};
        double xmaxy[4] = {15., 4.0, 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_cosThetaCS_phiCS", "", 4, varspTCS, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_y_cosThetaCS_phiCS", "", 4, varsrapCS, binsy, xminy, xmaxy, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-vp")) {
        int varspTVP[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCosPhiVP, VarManager::kPhiVP};
        int varsrapVP[4] = {VarManager::kMass, VarManager::kRap, VarManager::kCosPhiVP, VarManager::kPhiVP};
        int binspT[4] = {100, 20, 24, 24};
        int binsy[4] = {100, 10, 24, 24};
        double xminpT[4] = {1., 0., -1., 0.};
        double xmaxpT[4] = {5., 20., 1., +3.14};
        double xminy[4] = {1., 2.5, -1., 0.};
        double xmaxy[4] = {5., 4.0, 1., +3.14};
        hm->AddHistogram(histClass, "Mass_Pt_phiVP", "", 4, varspTVP, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_y_phiVP", "", 4, varsrapVP, binsy, xminy, xmaxy, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-rap")) {
        int vars[4] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kRap};
        int binspT[4] = {300, 200, 10, 6};
        double xminpT[4] = {2., 0., 0, 2.5};
        double xmaxpT[4] = {8., 20., 100, 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_Rap", "", 4, vars, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-he-pbpb")) {
        int varsHEpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binspT[5] = {150, 30, 10, 10, 10};
        double xminpT[5] = {2., 0., 0, -1., -3.14};
        double xmaxpT[5] = {5., 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaHE", "", 5, varsHEpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-lowmass-he-pbpb")) {
        int varsHEpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binspT[5] = {200, 30, 10, 10, 10};
        double xminpT[5] = {0.2, 0., 0, -1., -3.14};
        double xmaxpT[5] = {1.2, 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaHE_lowmass", "", 5, varsHEpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-cs-pbpb")) {
        int varsCSpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binspT[5] = {150, 30, 10, 10, 10};
        double xminpT[5] = {2., 0., 0, -1., -3.14};
        double xmaxpT[5] = {5., 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaCS", "", 5, varsCSpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-lowmass-cs-pbpb")) {
        int varsCSpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binspT[5] = {200, 30, 10, 10, 10};
        double xminpT[5] = {0.2, 0., 0, -1., -3.14};
        double xmaxpT[5] = {1.2, 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaCS_lowmass", "", 5, varsCSpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-vp-pbpb")) {
        int varsVPpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosPhiVP, VarManager::kPhiVP};
        int binspT[5] = {150, 30, 10, 24, 24};
        double xminpT[5] = {2., 0., 0, -1., 0.};
        double xmaxpT[5] = {5., 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_phiVP", "", 5, varsVPpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-polarization-lowmass-vp-pbpb")) {
        int varsVPpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosPhiVP, VarManager::kPhiVP};
        int binspT[5] = {200, 30, 10, 24, 24};
        double xminpT[5] = {0.2, 0., 0, -1., 0.};
        double xmaxpT[5] = {1.2, 3., 100, 1., 3.14};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_phiVP_lowmass", "", 5, varsVPpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-rap-polarization-he-pbpb")) {
        int varsHEpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaHE, VarManager::kRap};
        int binspT[5] = {150, 30, 10, 10, 6};
        double xminpT[5] = {2., 0., 0, -1., 2.5};
        double xmaxpT[5] = {5., 3., 100, 1., 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaHE_Rap", "", 5, varsHEpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-rap-polarization-cs-pbpb")) {
        int varsCSpbpb[5] = {VarManager::kMass, VarManager::kPt, VarManager::kCentFT0C, VarManager::kCosThetaCS, VarManager::kRap};
        int binspT[5] = {150, 30, 10, 10, 6};
        double xminpT[5] = {2., 0., 0, -1., 2.5};
        double xmaxpT[5] = {5., 3., 100, 1., 4.0};
        hm->AddHistogram(histClass, "Mass_Pt_Cent_cosThetaCS_Rap", "", 5, varsCSpbpb, binspT, xminpT, xmaxpT, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-midmult-polarization-he")) {
        int varsITSTPCMulHE[4] = {VarManager::kMass, VarManager::kMultNTracksITSTPC, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsITSMulHE[4] = {VarManager::kMass, VarManager::kMultNTracksHasITS, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 120., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_ITSTPCMult_cosThetaHE_phiHE", "", 4, varsITSTPCMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_ITSMult_cosThetaHE_phiHE", "", 4, varsITSMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-midmult-polarization-cs")) {
        int varsITSTPCMulCS[4] = {VarManager::kMass, VarManager::kMultNTracksITSTPC, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int varsITSMulCS[4] = {VarManager::kMass, VarManager::kMultNTracksHasITS, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 120., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_ITSTPCMult_cosThetaCS_phiCS", "", 4, varsITSTPCMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_ITSMult_cosThetaCS_phiCS", "", 4, varsITSMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-fwdmult-polarization-he")) {
        int varsFT0AMulHE[4] = {VarManager::kMass, VarManager::kMultFT0A, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int varsFV0AMulHE[4] = {VarManager::kMass, VarManager::kMultFV0A, VarManager::kCosThetaHE, VarManager::kPhiHE};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 3000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_FT0AMult_cosThetaHE_phiHE", "", 4, varsFT0AMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_FV0AMult_cosThetaHE_phiHE", "", 4, varsFV0AMulHE, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("dimuon-fwdmult-polarization-cs")) {
        int varsFT0AMulCS[4] = {VarManager::kMass, VarManager::kMultFT0A, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int varsFV0AMulCS[4] = {VarManager::kMass, VarManager::kMultFV0A, VarManager::kCosThetaCS, VarManager::kPhiCS};
        int binsMul[4] = {100, 20, 20, 20};
        double xminMul[4] = {1., 0., -1., -3.14};
        double xmaxMul[4] = {5., 3000., 1., +3.14};
        hm->AddHistogram(histClass, "Mass_FT0AMult_cosThetaCS_phiCS", "", 4, varsFT0AMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
        hm->AddHistogram(histClass, "Mass_FV0AMult_cosThetaCS_phiCS", "", 4, varsFV0AMulCS, binsMul, xminMul, xmaxMul, 0, -1, kFALSE);
      }
      if (subGroupStr.Contains("vertexing-forward")) {
        hm->AddHistogram(histClass, "Lxyz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyz);
        hm->AddHistogram(histClass, "Lz", "", false, 100, 0.0, 10.0, VarManager::kVertexingLz);
        hm->AddHistogram(histClass, "Tauz", "", false, 100, -0.01, 0.01, VarManager::kVertexingTauz);
        hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
        hm->AddHistogram(histClass, "LzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzErr);
        hm->AddHistogram(histClass, "TauzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauzErr);
        hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
        hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
        hm->AddHistogram(histClass, "Pz", "", false, 200, 0.0, 100.0, VarManager::kVertexingPz);
        hm->AddHistogram(histClass, "Secondary_Vertexing_z", "", false, 200, -10.0, 10.0, VarManager::kVertexingSV);
        hm->AddHistogram(histClass, "Primary_Vertexing_z", "", false, 200, -10.0, 10.0, VarManager::kVtxZ);
      }
      if (subGroupStr.Contains("pbpb")) {
        hm->AddHistogram(histClass, "Mass_CentVZERO", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Pt_CentVZERO", "", false, 120, 0.0, 30.0, VarManager::kPt, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Rapidity_CentVZERO", "", false, 200, 2.5, 4.0, VarManager::kRap, 100, 0., 100., VarManager::kCentVZERO);
        hm->AddHistogram(histClass, "Mass_CentFT0C", "", false, 750, 0.0, 15.0, VarManager::kMass, 100, 0., 100., VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Pt_CentFT0C", "", false, 120, 0.0, 30.0, VarManager::kPt, 100, 0., 100., VarManager::kCentFT0C);
        hm->AddHistogram(histClass, "Rapidity_CentFT0C", "", false, 200, 2.5, 4.0, VarManager::kRap, 100, 0., 100., VarManager::kCentFT0C);
      }
      if (subGroupStr.Contains("lowmass")) {
        hm->AddHistogram(histClass, "MassLow", "", false, 400, 0.0, 2.0, VarManager::kMass);
      }
      if (subGroupStr.Contains("lmmumu")) {
        hm->AddHistogram(histClass, "Mass_QuadDCAabsXY", "", false, 250, 0.0, 5.0, VarManager::kMass, 900, 0.0, 3, VarManager::kQuadDCAabsXY);
        hm->AddHistogram(histClass, "Mass_Lxyz", "", false, 250, 0.0, 5.0, VarManager::kMass, 1000, 0.0, 5, VarManager::kVertexingLxyz);
        hm->AddHistogram(histClass, "Mass_OpeningAngle", "", false, 250, 0.0, 5.0, VarManager::kMass, 800, 0, 0.8, VarManager::kOpeningAngle);
      }
      if (subGroupStr.Contains("flow-dimuon-high-mass")) {
        int varV2[6] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C, VarManager::kU2Q2, VarManager::kCos2DeltaPhi};

        int bins[6] = {50, 30, 6, 18, 200, 40};
        double minBins[6] = {7.0, 0.0, 2.5, 0.0, -10.0, -2.0};
        double maxBins[6] = {12.0, 30.0, 4.0, 90.0, 10.0, 2.0};
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_V2", "", 6, varV2, bins, minBins, maxBins, 0, -1, kTRUE);
      }
      if (subGroupStr.Contains("flow-dimuon")) {
        int varV2[6] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C, VarManager::kU2Q2, VarManager::kCos2DeltaPhi};
        // int varV3[6] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C, VarManager::kU3Q3, VarManager::kCos3DeltaPhi}; // removed temporarily

        int bins[6] = {250, 60, 6, 18, 200, 40};
        double minBins[6] = {0.0, 0.0, 2.5, 0.0, -10.0, -2.0};
        double maxBins[6] = {5.0, 30.0, 4.0, 90.0, 10.0, 2.0};
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_V2", "", 6, varV2, bins, minBins, maxBins, 0, -1, kTRUE);
        // hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_V3", "", 6, varV3, bins, minBins, maxBins, 0, -1, kTRUE); // removed temporarily
      }
      if (subGroupStr.Contains("flow-ccdb")) {
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V2SPwR", "Mass_Pt_CentFT0C_V2SPwR", true, 250, 0.0, 5.0, VarManager::kMass, 200, 0.0, 20.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV2SP, VarManager::kWV2SP);
        hm->AddHistogram(histClass, "Mass_Pt_CentFT0C_V2EPwR", "Mass_Pt_CentFT0C_V2EPwR", true, 250, 0.0, 5.0, VarManager::kMass, 200, 0.0, 20.0, VarManager::kPt, 90, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kV2EP, VarManager::kWV2EP);
      }
      if (subGroupStr.Contains("cumulant")) {
        int var[4] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C};
        int bins[4] = {250, 60, 6, 18};
        double minBins[4] = {0.0, 0.0, 2.5, 0.0};
        double maxBins[4] = {5.0, 30.0, 4.0, 90.0};
        hm->AddHistogram(histClass, "centrFT0C_M11REFoverMp_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 1000000.0, VarManager::kM11REFoverMp);
        hm->AddHistogram(histClass, "centrFT0C_M1111REFoverMp_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 100000000000000.0, VarManager::kM1111REFoverMp);
        hm->AddHistogram(histClass, "centrFT0C_M11M1111REFoverMp_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 1000, 0.0, 10000000000000000.0, VarManager::kM11M1111REFoverMp);
        hm->AddHistogram(histClass, "centrFT0C_Corr2REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR2REFbydimuons, VarManager::kM11REFoverMp);
        hm->AddHistogram(histClass, "centrFT0C_Corr4REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR4REFbydimuons, VarManager::kM1111REFoverMp);
        hm->AddHistogram(histClass, "centrFT0C_Corr2Corr4REF_ev", "", true, 100, 0.0, 100.0, VarManager::kCentFT0C, 250, -1.0, 1.0, VarManager::kCORR2CORR4REF, VarManager::kM11M1111REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M11REFoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM11REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M1111REFoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM1111REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M01POIoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM01POIoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M0111POIoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM0111POIoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M11M1111REFoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM11M1111REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M01M0111overMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM01M0111overMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M11M0111overMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM11M0111overMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_M11M01REFoverMp", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kM11M01overMp);
        hm->AddHistogram(histClass, "Mass_Pt_Rapidity_CentFT0C", "", 4, var, bins, minBins, maxBins, 0, -1, kTRUE);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REF", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2REFbydimuons, VarManager::kM11REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4REF", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR4REFbydimuons, VarManager::kM1111REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2POI", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2POI, VarManager::kM01POIoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4POI", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR4POI, VarManager::kM0111POIoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFCorr4REF", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2CORR4REF, VarManager::kM11M1111REFoverMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2POICorr4POI", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2POICORR4POI, VarManager::kM01M0111overMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFCorr4POI", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2REFCORR4POI, VarManager::kM11M0111overMp);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFCorr2POI", "", true, 250, 0.0, 5.0, VarManager::kMass, 60, 0.0, 30.0, VarManager::kPt, 18, 0.0, 90.0, VarManager::kCentFT0C, "", "", "", VarManager::kCORR2REFCORR2POI, VarManager::kM11M01overMp);
      }
      if (subGroupStr.Contains("singlecumulant")) {
        int var[4] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C};
        int bins[4] = {250, 60, 6, 18};
        double minBins[4] = {0.0, 0.0, 2.5, 0.0};
        double maxBins[4] = {5.0, 30.0, 4.0, 90.0};
        hm->AddHistogram(histClass, "Mass_Pt_Rapidity_CentFT0C", "", 4, var, bins, minBins, maxBins, 0, -1, kTRUE);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFminus", "", true, 60, 0.0, 30.0, VarManager::kPt2, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR2REFbydimuons, "", "", "", VarManager::kNothing, VarManager::kM11REFoverMpminus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4REFminus", "", true, 60, 0.0, 30.0, VarManager::kPt2, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR4REFbydimuons, "", "", "", VarManager::kNothing, VarManager::kM1111REFoverMpminus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2REFplus", "", true, 60, 0.0, 30.0, VarManager::kPt1, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR2REFbydimuons, "", "", "", VarManager::kNothing, VarManager::kM11REFoverMpplus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4REFplus", "", true, 60, 0.0, 30.0, VarManager::kPt1, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR4REFbydimuons, "", "", "", VarManager::kNothing, VarManager::kM1111REFoverMpplus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2POIminus", "", true, 60, 0.0, 30.0, VarManager::kPt2, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR2POIminus, "", "", "", VarManager::kNothing, VarManager::kM01POIoverMpminus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4POIminus", "", true, 60, 0.0, 30.0, VarManager::kPt2, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR4POIminus, "", "", "", VarManager::kNothing, VarManager::kM0111POIoverMpminus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr2POIplus", "", true, 60, 0.0, 30.0, VarManager::kPt1, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR2POIplus, "", "", "", VarManager::kNothing, VarManager::kM01POIoverMpplus);
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_Corr4POIplus", "", true, 60, 0.0, 30.0, VarManager::kPt1, 18, 0.0, 90.0, VarManager::kCentFT0C, 0, 0.0, 1.0, VarManager::kCORR4POIplus, "", "", "", VarManager::kNothing, VarManager::kM0111POIoverMpplus);
      }
      if (subGroupStr.Contains("res-flow-dimuon")) {
        int varV2[6] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C, VarManager::kR2SP_AB, VarManager::kR2EP_AB};
        // int varV3[6] = {VarManager::kMass, VarManager::kPt, VarManager::kRap, VarManager::kCentFT0C, VarManager::kR3SP, VarManager::kR3EP}; // removed temporarily

        int bins[6] = {125, 60, 6, 18, 200, 40};
        double minBins[6] = {0.0, 0.0, 2.5, 0.0, -10.0, -2.0};
        double maxBins[6] = {5.0, 30.0, 4.0, 90.0, 10.0, 2.0};
        hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_R2", "", 6, varV2, bins, minBins, maxBins, 0, -1, kTRUE);
        // hm->AddHistogram(histClass, "Mass_Pt_centrFT0C_R3", "", 6, varV3, bins, minBins, maxBins, 0, -1, kTRUE); // removed temporarily
      }
      if (subGroupStr.Contains("z-boson")) {
        hm->AddHistogram(histClass, "MassZboson", "", false, 240, 20.0, 140.0, VarManager::kMass);
      }
      if (subGroupStr.Contains("dimuon-vtxncontrib")) {
        hm->AddHistogram(histClass, "MassMult", "", false, 750, 0.0, 15.0, VarManager::kMass, 301, -0.5, 300.5, VarManager::kVtxNcontrib);
        hm->AddHistogram(histClass, "MassVtxZMult", "", false, 300, 0.0, 6.0, VarManager::kMass, 60, -15.0, 15.0, VarManager::kVtxZ, 100, 0.0, 100.0, VarManager::kVtxNcontrib);
      }
    } else if (subGroupStr.Contains("electronmuon")) {
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kMass);
      hm->AddHistogram(histClass, "Pt", "", false, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Rapidity", "", false, 500, -1.0, 4.0, VarManager::kRap);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 750, 0.0, 30.0, VarManager::kMass, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass_Rapidity", "", false, 750, 0.0, 30.0, VarManager::kMass, 500, -1.0, 4.0, VarManager::kRap);
      hm->AddHistogram(histClass, "Mass_VtxZ", "", true, 30, -15.0, 15.0, VarManager::kVtxZ, 750, 0.0, 30.0, VarManager::kMass);
      hm->AddHistogram(histClass, "DeltaPhiPair", "", false, 130, -6.5, 6.5, VarManager::kDeltaPhiPair);
    }
    if (subGroupStr.Contains("correlation-emu")) {
      hm->AddHistogram(histClass, "DeltaPhiPair2_DeltaEtaPair2", "", false, 600, -o2::constants::math::PIHalf, 1.5 * o2::constants::math::PI, VarManager::kDeltaPhiPair2, 350, 1.5, 5.0, VarManager::kDeltaEtaPair2);
      hm->AddHistogram(histClass, "DeltaPhiPair2_Pt", "", false, 600, -o2::constants::math::PIHalf, 1.5 * o2::constants::math::PI, VarManager::kDeltaPhiPair2, 200, 0.0, 20.0, VarManager::kPt);
    }
    if (subGroupStr.Contains("dielectrons")) {
      if (subGroupStr.Contains("prefilter")) {
        hm->AddHistogram(histClass, "MassLow_OpeningAngle", "", false, 150, 0., 0.15, VarManager::kMass, 80, 0., 0.8, VarManager::kOpeningAngle);
      }
      if (subGroupStr.Contains("phiv")) {
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, 20, 0.0, 0.2, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 100, 0.0, o2::constants::math::PI, VarManager::kPairPhiv);
      }
      if (subGroupStr.Contains("double-phi-v")) {
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, 20, 0.0, 0.2, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 100, 0.0, o2::constants::math::PI, VarManager::kPairPhiv, "", "", "", -1, -1, true);
      }
      if (subGroupStr.Contains("largemass-phi-v")) {
        // binning for mee at large scales:
        // every 10 MeV from 0 to 0.2 GeV/c2
        // every 100 MeV from 0.2 to 1. GeV/c2
        // every 500 GeV from 1 to 5 GeV/c2
        double mee_bins[37];
        for (int i = 0; i <= 20; i++)
          mee_bins[i] = 0.01 * i;
        for (int i = 1; i <= 8; i++)
          mee_bins[20 + i] = 0.2 + 0.1 * i;
        for (int i = 1; i <= 8; i++)
          mee_bins[28 + i] = 1. + 0.5 * i;
        int nbins_mee = sizeof(mee_bins) / sizeof(*mee_bins) - 1;

        // binning for ptee at large scales:
        // every 0.2 GeV/c from 0 to 10 GeV/c
        double ptee_bins[51];
        for (int i = 0; i <= 50; i++)
          ptee_bins[i] = 0.2 * i;
        int nbins_ptee = sizeof(ptee_bins) / sizeof(*ptee_bins) - 1;

        // binning for phiv:
        // steps of size pi/100
        double phiv_bins[101];
        for (int i = 0; i <= 100; i++)
          phiv_bins[i] = o2::constants::math::PI / 100. * i;
        int nbins_phiv = sizeof(phiv_bins) / sizeof(*phiv_bins) - 1;

        // 3D histo
        hm->AddHistogram(histClass, "Mass_Pt_PhiV", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_ptee, ptee_bins, VarManager::kPt, nbins_phiv, phiv_bins, VarManager::kPairPhiv);
      }
      if (subGroupStr.Contains("meeptee")) {
        hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt);
      }
      if (subGroupStr.Contains("double-mee-ptee")) {
        hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 0, 0, 0, -1, "", "", "", -1, -1, true);
      }
      if (subGroupStr.Contains("lmee")) {
        hm->AddHistogram(histClass, "Mass_QuadDCAsigXY", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 20.0, VarManager::kQuadDCAsigXY);
        hm->AddHistogram(histClass, "Mass_QuadDCAsigZ", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 20.0, VarManager::kQuadDCAsigZ);
        hm->AddHistogram(histClass, "Mass_Pt_QuadDCAsigXYZ", "", false, 500, 0.0, 5.0, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPt, 50, 0.0, 20.0, VarManager::kQuadDCAsigXYZ);
        hm->AddHistogram(histClass, "Mass_Pt_SignQuadDCAsigXY", "", false, 50, 0.0, 5.0, VarManager::kMass, 50, 0.0, 10.0, VarManager::kPt, 100, -20.0, 20.0, VarManager::kSignQuadDCAsigXY);
      }
      if (subGroupStr.Contains("largescale")) {
        // binning for mee at large scales:
        // every 10 MeV from 0 to 1.1 GeV/c2
        // every 50 MeV from 1.1 to 2.7 GeV/c2
        // every 10 MeV from 2.7 to 3.2 GeV/c2
        // every 50 MeV from 3.2 to 12 GeV/c2
        double mee_bins[369];
        for (int i = 0; i <= 110; i++)
          mee_bins[i] = 0.01 * i;
        for (int i = 1; i <= 32; i++)
          mee_bins[110 + i] = 1.1 + 0.05 * i;
        for (int i = 1; i <= 50; i++)
          mee_bins[142 + i] = 2.7 + 0.01 * i;
        for (int i = 1; i <= 176; i++)
          mee_bins[192 + i] = 3.2 + 0.05 * i;
        int nbins_mee = sizeof(mee_bins) / sizeof(*mee_bins) - 1;

        // binning for ptee at large scales:
        // every 0.1 GeV/c from 0 to 10 GeV/c
        // every 0.5 GeV/c from 10 to 30 GeV/c
        double ptee_bins[201];
        for (int i = 0; i <= 100; i++)
          ptee_bins[i] = 0.1 * i;
        for (int i = 1; i <= 100; i++)
          ptee_bins[100 + i] = 10 + 0.2 * i;
        int nbins_ptee = sizeof(ptee_bins) / sizeof(*ptee_bins) - 1;

        // binning for dca at large scales:
        // every 0.1 sigma from 0 to 5 sigma
        // every 0.5 sigma from 5 to 10 sigma
        // every 1.0 sigma from 10 to 40 sigma
        double dca_bins[91];
        for (int i = 0; i <= 50; i++)
          dca_bins[i] = 0.1 * i;
        for (int i = 1; i <= 10; i++)
          dca_bins[50 + i] = 5 + 0.5 * i;
        for (int i = 1; i <= 30; i++)
          dca_bins[60 + i] = 10 + 1 * i;
        int nbins_dca = sizeof(dca_bins) / sizeof(*dca_bins) - 1;

        hm->AddHistogram(histClass, "Mass_QuadDCAsigXY", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_dca, dca_bins, VarManager::kQuadDCAsigXY);
        hm->AddHistogram(histClass, "Mass_QuadDCAsigZ", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_dca, dca_bins, VarManager::kQuadDCAsigZ);
        hm->AddHistogram(histClass, "Mass_Pt_QuadDCAsigXYZ", "", false, nbins_mee, mee_bins, VarManager::kMass, nbins_ptee, ptee_bins, VarManager::kPt, nbins_dca, dca_bins, VarManager::kQuadDCAsigXYZ);
      }
    }
    if (subGroupStr.Contains("opencharm")) {
      if (subGroupStr.Contains("dmeson")) {
        hm->AddHistogram(histClass, "MassD0region", "", false, 140, 1.5, 2.2, VarManager::kMass);
        hm->AddHistogram(histClass, "MassD0region_Pt", "", false, 70, 1.5, 2.2, VarManager::kMass, 160, 0., 20., VarManager::kPt);
        hm->AddHistogram(histClass, "MassD0region_Rapidity", "", false, 140, 1.5, 2.2, VarManager::kMass, 10, -0.8, 0.8, VarManager::kRap);
        hm->AddHistogram(histClass, "MassD0region_eta", "", false, 140, 1.5, 2.2, VarManager::kMass, 40, -2., 2., VarManager::kEta);
        hm->AddHistogram(histClass, "MassD0region_TauxyzProj", "", false, 140, 1.5, 2.2, VarManager::kMass, 200, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
        hm->AddHistogram(histClass, "MassD0region_TauxyProj", "", false, 140, 1.5, 2.2, VarManager::kMass, 200, -0.03, 0.03, VarManager::kVertexingTauxyProjected);
        hm->AddHistogram(histClass, "MassD0region_CosPointing", "", false, 140, 1.5, 2.2, VarManager::kMass, 200, -1.0, 1.0, VarManager::kCosPointingAngle);
        hm->AddHistogram(histClass, "MassD0region_VtxNContribReal", "", false, 140, 1.5, 2.2, VarManager::kMass, 50, 0, 50, VarManager::kVtxNcontribReal);
      }
      if (subGroupStr.Contains("3d-mass-histograms")) {
        hm->AddHistogram(histClass, "MassD0region_Pt_TauxyzProj", "", false, 140, 1.5, 2.2, VarManager::kMass, 80, 0., 20., VarManager::kPt, 300, 0., 0.03, VarManager::kVertexingTauxyzProjected);
        hm->AddHistogram(histClass, "MassD0region_Pt_CosPointing", "", false, 140, 1.5, 2.2, VarManager::kMass, 80, 0., 20., VarManager::kPt, 100, -1., 0., VarManager::kCosPointingAngle);
        hm->AddHistogram(histClass, "MassD0region_Rapidity_AveragePt", "", true, 140, 1.5, 2.2, VarManager::kMass, 10, -0.8, 0.8, VarManager::kRap, 150, 0.0, 30.0, VarManager::kPt);
        hm->AddHistogram(histClass, "MassD0region_Pt_ITStrackOccupancy", "Pair mass vs pair Pt vs event ITS occupancy", false, 70, 1.5, 2.2, VarManager::kMass, 160, 0., 20., VarManager::kPt, 200, 0., 20000., VarManager::kTrackOccupancyInTimeRange);
        hm->AddHistogram(histClass, "MassD0region_TPCnSigKa_pIN", "Pair mass vs kaon cand. pIN vs kaon cand. TPC n-#sigma(K)", false, 140, 1.5, 2.2, VarManager::kMass, 100, 0.0, 10.0, VarManager::kPin_leg1, 20, -5.0, 5.0, VarManager::kTPCnSigmaKa_leg1);
      }
      if (subGroupStr.Contains("lambdac")) {
        hm->AddHistogram(histClass, "MassLambdacRegion", "", false, 50, 2.15, 2.4, VarManager::kMass);
        hm->AddHistogram(histClass, "MassLambdacRegion_Pt", "", false, 50, 2.15, 2.4, VarManager::kMass, 40, 0.0, 20.0, VarManager::kPt);
        hm->AddHistogram(histClass, "MassLambdacRegion_eta", "", false, 50, 2.15, 2.4, VarManager::kMass, 40, -2., 2., VarManager::kEta);
        hm->AddHistogram(histClass, "MassLambdacRegion_TauxyzProj", "", false, 50, 2.15, 2.4, VarManager::kMass, 1000, -0.03, 0.03, VarManager::kVertexingTauxyzProjected);
      }
    }
  }

  if (!groupStr.CompareTo("dilepton-track")) {
    if (subGroupStr.Contains("mixedevent")) { // for mixed event
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 40, 0.0, 20.0, VarManager::kPairMass, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kPairMass);
      hm->AddHistogram(histClass, "Pt", "", false, 750, 0.0, 30.0, VarManager::kPairPt);
    }
    if (subGroupStr.Contains("invmass")) {
      hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 125, 0.0, 5.0, VarManager::kPairMassDau);
      hm->AddHistogram(histClass, "Mass_Hadron", "", false, 125, 0.0, 5.0, VarManager::kMassDau);
      hm->AddHistogram(histClass, "Delta_Mass", "", false, 125, 0.0, 5.0, VarManager::kDeltaMass);
      hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Hadron", "", false, 125, 0.0, 5.0, VarManager::kPairMassDau, 125, 0.0, 5.0, VarManager::kMassDau);
      hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 120, 0.0, 30.0, VarManager::kPairPtDau);
      hm->AddHistogram(histClass, "Pt_Track", "", false, 120, 0.0, 30.0, VarManager::kPt);
      hm->AddHistogram(histClass, "Mass", "", false, 750, 0.0, 30.0, VarManager::kPairMass);
      hm->AddHistogram(histClass, "Pt", "", false, 750, 0.0, 30.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Mass_Pt", "", false, 100, 0.0, 20.0, VarManager::kPairMass, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Pt_Dilepton__Pt", "", false, 40, 0.0, 20.0, VarManager::kPairPtDau, 40, 0.0, 20.0, VarManager::kPairPt);
      hm->AddHistogram(histClass, "Pt_Track__Pt", "", false, 40, 0.0, 20.0, VarManager::kPt, 40, 0.0, 20.0, VarManager::kPairPt);
    }
    if (subGroupStr.Contains("vertexing")) {
      hm->AddHistogram(histClass, "UsedKF", "", false, 2, -0.5, 1.5, VarManager::kUsedKF);
      hm->AddHistogram(histClass, "KFMass", "", false, 750, 0.0, 30.0, VarManager::kKFMass);
      hm->AddHistogram(histClass, "Lz", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLz);
      hm->AddHistogram(histClass, "Lxy", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxy);
      hm->AddHistogram(histClass, "Lxyz", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyz);
      hm->AddHistogram(histClass, "Tauz", "", false, 4000, -0.01, 0.01, VarManager::kVertexingTauz);
      hm->AddHistogram(histClass, "Tauxy", "", false, 4000, -0.01, 0.01, VarManager::kVertexingTauxy);
      hm->AddHistogram(histClass, "LxyzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLxyzErr);
      hm->AddHistogram(histClass, "LzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingLzErr);
      hm->AddHistogram(histClass, "TauzErr", "", false, 100, 0.0, 10.0, VarManager::kVertexingTauzErr);
      hm->AddHistogram(histClass, "VtxingProcCode", "", false, 10, 0.0, 10.0, VarManager::kVertexingProcCode);
      hm->AddHistogram(histClass, "VtxingChi2PCA", "", false, 100, 0.0, 10.0, VarManager::kVertexingChi2PCA);
      hm->AddHistogram(histClass, "LzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLzProjected);
      hm->AddHistogram(histClass, "LxyProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyProjected);
      hm->AddHistogram(histClass, "LxyzProj", "", false, 1000, -2.0, 2.0, VarManager::kVertexingLxyzProjected);
      hm->AddHistogram(histClass, "TauzProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauzProjected);
      hm->AddHistogram(histClass, "TauxyProj", "", false, 4000, -0.5, 0.5, VarManager::kVertexingTauxyProjected);
      hm->AddHistogram(histClass, "CosPointingAngle", "", false, 100, 0.0, 1.0, VarManager::kCosPointingAngle);
      hm->AddHistogram(histClass, "DCAxyzBetweenProngs", "", false, 100, 0.0, 1.0, VarManager::kKFDCAxyzBetweenProngs);
    }
    if (subGroupStr.Contains("multidimentional-vertexing-histograms")) {
      hm->AddHistogram(histClass, "Mass_Tauxy", "", false, 75, 4.0, 7.0, VarManager::kPairMass, 40, -0.0, 0.02, VarManager::kVertexingTauxy);
      hm->AddHistogram(histClass, "Mass_cosPointing", "", false, 75, 4.0, 7.0, VarManager::kPairMass, 40, 0.0, 1.0, VarManager::kCosPointingAngle);

      const int kNvarsTripletCuts = 4;
      const int kInvMassNbins = 100;
      double InvMassBinLims[kInvMassNbins + 1];
      for (int i = 0; i <= kInvMassNbins; ++i)
        InvMassBinLims[i] = 4.0 + 0.02 * i;

      const int kPtNbins = 6;
      double PtBinLims[kPtNbins + 1] = {0., 2., 4., 6., 8., 10., 20.};
      const int kCosPointingAngleNbins = 5;
      double CosPointingAngleBinLims[kCosPointingAngleNbins + 1] = {0., 0.86, 0.90, 0.94, 0.98, 1.0};

      const int kTauNBins = 6;
      double TauBinLims[kTauNBins + 1] = {0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.3};

      TArrayD nCutsBinLimits[kNvarsTripletCuts];
      nCutsBinLimits[0] = TArrayD(kInvMassNbins + 1, InvMassBinLims);
      nCutsBinLimits[1] = TArrayD(kPtNbins + 1, PtBinLims);
      nCutsBinLimits[2] = TArrayD(kCosPointingAngleNbins + 1, CosPointingAngleBinLims);
      nCutsBinLimits[3] = TArrayD(kTauNBins + 1, TauBinLims);

      int varsTripletCuts[kNvarsTripletCuts] = {VarManager::kPairMass, VarManager::kPairPt, VarManager::kCosPointingAngle, VarManager::kVertexingTauxyProjected};
      hm->AddHistogram(histClass, "multidimentional-vertexing", "Invariant mass vs. pT vs. cosine of pointing angle vs. tau", kNvarsTripletCuts, varsTripletCuts, nCutsBinLimits);
    }
    if (subGroupStr.Contains("correlation")) {
      hm->AddHistogram(histClass, "DeltaEta_DeltaPhi", "", false, 20, -2.0, 2.0, VarManager::kDeltaEta, 50, -8.0, 8.0, VarManager::kDeltaPhi);
      hm->AddHistogram(histClass, "DeltaEta_DeltaPhiSym", "", false, 20, -2.0, 2.0, VarManager::kDeltaEta, 50, -8.0, 8.0, VarManager::kDeltaPhiSym);
    }
    if (subGroupStr.Contains("dilepton-hadron-array-correlation")) {
      const int kInvMassBins = 500;
      double InvMassBinLims[kInvMassBins + 1];
      for (int i = 0; i <= kInvMassBins; i++)
        InvMassBinLims[i] = 0 + i * 0.01;

      const int kDelEtaBins = 20;
      double DelEtaBinLims[kDelEtaBins + 1];
      for (int i = 0; i <= kDelEtaBins; i++)
        DelEtaBinLims[i] = -2 + i * 0.2;

      const int kDelPhiBins = 52;
      double DelPhiBinLims[] = {-1.69647, -1.57080, -1.44513, -1.31947, -1.19381, -1.06814, -0.94248, -0.81681, -0.69115, -0.56549, -0.43982, -0.31416, -0.18850, -0.06283, 0.06283, 0.18850, 0.31416, 0.43982, 0.56549, 0.69115, 0.81681, 0.94248, 1.06814, 1.19381, 1.31947, 1.44513, 1.57080, 1.69646, 1.82212, 1.94779, 2.07345, 2.19911, 2.32478, 2.45044, 2.57611, 2.70177, 2.82743, 2.95310, 3.07876, 3.20442, 3.33009, 3.45575, 3.58142, 3.70708, 3.83274, 3.95841, 4.08407, 4.20973, 4.33540, 4.46106, 4.58673, 4.71239, 4.8380600};

      const int kPtBins = 45;
      double PtBinLims[kPtBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.5, 5, 7.5, 10, 20};

      TArrayD nJPsiHadCorr[4];
      nJPsiHadCorr[0] = TArrayD(kInvMassBins + 1, InvMassBinLims);
      nJPsiHadCorr[1] = TArrayD(kDelEtaBins + 1, DelEtaBinLims);
      nJPsiHadCorr[2] = TArrayD(kDelPhiBins + 1, DelPhiBinLims);
      nJPsiHadCorr[3] = TArrayD(kPtBins + 1, PtBinLims);

      int varsJPsiHadCorr[4] = {VarManager::kPairMassDau, VarManager::kDeltaEta, VarManager::kDeltaPhi, VarManager::kPairPtDau};
      hm->AddHistogram(histClass, "InvMass_DelEta_DelPhi", "", 4, varsJPsiHadCorr, nJPsiHadCorr); // Without efficiency
      // hm->AddHistogram(histClass, "InvMass_DelEta_DelPhi", "", 4, varsJPsiHadCorr, nJPsiHadCorr, nullptr, VarManager::kJpsiHadronEff);
    }
    if (subGroupStr.Contains("opencharm")) {
      hm->AddHistogram(histClass, "Delta_Mass_DstarD0region", "", false, 50, 0.14, 0.16, VarManager::kDeltaMass);
    }
  }

  if (!groupStr.CompareTo("dilepton-charmhadron")) {
    if (subGroupStr.EqualTo("jpsitomumu")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 50, -4.5f, -2.0f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
    } else if (subGroupStr.EqualTo("jpsitoee")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 60, -1.5f, 1.5f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsi", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
    } else if (subGroupStr.EqualTo("dmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmeson", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    } else if (subGroupStr.EqualTo("jpsitomumudmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 50, -4.5f, -2.0f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    } else if (subGroupStr.EqualTo("jpsitoeedmeson")) {
      hm->AddHistogram(histClass, "hMassVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 300, 2.f, 5.f, VarManager::kMass);
      hm->AddHistogram(histClass, "hRapVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 60, -1.5f, 1.5f, VarManager::kRap);
      hm->AddHistogram(histClass, "hPhiVsPtJPsiWithDmeson", "", false, 100, 0.f, 50.f, VarManager::kPt, 180, -constants::math::PI, constants::math::PI, VarManager::kPhi);
      hm->AddHistogram(histClass, "hMassVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 200, 1.7f, 2.1f, VarManager::kMassCharmHadron);
      hm->AddHistogram(histClass, "hRapVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 60, -1.5f, 1.5f, VarManager::kRapCharmHadron);
      hm->AddHistogram(histClass, "hPhiVsPtVsBdtDmesonWithJPsi", "", false, 100, 0.f, 1.f, VarManager::kBdtCharmHadron, 100, 0.f, 50.f, VarManager::kPtCharmHadron, 180, 0., 2 * constants::math::PI, VarManager::kPhiCharmHadron);
    }
  }
  if (!groupStr.CompareTo("dilepton-dihadron")) {
    if (subGroupStr.EqualTo("xtojpsipipi")) {
      hm->AddHistogram(histClass, "hMass_X3872", "", false, 1000, 3.0, 5.0, VarManager::kQuadMass);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_X3872", "", false, 1000, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass);
      hm->AddHistogram(histClass, "hPt_X3872", "", false, 150, 0.0, 15.0, VarManager::kQuadPt);
      hm->AddHistogram(histClass, "hMass_Pt_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 150, 0.0, 15.0, VarManager::kQuadPt);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_Pt_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 150, 0.0, 15.0, VarManager::kQuadPt);
      hm->AddHistogram(histClass, "hCostheta_Jpsi_Dihadron", "", false, 100, -1.0, 1.0, VarManager::kCosthetaDileptonDitrack);
      hm->AddHistogram(histClass, "hPtDilepton_PtDihadron", "", false, 150, 0, 15.0, VarManager::kPairPt, 100, 0, 10, VarManager::kDitrackPt);
      hm->AddHistogram(histClass, "hPtDilepton_MassDihadron", "", false, 150, 0, 15.0, VarManager::kPairPt, 150, 0.0, 3.0, VarManager::kDitrackMass);
      hm->AddHistogram(histClass, "hQ_X3872", "", false, 150, 0.0, 3.0, VarManager::kQ);
      hm->AddHistogram(histClass, "hDeltaR1_X3872", "", false, 100, 0.0, 10.0, VarManager::kDeltaR1);
      hm->AddHistogram(histClass, "hDeltaR2_X3872", "", false, 100, 0.0, 10.0, VarManager::kDeltaR2);
      hm->AddHistogram(histClass, "hDeltaR_X3872", "", false, 100, 0.0, 10.0, VarManager::kDeltaR);
      hm->AddHistogram(histClass, "hMass_Q_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 150, 0.0, 3.0, VarManager::kQ);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_Q_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 150, 0.0, 3.0, VarManager::kQ);
      hm->AddHistogram(histClass, "hMass_DeltaR1_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 100, 0.0, 10.0, VarManager::kDeltaR1);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_DeltaR1_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 100, 0.0, 10.0, VarManager::kDeltaR1);
      hm->AddHistogram(histClass, "hMass_DeltaR2_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 100, 0.0, 10.0, VarManager::kDeltaR2);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_DeltaR2_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 100, 0.0, 10.0, VarManager::kDeltaR2);
      hm->AddHistogram(histClass, "hMass_X3872_MassDihadron", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 150, 0.0, 3.0, VarManager::kDitrackMass);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_X3872_MassDihadron", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 150, 0.0, 3.0, VarManager::kDitrackMass);
      hm->AddHistogram(histClass, "hRap_X3872", "", false, 1000, 0.0, 5.0, VarManager::kRap);
      hm->AddHistogram(histClass, "hMass_Rap_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 1000, 0.0, 5.0, VarManager::kRap);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_Rap_X3872", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 1000, 0.0, 5.0, VarManager::kRap);
      hm->AddHistogram(histClass, "hDCAxyTrack1", "", false, 100, -0.1, 0.1, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "hDCAzTrack1", "", false, 100, -0.1, 0.1, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_DCAxyTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 100, -0.1, 0.1, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_DCAzTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 100, -0.1, 0.1, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "hMass_DCAxyTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 100, -0.1, 0.1, VarManager::kTrackDCAxy);
      hm->AddHistogram(histClass, "hMass_DCAzTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 100, -0.1, 0.1, VarManager::kTrackDCAz);
      hm->AddHistogram(histClass, "hPtTrack1", "", false, 100, 0.0, 10.0, VarManager::kPt);
      hm->AddHistogram(histClass, "hMass_defaultDileptonMass_PtTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadDefaultDileptonMass, 100, 0.0, 10.0, VarManager::kPt);
      hm->AddHistogram(histClass, "hMass_PtTrack1", "", false, 100, 3.0, 5.0, VarManager::kQuadMass, 100, 0.0, 10.0, VarManager::kPt);
    }
  }
  if (!groupStr.CompareTo("dilepton-photon-mass")) {
    hm->AddHistogram(histClass, "Mass_Dilepton", "", false, 500, 0.0, 5.0, VarManager::kPairMassDau);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Dilepton_Mass_Photon", "", false, 250, 0.0, 5.0, VarManager::kPairMassDau, 250, 0.0, 5.0, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Pt_Dilepton", "", false, 2000, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Pt_Photon", "", false, 500, 0.0, 5.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Mass_DileptonPhoton", "", false, 4500, 0.0, 4.5, VarManager::kPairMass);
    hm->AddHistogram(histClass, "Pt_DileptonPhoton", "", false, 2000, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kPairMass, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass", "", false, 1500, 0.0, 1.5, VarManager::kDeltaMass);
    hm->AddHistogram(histClass, "DeltaMass_pt", "", false, 1000, 0.0, 1.0, VarManager::kDeltaMass, 3000, 0.0, 30.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kPairEta, 200, 0.0, 20.0, VarManager::kPairPt);
    hm->AddHistogram(histClass, "DeltaMass_Jpsi", "", false, 1500, 3, 4.5, VarManager::kDeltaMass_jpsi);
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
    hm->AddHistogram(histClass, "Eta_Pt_Dilepton", "", false, 100, -2.0, 2.0, VarManager::kDeltaEta, 200, 0.0, 20.0, VarManager::kPairPtDau);
    hm->AddHistogram(histClass, "Eta_Pt_Photon", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton1", "", false, 100, -2.0, 2.0, VarManager::kEta1, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta_Pt_lepton2", "", false, 100, -2.0, 2.0, VarManager::kEta2, 200, 0.0, 20.0, VarManager::kPt);
  }

  if (!groupStr.CompareTo("photon")) {
    hm->AddHistogram(histClass, "Pt_Photon", "p_{T} distribution", false, 4500, 0.0, 4.5, VarManager::kPt);
    hm->AddHistogram(histClass, "Eta", "#eta distribution", false, 500, -5.0, 5.0, VarManager::kEta);
    hm->AddHistogram(histClass, "Eta_Pt", "", false, 100, -2.0, 2.0, VarManager::kEta, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Phi", "#varphi distribution", false, 500, -2. * o2::constants::math::PI, 2. * o2::constants::math::PI, VarManager::kPhi);
    hm->AddHistogram(histClass, "Mass_Photon", "", false, 500, 0.0, 0.1, VarManager::kMassDau);
    hm->AddHistogram(histClass, "Mass_Pt", "", false, 500, 0.0, 5.0, VarManager::kMassDau, 200, 0.0, 20.0, VarManager::kPt);
    hm->AddHistogram(histClass, "Rapidity", "", false, 400, -4.0, 4.0, VarManager::kRap);
  }
}

//__________________________________________________________________
template <typename T>
bool o2::aod::dqhistograms::ValidateJSONHistogram(T hist)
{
  //
  // Validate JSON entry for this histogram
  //

  // The fields histClass, title and type are compulsory
  if (!hist->HasMember("histClass") || !hist->HasMember("title") || !hist->HasMember("type")) {
    LOG(fatal) << "Missing histClass, title or type fields";
    return false;
  }

  TString histTypeStr = hist->FindMember("type")->value.GetString();
  bool isTH1 = (histTypeStr.CompareTo("TH1") == 0);
  bool isTH2 = (histTypeStr.CompareTo("TH2") == 0);
  bool isTH3 = (histTypeStr.CompareTo("TH3") == 0);
  bool isTHn = (histTypeStr.CompareTo("THn") == 0);
  if (!(isTH1 || isTH2 || isTH3 || isTHn)) {
    LOG(fatal) << "The type field must be one of the TH1, TH2, TH3 or THn";
    return false;
  }
  // Check if the histogram uses constant binning
  bool isConstantBinning = true;
  if (!(hist->HasMember("xmin") && hist->HasMember("xmax"))) {
    isConstantBinning = false;
  }

  if (!isTHn && (!hist->HasMember("isProfile") || !hist->HasMember("nXbins") || !hist->HasMember("varX"))) {
    LOG(fatal) << "Missing isProfile, nXbins or varX information for histogram";
    return false;
  }
  bool isProfile = (hist->HasMember("isProfile") ? hist->FindMember("isProfile")->value.GetBool() : false);

  if (isConstantBinning) {
    if (!hist->HasMember("xmin") || !hist->HasMember("xmax")) {
      LOG(fatal) << "Missing xmin or xmax information for histogram";
      return false;
    }
    if (isTHn) {
      if (!hist->FindMember("xmin")->value.IsArray()) {
        LOG(fatal) << "xmin field should be an array of arrays";
        return false;
      }
      if (!hist->FindMember("xmax")->value.IsArray()) {
        LOG(fatal) << "xmax field should be an array of arrays";
        return false;
      }
    }
  } else {
    if (isTHn && !hist->HasMember("binLimits")) {
      LOG(fatal) << "Missing binLimits information for histogram";
      return false;
    }
    if (!isTHn && !hist->HasMember("xbins")) {
      LOG(fatal) << "Missing xbins information for histogram";
      return false;
    }
    if (isTHn && !hist->FindMember("binLimits")->value.IsArray()) {
      LOG(fatal) << "binLimits field should be an array of arrays";
      return false;
    }
    if (!isTHn && !hist->FindMember("xbins")->value.IsArray()) {
      LOG(fatal) << "xbins field should be an array";
      return false;
    }
  }
  if (isProfile && !hist->HasMember("varY")) {
    LOG(fatal) << "Missing varY information for histogram";
    return false;
  }

  if (isTHn) {
    if (!hist->HasMember("nDimensions") || !hist->HasMember("vars")) {
      LOG(fatal) << "Missing nDimensions or vars fields for histogram";
      return false;
    }
    if (isConstantBinning) {
      if (!hist->HasMember("nBins")) {
        LOG(fatal) << "Missing nBins field for histogram";
        return false;
      } else {
        if (!hist->FindMember("nBins")->value.IsArray()) {
          LOG(fatal) << "nBins field should be an array";
          return false;
        }
      }
    }
    if (hist->HasMember("axLabels") && !hist->FindMember("axLabels")->value.IsArray()) {
      LOG(fatal) << "axLabels field should be an array of strings";
      return false;
    }
  }

  if (isTH2 || isTH3) {
    if (!hist->HasMember("nYbins") || !hist->HasMember("varY")) {
      LOG(fatal) << "Missing nYbins or varY information for histogram";
      return false;
    }
    if (isConstantBinning && (!hist->HasMember("ymin") || !hist->HasMember("ymax"))) {
      LOG(fatal) << "Missing ymin or ymax information for histogram";
      return false;
    }
    if (!isConstantBinning && !hist->HasMember("ybins")) {
      LOG(fatal) << "Missing ybins information for histogram";
      return false;
    }
    if (!isConstantBinning && !hist->FindMember("xbins")->value.IsArray()) {
      LOG(fatal) << "ybins field should be an array";
    }

    if (isTH3) {
      if (!hist->HasMember("nZbins") || !hist->HasMember("varZ")) {
        LOG(fatal) << "Missing nZbins or varZ information for histogram";
        return false;
      }
      if (isConstantBinning && (!hist->HasMember("zmin") || !hist->HasMember("zmax"))) {
        LOG(fatal) << "Missing zmin or zmax information for histogram";
        return false;
      }
      if (!isConstantBinning && !hist->HasMember("zbins")) {
        LOG(fatal) << "Missing zbins information for histogram";
        return false;
      }
      if (!isConstantBinning && !hist->FindMember("zbins")->value.IsArray()) {
        LOG(fatal) << "zbins field should be an array";
      }
    }
  }
  if (isTH2 && isProfile && !hist->HasMember("varZ")) {
    LOG(fatal) << "Missing varZ information for histogram";
    return false;
  }
  if (isTH3 && isProfile && !hist->HasMember("varT")) {
    LOG(fatal) << "Missing varT information for histogram";
    return false;
  }

  if (!isTHn) {
    TString varX = hist->FindMember("varX")->value.GetString();
    if (VarManager::fgVarNamesMap.find(varX) == VarManager::fgVarNamesMap.end()) {
      LOG(fatal) << "Bad varX variable (" << hist->FindMember("varX")->value.GetString() << ") specified for histogram";
      return false;
    }
    if (hist->HasMember("varY") && (VarManager::fgVarNamesMap.find(hist->FindMember("varY")->value.GetString()) == VarManager::fgVarNamesMap.end())) {
      LOG(fatal) << "Bad varY variable (" << hist->FindMember("varY")->value.GetString() << ") specified for histogram";
      return false;
    }
    if (hist->HasMember("varZ") && (VarManager::fgVarNamesMap.find(hist->FindMember("varZ")->value.GetString()) == VarManager::fgVarNamesMap.end())) {
      LOG(fatal) << "Bad varZ variable (" << hist->FindMember("varZ")->value.GetString() << ") specified for histogram";
      return false;
    }
    if (hist->HasMember("varT") && (VarManager::fgVarNamesMap.find(hist->FindMember("varT")->value.GetString()) == VarManager::fgVarNamesMap.end())) {
      LOG(fatal) << "Bad varT variable (" << hist->FindMember("varT")->value.GetString() << ") specified for histogram";
      return false;
    }
    if (hist->HasMember("varW") && (VarManager::fgVarNamesMap.find(hist->FindMember("varW")->value.GetString()) == VarManager::fgVarNamesMap.end())) {
      LOG(fatal) << "Bad varW variable (" << hist->FindMember("varW")->value.GetString() << ") specified for histogram";
      return false;
    }
  }
  if (isTHn) {
    for (auto& v : hist->FindMember("vars")->value.GetArray()) {
      if (VarManager::fgVarNamesMap.find(v.GetString()) == VarManager::fgVarNamesMap.end()) {
        LOG(fatal) << "Bad variable in vars (" << v.GetString() << ") specified for histogram";
        return false;
      }
    }
  }

  return true;
}

//__________________________________________________________________
void o2::aod::dqhistograms::AddHistogramsFromJSON(HistogramManager* hm, const char* json)
{
  //
  // Add histograms to already existing histogram classes from a JSON formatted string
  //   The JSON is expected to contain a list of objects, with each object containing the fields needed
  //    to define a histogram via the HistogramManager::AddHistogram() functions

  LOG(info) << "========================================== interpreting JSON for adding histograms";
  LOG(info) << "      json string is: " << json;

  TString jsonStr = json;
  if (jsonStr == "") {
    // No histograms to add
    return;
  }

  rapidjson::Document document;
  rapidjson::ParseResult ok = document.Parse(json);
  if (!ok) {
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")";
    TString str = "";
    for (int i = ok.Offset() - 30; i < static_cast<int>(ok.Offset()) + 50; i++) {
      if ((i >= 0) && (i < static_cast<int>(strlen(json)))) {
        str += json[i];
      }
    }
    LOG(fatal) << "**** Parsing error is somewhere here: " << str.Data() << endl;
    return;
  }

  for (rapidjson::Value::ConstMemberIterator it = document.MemberBegin(); it != document.MemberEnd(); it++) {

    const char* histName = it->name.GetString();
    LOG(info) << "Configuring histogram " << histName;
    const auto& hist = it->value;
    if (!ValidateJSONHistogram(&hist)) {
      LOG(fatal) << "Histogram not properly defined in the JSON file. Skipping it";
      continue;
    }

    TString histTypeStr = hist.FindMember("type")->value.GetString();
    bool isTH2 = (histTypeStr.CompareTo("TH2") == 0);
    bool isTH3 = (histTypeStr.CompareTo("TH3") == 0);
    bool isTHn = (histTypeStr.CompareTo("THn") == 0);
    bool isConstantBinning = true;
    if (!(hist.HasMember("xmin") && hist.HasMember("xmax"))) {
      isConstantBinning = false;
    }

    const char* histClass = hist.FindMember("histClass")->value.GetString();
    const char* title = hist.FindMember("title")->value.GetString();

    if (isTHn) {
      int nDimensions = hist.FindMember("nDimensions")->value.GetInt();
      LOG(debug) << "nDimensions: " << nDimensions;

      int* vars = new int[nDimensions];
      int iDim = 0;
      for (auto& v : hist.FindMember("vars")->value.GetArray()) {
        LOG(debug) << "iDim " << iDim << ": " << v.GetString();
        vars[iDim++] = VarManager::fgVarNamesMap[v.GetString()];
      }

      int* nBins = nullptr;
      double* xmin = nullptr;
      double* xmax = nullptr;
      TArrayD* binLimits = nullptr;
      if (isConstantBinning) {
        nBins = new int[nDimensions];
        xmin = new double[nDimensions];
        xmax = new double[nDimensions];
        int iDim = 0;
        for (auto& v : hist.FindMember("nBins")->value.GetArray()) {
          nBins[iDim++] = v.GetInt();
          LOG(debug) << "nBins " << iDim << ": " << nBins[iDim - 1];
        }
        iDim = 0;
        for (auto& v : hist.FindMember("xmin")->value.GetArray()) {
          xmin[iDim++] = v.GetDouble();
          LOG(debug) << "xmin " << iDim << ": " << xmin[iDim - 1];
        }
        iDim = 0;
        for (auto& v : hist.FindMember("xmax")->value.GetArray()) {
          xmax[iDim++] = v.GetDouble();
          LOG(debug) << "xmax " << iDim << ": " << xmax[iDim - 1];
        }
      } else {
        int iDim = 0;
        binLimits = new TArrayD[nDimensions];
        for (auto& v : hist.FindMember("binLimits")->value.GetArray()) {
          double* lims = new double[v.GetArray().Size()];
          int iElem = 0;
          for (auto& lim : v.GetArray()) {
            lims[iElem++] = lim.GetDouble();
          }
          binLimits[iDim++] = TArrayD(v.GetArray().Size(), lims);
        }
      }

      TString* axLabels = nullptr;
      if (hist.HasMember("axLabels")) {
        axLabels = new TString[hist.FindMember("axLabels")->value.GetArray().Size()];
        int iDim = 0;
        for (auto& v : hist.FindMember("axLabels")->value.GetArray()) {
          axLabels[iDim++] = v.GetString();
        }
      }

      int varW = (hist.HasMember("varW") ? VarManager::fgVarNamesMap[hist.FindMember("varW")->value.GetString()] : -1);
      bool useSparse = (hist.HasMember("useSparse") ? hist.FindMember("useSparse")->value.GetBool() : false);
      bool isDouble = (hist.HasMember("isDouble") ? hist.FindMember("isDouble")->value.GetBool() : false);

      if (isConstantBinning) {
        hm->AddHistogram(histClass, histName, title, nDimensions, vars, nBins, xmin, xmax, axLabels, varW, useSparse, isDouble);
      } else {
        hm->AddHistogram(histClass, histName, title, nDimensions, vars, binLimits, axLabels, varW, useSparse, isDouble);
      }

    } else { // TH1, TH2 or TH3

      LOG(debug) << "is TH1, TH2 or TH3 ";

      bool isProfile = hist.FindMember("isProfile")->value.GetBool();
      LOG(debug) << "isProfile: " << isProfile;

      int nXbins = hist.FindMember("nXbins")->value.GetInt();
      LOG(debug) << "nXbins: " << nXbins;

      const char* varX = hist.FindMember("varX")->value.GetString();
      LOG(debug) << "varX: " << varX;

      double xmin = (hist.HasMember("xmin") ? hist.FindMember("xmin")->value.GetDouble() : 0.0);
      LOG(debug) << "xmin: " << xmin;

      double xmax = (hist.HasMember("xmax") ? hist.FindMember("xmax")->value.GetDouble() : 0.0);
      LOG(debug) << "xmax: " << xmax;

      std::vector<double> xbinsVec;
      if (hist.HasMember("xbins")) {
        LOG(debug) << "xbins: ";
        for (auto& v : hist.FindMember("xbins")->value.GetArray()) {
          xbinsVec.push_back(v.GetDouble());
          LOG(debug) << v.GetDouble();
        }
      }

      const char* varY = (hist.HasMember("varY") ? hist.FindMember("varY")->value.GetString() : "kNothing");
      LOG(debug) << "varY: " << varY;

      int nYbins = (hist.HasMember("nYbins") ? hist.FindMember("nYbins")->value.GetInt() : 0);
      LOG(debug) << "nYbins: " << nYbins;

      double ymin = (hist.HasMember("ymin") ? hist.FindMember("ymin")->value.GetDouble() : 0.0);
      LOG(debug) << "ymin: " << ymin;

      double ymax = (hist.HasMember("ymax") ? hist.FindMember("ymax")->value.GetDouble() : 0.0);
      LOG(debug) << "ymax: " << ymax;

      std::vector<double> ybinsVec;
      if (hist.HasMember("ybins")) {
        LOG(debug) << "ybins: ";
        for (auto& v : hist.FindMember("ybins")->value.GetArray()) {
          ybinsVec.push_back(v.GetDouble());
          LOG(debug) << v.GetDouble();
        }
      }

      const char* varZ = (hist.HasMember("varZ") ? hist.FindMember("varZ")->value.GetString() : "kNothing");
      LOG(debug) << "varZ: " << varZ;

      int nZbins = (hist.HasMember("nZbins") ? hist.FindMember("nZbins")->value.GetInt() : 0);
      LOG(debug) << "nZbins: " << nZbins;

      double zmin = (hist.HasMember("zmin") ? hist.FindMember("zmin")->value.GetDouble() : 0.0);
      LOG(debug) << "zmin: " << zmin;

      double zmax = (hist.HasMember("zmax") ? hist.FindMember("zmax")->value.GetDouble() : 0.0);
      LOG(debug) << "zmax: " << zmax;

      std::vector<double> zbinsVec;
      if (hist.HasMember("zbins")) {
        LOG(debug) << "zbins: ";
        for (auto& v : hist.FindMember("zbins")->value.GetArray()) {
          zbinsVec.push_back(v.GetDouble());
          LOG(debug) << v.GetDouble();
        }
      }

      const char* xLabels = (hist.HasMember("xLabels") ? hist.FindMember("xLabels")->value.GetString() : "");
      LOG(debug) << "xLabels: " << xLabels;

      const char* yLabels = (hist.HasMember("yLabels") ? hist.FindMember("yLabels")->value.GetString() : "");
      LOG(debug) << "yLabels: " << yLabels;

      const char* zLabels = (hist.HasMember("zLabels") ? hist.FindMember("zLabels")->value.GetString() : "");
      LOG(debug) << "zLabels: " << zLabels;

      const char* varT = (hist.HasMember("varT") ? hist.FindMember("varT")->value.GetString() : "kNothing");
      LOG(debug) << "varT: " << varT;

      const char* varW = (hist.HasMember("varW") ? hist.FindMember("varW")->value.GetString() : "kNothing");
      LOG(debug) << "varW: " << varW;

      bool isdouble = (hist.HasMember("isdouble") ? hist.FindMember("isdouble")->value.GetBool() : false);
      LOG(debug) << "isdouble: " << isdouble;

      bool isFillLabelx = (hist.HasMember("isFillLabelx") ? hist.FindMember("isFillLabelx")->value.GetBool() : false);
      LOG(debug) << "isFillLabelx: " << isFillLabelx;

      if (isConstantBinning) {
        hm->AddHistogram(histClass, histName, title, isProfile,
                         nXbins, xmin, xmax, VarManager::fgVarNamesMap[varX],
                         nYbins, ymin, ymax, VarManager::fgVarNamesMap[varY],
                         nZbins, zmin, zmax, VarManager::fgVarNamesMap[varZ],
                         xLabels, yLabels, zLabels,
                         VarManager::fgVarNamesMap[varT], VarManager::fgVarNamesMap[varW], isdouble, isFillLabelx);
      } else {
        int xBinsSize = xbinsVec.size();
        if (xBinsSize != (nXbins + 1)) {
          LOG(fatal) << "Histogram not properly defined in the JSON file. Wrong x binning for histogram";
          continue;
        }
        double* xbins = new double[xbinsVec.size()];
        std::copy(xbinsVec.begin(), xbinsVec.end(), xbins);

        double* ybins = nullptr;
        if (isTH2 || isTH3) {
          if (static_cast<int>(ybinsVec.size()) != (nYbins + 1)) {
            LOG(fatal) << "Histogram not properly defined in the JSON file. Wrong y binning for histogram";
            continue;
          }
          ybins = new double[ybinsVec.size()];
          std::copy(ybinsVec.begin(), ybinsVec.end(), ybins);
        }

        double* zbins = nullptr;
        if (isTH3) {
          if (static_cast<float>(zbinsVec.size()) != (nZbins + 1)) {
            LOG(fatal) << "Histogram not properly defined in the JSON file. Wrong z binning for histogram";
            continue;
          }
          zbins = new double[zbinsVec.size()];
          std::copy(zbinsVec.begin(), zbinsVec.end(), zbins);
        }
        hm->AddHistogram(histClass, histName, title, isProfile,
                         nXbins, xbins, VarManager::fgVarNamesMap[varX],
                         nYbins, ybins, VarManager::fgVarNamesMap[varY],
                         nZbins, zbins, VarManager::fgVarNamesMap[varZ],
                         xLabels, yLabels, zLabels,
                         VarManager::fgVarNamesMap[varT], VarManager::fgVarNamesMap[varW], isdouble, isFillLabelx);
      } // end if (!isTHn)
    }
  }
}
