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

/// \header file for histograms
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EVENTHISTOGRAMS_H_
#define PWGEM_DILEPTON_UTILS_EVENTHISTOGRAMS_H_

#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::aod::pwgem::dilepton::utils::eventhistogram
{
const int nbin_ev = 21;
template <const int nmod = -1>
void addEventHistograms(HistogramRegistry* fRegistry)
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1D, {{nbin_ev, 0.5, nbin_ev + 0.5}}, false);
  hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
  hCollisionCounter->GetXaxis()->SetBinLabel(2, "FT0AND");
  hCollisionCounter->GetXaxis()->SetBinLabel(3, "No TF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(4, "No ITS ROF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(5, "No Same Bunch Pileup");
  hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
  hCollisionCounter->GetXaxis()->SetBinLabel(7, "Is Vertex ITS-TPC");
  hCollisionCounter->GetXaxis()->SetBinLabel(8, "Is Vertex ITS-TPC-TRD");
  hCollisionCounter->GetXaxis()->SetBinLabel(9, "Is Vertex ITS-TPC-TOF");
  hCollisionCounter->GetXaxis()->SetBinLabel(10, "sel8");
  hCollisionCounter->GetXaxis()->SetBinLabel(11, "|Z_{vtx}| < 10 cm");
  hCollisionCounter->GetXaxis()->SetBinLabel(12, "NoCollInTimeRangeStandard");
  hCollisionCounter->GetXaxis()->SetBinLabel(13, "NoCollInTimeRangeStrict");
  hCollisionCounter->GetXaxis()->SetBinLabel(14, "NoCollInRofStandard");
  hCollisionCounter->GetXaxis()->SetBinLabel(15, "NoCollInRofStrict");
  hCollisionCounter->GetXaxis()->SetBinLabel(16, "NoHighMultCollInPrevRof");
  hCollisionCounter->GetXaxis()->SetBinLabel(17, "IsGoodITSLayer3");
  hCollisionCounter->GetXaxis()->SetBinLabel(18, "IsGoodITSLayer0123");
  hCollisionCounter->GetXaxis()->SetBinLabel(19, "IsGoodITSLayersAll");
  hCollisionCounter->GetXaxis()->SetBinLabel(20, "Calibrated Q vector");
  hCollisionCounter->GetXaxis()->SetBinLabel(21, "accepted");

  const AxisSpec axis_cent_ft0m{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                "centrality FT0M (%)"};

  const AxisSpec axis_cent_ft0a{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                "centrality FT0A (%)"};

  const AxisSpec axis_cent_ft0c{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                "centrality FT0C (%)"};

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1D, {{100, -50, +50}}, false);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1D, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1D, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2D, {{200, 0, 200000}, {60, 0, 60000}}, false);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1D, {{axis_cent_ft0a}}, false);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1D, {{axis_cent_ft0c}}, false);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1D, {{axis_cent_ft0m}}, false);
  fRegistry->add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2D, {{110, 0, 110}, {600, 0, 6000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2D, {{60, 0, 60000}, {600, 0, 6000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", kTH2D, {{60, 0, 60000}, {200, 0, 20000}}, false);
  fRegistry->add("Event/before/hNTracksPVvsOccupancy", "hNTracksPVvsOccupancy;N_{track} to PV;N_{track} in time range", kTH2D, {{600, 0, 6000}, {200, 0, 20000}}, false);
  fRegistry->add("Event/before/hCorrOccupancy", "occupancy correlation;FT0C occupancy;track occupancy", kTH2D, {{200, 0, 200000}, {200, 0, 20000}}, false);

  if constexpr (nmod == 2) {
    fRegistry->add("Event/before/hQ2xFT0M_CentFT0C", "hQ2xFT0M_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0M}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0M_CentFT0C", "hQ2yFT0M_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0M}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0A_CentFT0C", "hQ2xFT0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0A_CentFT0C", "hQ2yFT0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0C_CentFT0C", "hQ2xFT0C_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0C}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0C_CentFT0C", "hQ2yFT0C_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0C}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFV0A_CentFT0C", "hQ2xFV0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FV0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFV0A_CentFT0C", "hQ2yFV0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FV0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBPos_CentFT0C", "hQ2xBPos_CentFT0C;centrality FT0C (%);Q_{2,x}^{BPos}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBPos_CentFT0C", "hQ2yBPos_CentFT0C;centrality FT0C (%);Q_{2,y}^{BPos}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBNeg_CentFT0C", "hQ2xBNeg_CentFT0C;centrality FT0C (%);Q_{2,x}^{BNeg}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBNeg_CentFT0C", "hQ2yBNeg_CentFT0C;centrality FT0C (%);Q_{2,y}^{BNeg}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBTot_CentFT0C", "hQ2xBTot_CentFT0C;centrality FT0C (%);Q_{2,x}^{BTot}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBTot_CentFT0C", "hQ2yBTot_CentFT0C;centrality FT0C (%);Q_{2,y}^{BTot}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP2FT0M_CentFT0C", "2nd harmonics event plane FT0M;centrality FT0C (%);#Psi_{2}^{FT0M} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0A_CentFT0C", "2nd harmonics event plane FT0A;centrality FT0C (%);#Psi_{2}^{FT0A} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0C_CentFT0C", "2nd harmonics event plane FT0C;centrality FT0C (%);#Psi_{2}^{FT0C} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FV0A_CentFT0C", "2nd harmonics event plane FV0A;centrality FT0C (%);#Psi_{2}^{FV0A} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BPos_CentFT0C", "2nd harmonics event plane BPos;centrality FT0C (%);#Psi_{2}^{BPos} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BNeg_CentFT0C", "2nd harmonics event plane BNeg;centrality FT0C (%);#Psi_{2}^{BNeg} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BTot_CentFT0C", "2nd harmonics event plane BTot;centrality FT0C (%);#Psi_{2}^{BTot} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hPrfQ2FT0MQ2BPos_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0MQ2BNeg_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2BPosQ2BNeg_CentFT0C", "Q_{2}^{BPos} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{BPos} #upoint Q_{2}^{BNeg}", kTProfile, {{100, 0, 100}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.
    fRegistry->add("Event/before/hPrfQ2FT0CQ2BPos_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0CQ2BNeg_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0CQ2BTot_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BTot};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2BPos_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2BNeg_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2BTot_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BTot};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2FT0C_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{FT0C};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{FT0C}", kTProfile, {{100, 0, 100}}, false); // this is necessary for dimuons
    fRegistry->add("Event/before/hPrfQ2FV0AQ2BPos_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FV0AQ2BNeg_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FV0AQ2BTot_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BTot};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ2FV0AQ2FT0C_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{FT0C};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{FT0C}", kTProfile, {{100, 0, 100}}, false); // this is necessary for dimuons
  } else if constexpr (nmod == 3) {
    fRegistry->add("Event/before/hQ3xFT0M_CentFT0C", "hQ3xFT0M_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0M}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0M_CentFT0C", "hQ3yFT0M_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0M}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xFT0A_CentFT0C", "hQ3xFT0A_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0A_CentFT0C", "hQ3yFT0A_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xFT0C_CentFT0C", "hQ3xFT0C_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0C}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0C_CentFT0C", "hQ3yFT0C_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0C}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xFV0A_CentFT0C", "hQ3xFV0A_CentFT0C;centrality FT0C (%);Q_{3,x}^{FV0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFV0A_CentFT0C", "hQ3yFV0A_CentFT0C;centrality FT0C (%);Q_{3,y}^{FV0A}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xBPos_CentFT0C", "hQ3xBPos_CentFT0C;centrality FT0C (%);Q_{3,x}^{BPos}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yBPos_CentFT0C", "hQ3yBPos_CentFT0C;centrality FT0C (%);Q_{3,y}^{BPos}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xBNeg_CentFT0C", "hQ3xBNeg_CentFT0C;centrality FT0C (%);Q_{3,x}^{BNeg}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yBNeg_CentFT0C", "hQ3yBNeg_CentFT0C;centrality FT0C (%);Q_{3,y}^{BNeg}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xBTot_CentFT0C", "hQ3xBTot_CentFT0C;centrality FT0C (%);Q_{3,x}^{BTot}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yBTot_CentFT0C", "hQ3yBTot_CentFT0C;centrality FT0C (%);Q_{3,y}^{BTot}", kTH2D, {{100, 0, 100}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP3FT0M_CentFT0C", "3rd harmonics event plane FT0M;centrality FT0C (%);#Psi_{3}^{FT0M} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3FT0A_CentFT0C", "3rd harmonics event plane FT0A;centrality FT0C (%);#Psi_{3}^{FT0A} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3FT0C_CentFT0C", "3rd harmonics event plane FT0C;centrality FT0C (%);#Psi_{3}^{FT0C} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3FV0A_CentFT0C", "3rd harmonics event plane FV0A;centrality FT0C (%);#Psi_{3}^{FV0A} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3BPos_CentFT0C", "3rd harmonics event plane BPos;centrality FT0C (%);#Psi_{3}^{BPos} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3BNeg_CentFT0C", "3rd harmonics event plane BNeg;centrality FT0C (%);#Psi_{3}^{BNeg} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3BTot_CentFT0C", "3rd harmonics event plane BTot;centrality FT0C (%);#Psi_{3}^{BTot} (rad.)", kTH2D, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hPrfQ3FT0MQ3BPos_CentFT0C", "Q_{3}^{FT0M} #upoint Q_{3}^{BPos};centrality FT0C (%);Q_{3}^{FT0M} #upoint Q_{3}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0MQ3BNeg_CentFT0C", "Q_{3}^{FT0M} #upoint Q_{3}^{BNeg};centrality FT0C (%);Q_{3}^{FT0M} #upoint Q_{3}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3BPosQ3BNeg_CentFT0C", "Q_{3}^{BPos} #upoint Q_{3}^{BNeg};centrality FT0C (%);Q_{3}^{BPos} #upoint Q_{3}^{BNeg}", kTProfile, {{100, 0, 100}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.
    fRegistry->add("Event/before/hPrfQ3FT0CQ3BPos_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{BPos};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0CQ3BNeg_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{BNeg};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0CQ3BTot_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{BTot};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3BPos_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{BPos};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3BNeg_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{BNeg};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3BTot_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{BTot};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3FT0C_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{FT0C};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{FT0C}", kTProfile, {{100, 0, 100}}, false); // this is necessary for dimuons
    fRegistry->add("Event/before/hPrfQ3FV0AQ3BPos_CentFT0C", "Q_{3}^{FV0A} #upoint Q_{3}^{BPos};centrality FT0C (%);Q_{3}^{FV0A} #upoint Q_{3}^{BPos}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FV0AQ3BNeg_CentFT0C", "Q_{3}^{FV0A} #upoint Q_{3}^{BNeg};centrality FT0C (%);Q_{3}^{FV0A} #upoint Q_{3}^{BNeg}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FV0AQ3BTot_CentFT0C", "Q_{3}^{FV0A} #upoint Q_{3}^{BTot};centrality FT0C (%);Q_{3}^{FV0A} #upoint Q_{3}^{BTot}", kTProfile, {{100, 0, 100}}, false);
    fRegistry->add("Event/before/hPrfQ3FV0AQ3FT0C_CentFT0C", "Q_{3}^{FV0A} #upoint Q_{3}^{FT0C};centrality FT0C (%);Q_{3}^{FV0A} #upoint Q_{3}^{FT0C}", kTProfile, {{100, 0, 100}}, false); // this is necessary for dimuons
  }
  fRegistry->addClone("Event/before/", "Event/after/");
}

template <const int ev_id, const int nmod = -1, typename TCollision>
void fillEventInfo(HistogramRegistry* fRegistry, TCollision const& collision, const float /*weight*/ = 1.f)
{
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
  }
  if (collision.sel8()) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 10.0);
  }
  if (abs(collision.posZ()) < 10.0) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 11.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 12.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 13.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 14.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 15.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 16.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 17.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 18.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 19.0);
  }
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hNTracksPVvsOccupancy"), collision.multNTracksPV(), collision.trackOccupancyInTimeRange());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCorrOccupancy"), collision.ft0cOccupancyInTimeRange(), collision.trackOccupancyInTimeRange());

  if constexpr (nmod == 2) {
    std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
    std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
    std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
    std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};
    std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
    std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
    std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0M_CentFT0C"), collision.centFT0C(), collision.q2xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0M_CentFT0C"), collision.centFT0C(), collision.q2yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0A_CentFT0C"), collision.centFT0C(), collision.q2xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0A_CentFT0C"), collision.centFT0C(), collision.q2yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0C_CentFT0C"), collision.centFT0C(), collision.q2xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0C_CentFT0C"), collision.centFT0C(), collision.q2yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFV0A_CentFT0C"), collision.centFT0C(), collision.q2xfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFV0A_CentFT0C"), collision.centFT0C(), collision.q2yfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBPos_CentFT0C"), collision.centFT0C(), collision.q2xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBPos_CentFT0C"), collision.centFT0C(), collision.q2ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBNeg_CentFT0C"), collision.centFT0C(), collision.q2xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBNeg_CentFT0C"), collision.centFT0C(), collision.q2ybneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBTot_CentFT0C"), collision.centFT0C(), collision.q2xbtot());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBTot_CentFT0C"), collision.centFT0C(), collision.q2ybtot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0M_CentFT0C"), collision.centFT0C(), collision.ep2ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0A_CentFT0C"), collision.centFT0C(), collision.ep2ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0C_CentFT0C"), collision.centFT0C(), collision.ep2ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FV0A_CentFT0C"), collision.centFT0C(), collision.ep2fv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BPos_CentFT0C"), collision.centFT0C(), collision.ep2bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BNeg_CentFT0C"), collision.centFT0C(), collision.ep2bneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BTot_CentFT0C"), collision.centFT0C(), collision.ep2btot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0MQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0MQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2BPosQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2bpos, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2ft0c));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FV0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FV0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FV0AQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FV0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2ft0c));
  } else if constexpr (nmod == 3) {
    std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
    std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
    std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
    std::array<float, 2> q3fv0a = {collision.q3xfv0a(), collision.q3yfv0a()};
    std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
    std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};
    std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0M_CentFT0C"), collision.centFT0C(), collision.q3xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0M_CentFT0C"), collision.centFT0C(), collision.q3yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0A_CentFT0C"), collision.centFT0C(), collision.q3xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0A_CentFT0C"), collision.centFT0C(), collision.q3yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0C_CentFT0C"), collision.centFT0C(), collision.q3xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0C_CentFT0C"), collision.centFT0C(), collision.q3yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFV0A_CentFT0C"), collision.centFT0C(), collision.q3xfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFV0A_CentFT0C"), collision.centFT0C(), collision.q3yfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xBPos_CentFT0C"), collision.centFT0C(), collision.q3xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yBPos_CentFT0C"), collision.centFT0C(), collision.q3ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xBNeg_CentFT0C"), collision.centFT0C(), collision.q3xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yBNeg_CentFT0C"), collision.centFT0C(), collision.q3ybneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xBTot_CentFT0C"), collision.centFT0C(), collision.q3xbtot());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yBTot_CentFT0C"), collision.centFT0C(), collision.q3ybtot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0M_CentFT0C"), collision.centFT0C(), collision.ep3ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0A_CentFT0C"), collision.centFT0C(), collision.ep3ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0C_CentFT0C"), collision.centFT0C(), collision.ep3ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FV0A_CentFT0C"), collision.centFT0C(), collision.ep3ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3BPos_CentFT0C"), collision.centFT0C(), collision.ep3bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3BNeg_CentFT0C"), collision.centFT0C(), collision.ep3bneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3BTot_CentFT0C"), collision.centFT0C(), collision.ep3btot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0MQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0m, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0MQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0m, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3BPosQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3bpos, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3ft0c));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FV0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3fv0a, q3ft0c));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FV0AQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3fv0a, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FV0AQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3fv0a, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FV0AQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3fv0a, q3btot));
  }
}
} // namespace o2::aod::pwgem::dilepton::utils::eventhistogram
#endif // PWGEM_DILEPTON_UTILS_EVENTHISTOGRAMS_H_
