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
using namespace o2::framework;

namespace o2::aod::pwgem::dilepton::utils::eventhistogram
{
const int nbin_ev = 13;
template <const int nmod = -1>
void addEventHistograms(HistogramRegistry* fRegistry)
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{nbin_ev, 0.5, nbin_ev + 0.5}}, false);
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
  hCollisionCounter->GetXaxis()->SetBinLabel(12, "Calibrated Q vector");
  hCollisionCounter->GetXaxis()->SetBinLabel(13, "accepted");

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0A_HMpp", "hCentFT0A for HM pp;centrality FT0A (%)", kTH1F, {{100, 0, 1}}, false);
  fRegistry->add("Event/before/hCentFT0C_HMpp", "hCentFT0C for HM pp;centrality FT0C (%)", kTH1F, {{100, 0, 1}}, false);
  fRegistry->add("Event/before/hCentFT0M_HMpp", "hCentFT0M for HM pp;centrality FT0M (%)", kTH1F, {{100, 0, 1}}, false);
  fRegistry->add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", kTH2F, {{60, 0, 60000}, {2000, 0, 20000}}, false);
  fRegistry->add("Event/before/hSpherocity", "hSpherocity;spherocity", kTH1F, {{100, 0, 1}}, false);

  if constexpr (nmod == 2) { // Q2
    fRegistry->add("Event/before/hQ2xFT0M_CentFT0C", "hQ2xFT0M_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0M_CentFT0C", "hQ2yFT0M_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0A_CentFT0C", "hQ2xFT0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0A_CentFT0C", "hQ2yFT0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0C_CentFT0C", "hQ2xFT0C_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0C_CentFT0C", "hQ2yFT0C_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xTPCpos_CentFT0C", "hQ2xTPCpos_CentFT0C;centrality FT0C (%);Q_{2,x}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yTPCpos_CentFT0C", "hQ2yTPCpos_CentFT0C;centrality FT0C (%);Q_{2,y}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xTPCneg_CentFT0C", "hQ2xTPCneg_CentFT0C;centrality FT0C (%);Q_{2,x}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yTPCneg_CentFT0C", "hQ2yTPCneg_CentFT0C;centrality FT0C (%);Q_{2,y}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xTPCall_CentFT0C", "hQ2xTPCall_CentFT0C;centrality FT0C (%);Q_{2,x}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yTPCall_CentFT0C", "hQ2yTPCall_CentFT0C;centrality FT0C (%);Q_{2,y}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP2FT0M_CentFT0C", "2nd harmonics event plane FT0M;centrality FT0C (%);#Psi_{2}^{FT0M} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0A_CentFT0C", "2nd harmonics event plane FT0A;centrality FT0C (%);#Psi_{2}^{FT0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0C_CentFT0C", "2nd harmonics event plane FT0C;centrality FT0C (%);#Psi_{2}^{FT0C} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2TPCpos_CentFT0C", "2nd harmonics event plane TPCpos;centrality FT0C (%);#Psi_{2}^{TPCpos} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2TPCneg_CentFT0C", "2nd harmonics event plane TPCneg;centrality FT0C (%);#Psi_{2}^{TPCneg} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2TPCall_CentFT0C", "2nd harmonics event plane TPCall;centrality FT0C (%);#Psi_{2}^{TPCall} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hPrfQ2FT0MQ2TPCpos_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{TPCpos};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0MQ2TPCneg_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{TPCneg};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2TPCposQ2TPCneg_CentFT0C", "Q_{2}^{TPCpos} #upoint Q_{2}^{TPCneg};centrality FT0C (%);Q_{2}^{TPCpos} #upoint Q_{2}^{TPCneg}", kTProfile, {{110, 0, 110}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.
    fRegistry->add("Event/before/hPrfQ2FT0CQ2TPCpos_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{TPCpos};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0CQ2TPCneg_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{TPCneg};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0CQ2TPCall_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{TPCall};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2TPCpos_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{TPCpos};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2TPCneg_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{TPCneg};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2TPCall_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{TPCall};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ2FT0AQ2FT0C_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{FT0C};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{FT0C}", kTProfile, {{110, 0, 110}}, false); // this is necessary for dimuons
  } else if constexpr (nmod == 3) {                                                                                                                                                         // Q3
    fRegistry->add("Event/before/hQ3xFT0M_CentFT0C", "hQ3xFT0M_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0M_CentFT0C", "hQ3yFT0M_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xFT0A_CentFT0C", "hQ3xFT0A_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0A_CentFT0C", "hQ3yFT0A_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xFT0C_CentFT0C", "hQ3xFT0C_CentFT0C;centrality FT0C (%);Q_{3,x}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yFT0C_CentFT0C", "hQ3yFT0C_CentFT0C;centrality FT0C (%);Q_{3,y}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xTPCpos_CentFT0C", "hQ3xTPCpos_CentFT0C;centrality FT0C (%);Q_{3,x}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yTPCpos_CentFT0C", "hQ3yTPCpos_CentFT0C;centrality FT0C (%);Q_{3,y}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xTPCneg_CentFT0C", "hQ3xTPCneg_CentFT0C;centrality FT0C (%);Q_{3,x}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yTPCneg_CentFT0C", "hQ3yTPCneg_CentFT0C;centrality FT0C (%);Q_{3,y}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3xTPCall_CentFT0C", "hQ3xTPCall_CentFT0C;centrality FT0C (%);Q_{3,x}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ3yTPCall_CentFT0C", "hQ3yTPCall_CentFT0C;centrality FT0C (%);Q_{3,y}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP3FT0M_CentFT0C", "3rd harmonics event plane FT0M;centrality FT0C (%);#Psi_{3}^{FT0M} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3FT0A_CentFT0C", "3rd harmonics event plane FT0A;centrality FT0C (%);#Psi_{3}^{FT0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3FT0C_CentFT0C", "3rd harmonics event plane FT0C;centrality FT0C (%);#Psi_{3}^{FT0C} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3TPCpos_CentFT0C", "3rd harmonics event plane TPCpos;centrality FT0C (%);#Psi_{3}^{TPCpos} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3TPCneg_CentFT0C", "3rd harmonics event plane TPCneg;centrality FT0C (%);#Psi_{3}^{TPCneg} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP3TPCall_CentFT0C", "3rd harmonics event plane TPCall;centrality FT0C (%);#Psi_{3}^{TPCall} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hPrfQ3FT0MQ3TPCpos_CentFT0C", "Q_{3}^{FT0M} #upoint Q_{3}^{TPCpos};centrality FT0C (%);Q_{3}^{FT0M} #upoint Q_{3}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0MQ3TPCneg_CentFT0C", "Q_{3}^{FT0M} #upoint Q_{3}^{TPCneg};centrality FT0C (%);Q_{3}^{FT0M} #upoint Q_{3}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3TPCposQ3TPCneg_CentFT0C", "Q_{3}^{TPCpos} #upoint Q_{3}^{TPCneg};centrality FT0C (%);Q_{3}^{TPCpos} #upoint Q_{3}^{TPCneg}", kTProfile, {{110, 0, 110}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.
    fRegistry->add("Event/before/hPrfQ3FT0CQ3TPCpos_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{TPCpos};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0CQ3TPCneg_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{TPCneg};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0CQ3TPCall_CentFT0C", "Q_{3}^{FT0C} #upoint Q_{3}^{TPCall};centrality FT0C (%);Q_{3}^{FT0C} #upoint Q_{3}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3TPCpos_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{TPCpos};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3TPCneg_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{TPCneg};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3TPCall_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{TPCall};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ3FT0AQ3FT0C_CentFT0C", "Q_{3}^{FT0A} #upoint Q_{3}^{FT0C};centrality FT0C (%);Q_{3}^{FT0A} #upoint Q_{3}^{FT0C}", kTProfile, {{110, 0, 110}}, false); // this is necessary for dimuons
  } else if constexpr (nmod == 4) {                                                                                                                                                         // Q4
    fRegistry->add("Event/before/hQ4xFT0M_CentFT0C", "hQ4xFT0M_CentFT0C;centrality FT0C (%);Q_{4,x}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yFT0M_CentFT0C", "hQ4yFT0M_CentFT0C;centrality FT0C (%);Q_{4,y}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4xFT0A_CentFT0C", "hQ4xFT0A_CentFT0C;centrality FT0C (%);Q_{4,x}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yFT0A_CentFT0C", "hQ4yFT0A_CentFT0C;centrality FT0C (%);Q_{4,y}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4xFT0C_CentFT0C", "hQ4xFT0C_CentFT0C;centrality FT0C (%);Q_{4,x}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yFT0C_CentFT0C", "hQ4yFT0C_CentFT0C;centrality FT0C (%);Q_{4,y}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4xTPCpos_CentFT0C", "hQ4xTPCpos_CentFT0C;centrality FT0C (%);Q_{4,x}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yTPCpos_CentFT0C", "hQ4yTPCpos_CentFT0C;centrality FT0C (%);Q_{4,y}^{TPCpos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4xTPCneg_CentFT0C", "hQ4xTPCneg_CentFT0C;centrality FT0C (%);Q_{4,x}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yTPCneg_CentFT0C", "hQ4yTPCneg_CentFT0C;centrality FT0C (%);Q_{4,y}^{TPCneg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4xTPCall_CentFT0C", "hQ4xTPCall_CentFT0C;centrality FT0C (%);Q_{4,x}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ4yTPCall_CentFT0C", "hQ4yTPCall_CentFT0C;centrality FT0C (%);Q_{4,y}^{TPCall}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP4FT0M_CentFT0C", "4rd harmonics event plane FT0M;centrality FT0C (%);#Psi_{4}^{FT0M} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP4FT0A_CentFT0C", "4rd harmonics event plane FT0A;centrality FT0C (%);#Psi_{4}^{FT0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP4FT0C_CentFT0C", "4rd harmonics event plane FT0C;centrality FT0C (%);#Psi_{4}^{FT0C} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP4TPCpos_CentFT0C", "4rd harmonics event plane TPCpos;centrality FT0C (%);#Psi_{4}^{TPCpos} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP4TPCneg_CentFT0C", "4rd harmonics event plane TPCneg;centrality FT0C (%);#Psi_{4}^{TPCneg} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP4TPCall_CentFT0C", "4rd harmonics event plane TPCall;centrality FT0C (%);#Psi_{4}^{TPCall} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hPrfQ4FT0MQ4TPCpos_CentFT0C", "Q_{4}^{FT0M} #upoint Q_{4}^{TPCpos};centrality FT0C (%);Q_{4}^{FT0M} #upoint Q_{4}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0MQ4TPCneg_CentFT0C", "Q_{4}^{FT0M} #upoint Q_{4}^{TPCneg};centrality FT0C (%);Q_{4}^{FT0M} #upoint Q_{4}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4TPCposQ4TPCneg_CentFT0C", "Q_{4}^{TPCpos} #upoint Q_{4}^{TPCneg};centrality FT0C (%);Q_{4}^{TPCpos} #upoint Q_{4}^{TPCneg}", kTProfile, {{110, 0, 110}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.
    fRegistry->add("Event/before/hPrfQ4FT0CQ4TPCpos_CentFT0C", "Q_{4}^{FT0C} #upoint Q_{4}^{TPCpos};centrality FT0C (%);Q_{4}^{FT0C} #upoint Q_{4}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0CQ4TPCneg_CentFT0C", "Q_{4}^{FT0C} #upoint Q_{4}^{TPCneg};centrality FT0C (%);Q_{4}^{FT0C} #upoint Q_{4}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0CQ4TPCall_CentFT0C", "Q_{4}^{FT0C} #upoint Q_{4}^{TPCall};centrality FT0C (%);Q_{4}^{FT0C} #upoint Q_{4}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0AQ4TPCpos_CentFT0C", "Q_{4}^{FT0A} #upoint Q_{4}^{TPCpos};centrality FT0C (%);Q_{4}^{FT0A} #upoint Q_{4}^{TPCpos}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0AQ4TPCneg_CentFT0C", "Q_{4}^{FT0A} #upoint Q_{4}^{TPCneg};centrality FT0C (%);Q_{4}^{FT0A} #upoint Q_{4}^{TPCneg}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0AQ4TPCall_CentFT0C", "Q_{4}^{FT0A} #upoint Q_{4}^{TPCall};centrality FT0C (%);Q_{4}^{FT0A} #upoint Q_{4}^{TPCall}", kTProfile, {{110, 0, 110}}, false);
    fRegistry->add("Event/before/hPrfQ4FT0AQ4FT0C_CentFT0C", "Q_{4}^{FT0A} #upoint Q_{4}^{FT0C};centrality FT0C (%);Q_{4}^{FT0A} #upoint Q_{4}^{FT0C}", kTProfile, {{110, 0, 110}}, false); // this is necessary for dimuons
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
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A_HMpp"), collision.centFT0A());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C_HMpp"), collision.centFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M_HMpp"), collision.centFT0M());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange());

  if constexpr (nmod == 2) { // Q2
    std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
    std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
    std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
    std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
    std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
    std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0M_CentFT0C"), collision.centFT0C(), collision.q2xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0M_CentFT0C"), collision.centFT0C(), collision.q2yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0A_CentFT0C"), collision.centFT0C(), collision.q2xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0A_CentFT0C"), collision.centFT0C(), collision.q2yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0C_CentFT0C"), collision.centFT0C(), collision.q2xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0C_CentFT0C"), collision.centFT0C(), collision.q2yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xTPCpos_CentFT0C"), collision.centFT0C(), collision.q2xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yTPCpos_CentFT0C"), collision.centFT0C(), collision.q2ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xTPCneg_CentFT0C"), collision.centFT0C(), collision.q2xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yTPCneg_CentFT0C"), collision.centFT0C(), collision.q2ybneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xTPCall_CentFT0C"), collision.centFT0C(), collision.q2xbtot());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yTPCall_CentFT0C"), collision.centFT0C(), collision.q2ybtot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0M_CentFT0C"), collision.centFT0C(), collision.ep2ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0A_CentFT0C"), collision.centFT0C(), collision.ep2ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0C_CentFT0C"), collision.centFT0C(), collision.ep2ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2TPCpos_CentFT0C"), collision.centFT0C(), collision.ep2bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2TPCneg_CentFT0C"), collision.centFT0C(), collision.ep2bneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2TPCall_CentFT0C"), collision.centFT0C(), collision.ep2btot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0MQ2TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0MQ2TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2TPCposQ2TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2bpos, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0CQ2TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ2FT0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2ft0c));
  } else if constexpr (nmod == 3) { // Q3
    std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
    std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
    std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
    std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
    std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};
    std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0M_CentFT0C"), collision.centFT0C(), collision.q3xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0M_CentFT0C"), collision.centFT0C(), collision.q3yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0A_CentFT0C"), collision.centFT0C(), collision.q3xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0A_CentFT0C"), collision.centFT0C(), collision.q3yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xFT0C_CentFT0C"), collision.centFT0C(), collision.q3xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yFT0C_CentFT0C"), collision.centFT0C(), collision.q3yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xTPCpos_CentFT0C"), collision.centFT0C(), collision.q3xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yTPCpos_CentFT0C"), collision.centFT0C(), collision.q3ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xTPCneg_CentFT0C"), collision.centFT0C(), collision.q3xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yTPCneg_CentFT0C"), collision.centFT0C(), collision.q3ybneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3xTPCall_CentFT0C"), collision.centFT0C(), collision.q3xbtot());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ3yTPCall_CentFT0C"), collision.centFT0C(), collision.q3ybtot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0M_CentFT0C"), collision.centFT0C(), collision.ep3ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0A_CentFT0C"), collision.centFT0C(), collision.ep3ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3FT0C_CentFT0C"), collision.centFT0C(), collision.ep3ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3TPCpos_CentFT0C"), collision.centFT0C(), collision.ep3bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3TPCneg_CentFT0C"), collision.centFT0C(), collision.ep3bneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP3TPCall_CentFT0C"), collision.centFT0C(), collision.ep3btot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0MQ3TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0m, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0MQ3TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0m, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3TPCposQ3TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3bpos, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0CQ3TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0c, q3btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ3FT0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q3ft0a, q3ft0c));
  } else if constexpr (nmod == 4) { // Q4
    std::array<float, 2> q4ft0m = {collision.q4xft0m(), collision.q4yft0m()};
    std::array<float, 2> q4ft0a = {collision.q4xft0a(), collision.q4yft0a()};
    std::array<float, 2> q4ft0c = {collision.q4xft0c(), collision.q4yft0c()};
    std::array<float, 2> q4bpos = {collision.q4xbpos(), collision.q4ybpos()};
    std::array<float, 2> q4bneg = {collision.q4xbneg(), collision.q4ybneg()};
    std::array<float, 2> q4btot = {collision.q4xbtot(), collision.q4ybtot()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xFT0M_CentFT0C"), collision.centFT0C(), collision.q4xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yFT0M_CentFT0C"), collision.centFT0C(), collision.q4yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xFT0A_CentFT0C"), collision.centFT0C(), collision.q4xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yFT0A_CentFT0C"), collision.centFT0C(), collision.q4yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xFT0C_CentFT0C"), collision.centFT0C(), collision.q4xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yFT0C_CentFT0C"), collision.centFT0C(), collision.q4yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xTPCpos_CentFT0C"), collision.centFT0C(), collision.q4xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yTPCpos_CentFT0C"), collision.centFT0C(), collision.q4ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xTPCneg_CentFT0C"), collision.centFT0C(), collision.q4xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yTPCneg_CentFT0C"), collision.centFT0C(), collision.q4ybneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4xTPCall_CentFT0C"), collision.centFT0C(), collision.q4xbtot());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ4yTPCall_CentFT0C"), collision.centFT0C(), collision.q4ybtot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4FT0M_CentFT0C"), collision.centFT0C(), collision.ep4ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4FT0A_CentFT0C"), collision.centFT0C(), collision.ep4ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4FT0C_CentFT0C"), collision.centFT0C(), collision.ep4ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4TPCpos_CentFT0C"), collision.centFT0C(), collision.ep4bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4TPCneg_CentFT0C"), collision.centFT0C(), collision.ep4bneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP4TPCall_CentFT0C"), collision.centFT0C(), collision.ep4btot());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0MQ4TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0m, q4bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0MQ4TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0m, q4bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4TPCposQ4TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4bpos, q4bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0AQ4TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0a, q4bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0AQ4TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0a, q4bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0AQ4TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0a, q4btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0CQ4TPCpos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0c, q4bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0CQ4TPCneg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0c, q4bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0CQ4TPCall_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0c, q4btot));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hPrfQ4FT0AQ4FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q4ft0a, q4ft0c));
  }
}

} // namespace o2::aod::pwgem::dilepton::utils::eventhistogram

#endif // PWGEM_DILEPTON_UTILS_EVENTHISTOGRAMS_H_
