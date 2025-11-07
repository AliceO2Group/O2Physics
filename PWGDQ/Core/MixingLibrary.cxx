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
#include "PWGDQ/Core/MixingLibrary.h"

void o2::aod::dqmixing::SetUpMixing(MixingHandler* mh, const char* mixingVarible)
{
  std::string nameStr = mixingVarible;
  if (!nameStr.compare("Centrality1")) {
    std::vector<float> fCentLimsHashing = {0.0f, 20.0f, 40.0f, 60.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("Centrality2")) {
    std::vector<float> fCentLimsHashing = {0.0f, 10.0f, 20.0f, 40.0f, 60.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("Centrality3")) {
    std::vector<float> fCentLimsHashing = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("Centrality4")) {
    std::vector<float> fCentLimsHashing = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("Centrality5")) {
    std::vector<float> fCentLimsHashing = {0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("Centrality6")) {
    std::vector<float> fCentLimsHashing = {0.0f, 2.5f, 5.0f, 7.5f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentVZERO, fCentLimsHashing.size(), fCentLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C1")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 20.0f, 40.0f, 60.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C2")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 10.0f, 20.0f, 40.0f, 60.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C3")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C4")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C5")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("CentralityFT0C6")) {
    std::vector<float> fCentFT0CLimsHashing = {0.0f, 2.5f, 5.0f, 7.5f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f};
    mh->AddMixingVariable(VarManager::kCentFT0C, fCentFT0CLimsHashing.size(), fCentFT0CLimsHashing);
  }
  if (!nameStr.compare("Mult1")) {
    std::vector<float> fMultLimsHashing = {0.0f, 10.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 120.0f, 160.0f, 350.0f};
    mh->AddMixingVariable(VarManager::kVtxNcontrib, fMultLimsHashing.size(), fMultLimsHashing);
  }
  if (!nameStr.compare("Mult2")) {
    std::vector<float> fMultLimsHashing = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 120.0f, 140.0f, 180.0f, 350.0f};
    mh->AddMixingVariable(VarManager::kVtxNcontrib, fMultLimsHashing.size(), fMultLimsHashing);
  }
  if (!nameStr.compare("Mult3")) {
    std::vector<float> fMultLimsHashing = {0.0f, 5.0f, 10.0f, 15.0f, 20.0f, 25.0f, 30.0f, 35.0f, 40.0f, 45.0f, 50.0f, 55.0f, 60.0f, 65.0f, 70.0f, 75.0f, 80.0f, 85.0f, 90.0f, 100.0f, 120.0f, 140.0f, 165.0f, 200.0f, 350.0f};
    mh->AddMixingVariable(VarManager::kVtxNcontrib, fMultLimsHashing.size(), fMultLimsHashing);
  }
  if (!nameStr.compare("Vtx1")) {
    std::vector<float> fZLimsHashing = {-10.0f, 0.0f, 10.0f};
    mh->AddMixingVariable(VarManager::kVtxZ, fZLimsHashing.size(), fZLimsHashing);
  }
  if (!nameStr.compare("Vtx2")) {
    std::vector<float> fZLimsHashing = {-10.0f, -5.0f, 0.0f, 5.0f, 10.0f};
    mh->AddMixingVariable(VarManager::kVtxZ, fZLimsHashing.size(), fZLimsHashing);
  }
  if (!nameStr.compare("Vtx3")) {
    std::vector<float> fZLimsHashing = {-10.0f, -7.5f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.5f, 10.0f};
    mh->AddMixingVariable(VarManager::kVtxZ, fZLimsHashing.size(), fZLimsHashing);
  }
  if (!nameStr.compare("Vtx4")) {
    std::vector<float> fZLimsHashing = {-10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f};
    mh->AddMixingVariable(VarManager::kVtxZ, fZLimsHashing.size(), fZLimsHashing);
  }
  if (!nameStr.compare("Vtx5")) {
    std::vector<float> fZLimsHashing = {-10.0f, -9.0f, -8.0f, -7.0f, -6.0f, -5.0f, -4.0f, -3.0f, -2.0f, -1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f};
    mh->AddMixingVariable(VarManager::kVtxZ, fZLimsHashing.size(), fZLimsHashing);
  }
  if (!nameStr.compare("Occupancy1")) {
    std::vector<float> fOccLimsHashing = {0.0f, 500.0f, 1000.0f, 2000.0f, 3000.0f, 6000.0f, 50000.0f};
    mh->AddMixingVariable(VarManager::kTrackOccupancyInTimeRange, fOccLimsHashing.size(), fOccLimsHashing);
  }
  if (!nameStr.compare("Occupancy2")) {
    std::vector<float> fOccLimsHashing = {0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 6000.0f, 50000.0f};
    mh->AddMixingVariable(VarManager::kTrackOccupancyInTimeRange, fOccLimsHashing.size(), fOccLimsHashing);
  }
  if (!nameStr.compare("Occupancy3")) {
    std::vector<float> fOccLimsHashing = {0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f};
    mh->AddMixingVariable(VarManager::kTrackOccupancyInTimeRange, fOccLimsHashing.size(), fOccLimsHashing);
  }
  if (!nameStr.compare("Psi2A1")) {
    std::vector<float> fPsi2A = {-TMath::Pi() / 2., 0.0f, TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2A, fPsi2A.size(), fPsi2A);
  }
  if (!nameStr.compare("Psi2A2")) {
    std::vector<float> fPsi2A = {-TMath::Pi() / 2., -TMath::Pi() / 4., 0.0f, TMath::Pi() / 4., TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2A, fPsi2A.size(), fPsi2A);
  }
  if (!nameStr.compare("Psi2A3")) {
    std::vector<float> fPsi2A = {-4 * TMath::Pi() / 8., -3 * TMath::Pi() / 8., -2 * TMath::Pi() / 8., -TMath::Pi() / 8., 0.0f, TMath::Pi() / 8., 2 * TMath::Pi() / 8., 3 * TMath::Pi() / 8., 4 * TMath::Pi() / 8.};
    mh->AddMixingVariable(VarManager::kPsi2A, fPsi2A.size(), fPsi2A);
  }
  if (!nameStr.compare("Psi2A4")) {
    std::vector<float> fPsi2A = {-8 * TMath::Pi() / 16., -7 * TMath::Pi() / 16., -6 * TMath::Pi() / 16., -5 * TMath::Pi() / 16., -4 * TMath::Pi() / 16., -3 * TMath::Pi() / 16., -2 * TMath::Pi() / 16., -TMath::Pi() / 16., 0.0f, TMath::Pi() / 16., 2 * TMath::Pi() / 16., 3 * TMath::Pi() / 16., 4 * TMath::Pi() / 16., 5 * TMath::Pi() / 16., 6 * TMath::Pi() / 16., 7 * TMath::Pi() / 16., 8 * TMath::Pi() / 16.};
    mh->AddMixingVariable(VarManager::kPsi2A, fPsi2A.size(), fPsi2A);
  }
  if (!nameStr.compare("Psi2A5")) {
    std::vector<float> fPsi2A = {-12 * TMath::Pi() / 24., -11 * TMath::Pi() / 24., -10 * TMath::Pi() / 24., -9 * TMath::Pi() / 24., -8 * TMath::Pi() / 24., -7 * TMath::Pi() / 24., -6 * TMath::Pi() / 24., -5 * TMath::Pi() / 24., -4 * TMath::Pi() / 24., -3 * TMath::Pi() / 24., -2 * TMath::Pi() / 24., -TMath::Pi() / 24., 0.0f, TMath::Pi() / 24., 2 * TMath::Pi() / 24., 3 * TMath::Pi() / 24., 4 * TMath::Pi() / 24., 5 * TMath::Pi() / 24., 6 * TMath::Pi() / 24., 7 * TMath::Pi() / 24., 8 * TMath::Pi() / 24., 9 * TMath::Pi() / 24., 10 * TMath::Pi() / 24., 11 * TMath::Pi() / 24., 12 * TMath::Pi() / 24.};
    mh->AddMixingVariable(VarManager::kPsi2A, fPsi2A.size(), fPsi2A);
  }
  if (!nameStr.compare("Psi2B1")) {
    std::vector<float> fPsi2B = {-TMath::Pi() / 2., 0.0f, TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2B, fPsi2B.size(), fPsi2B);
  }
  if (!nameStr.compare("Psi2B2")) {
    std::vector<float> fPsi2B = {-TMath::Pi() / 2., -TMath::Pi() / 4., 0.0f, TMath::Pi() / 4., TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2B, fPsi2B.size(), fPsi2B);
  }
  if (!nameStr.compare("Psi2B3")) {
    std::vector<float> fPsi2B = {-4 * TMath::Pi() / 8., -3 * TMath::Pi() / 8., -2 * TMath::Pi() / 8., -TMath::Pi() / 8., 0.0f, TMath::Pi() / 8., 2 * TMath::Pi() / 8., 3 * TMath::Pi() / 8., 4 * TMath::Pi() / 8.};
    mh->AddMixingVariable(VarManager::kPsi2B, fPsi2B.size(), fPsi2B);
  }
  if (!nameStr.compare("Psi2B4")) {
    std::vector<float> fPsi2B = {-8 * TMath::Pi() / 16., -7 * TMath::Pi() / 16., -6 * TMath::Pi() / 16., -5 * TMath::Pi() / 16., -4 * TMath::Pi() / 16., -3 * TMath::Pi() / 16., -2 * TMath::Pi() / 16., -TMath::Pi() / 16., 0.0f, TMath::Pi() / 16., 2 * TMath::Pi() / 16., 3 * TMath::Pi() / 16., 4 * TMath::Pi() / 16., 5 * TMath::Pi() / 16., 6 * TMath::Pi() / 16., 7 * TMath::Pi() / 16., 8 * TMath::Pi() / 16.};
    mh->AddMixingVariable(VarManager::kPsi2B, fPsi2B.size(), fPsi2B);
  }
  if (!nameStr.compare("Psi2B5")) {
    std::vector<float> fPsi2B = {-12 * TMath::Pi() / 24., -11 * TMath::Pi() / 24., -10 * TMath::Pi() / 24., -9 * TMath::Pi() / 24., -8 * TMath::Pi() / 24., -7 * TMath::Pi() / 24., -6 * TMath::Pi() / 24., -5 * TMath::Pi() / 24., -4 * TMath::Pi() / 24., -3 * TMath::Pi() / 24., -2 * TMath::Pi() / 24., -TMath::Pi() / 24., 0.0f, TMath::Pi() / 24., 2 * TMath::Pi() / 24., 3 * TMath::Pi() / 24., 4 * TMath::Pi() / 24., 5 * TMath::Pi() / 24., 6 * TMath::Pi() / 24., 7 * TMath::Pi() / 24., 8 * TMath::Pi() / 24., 9 * TMath::Pi() / 24., 10 * TMath::Pi() / 24., 11 * TMath::Pi() / 24., 12 * TMath::Pi() / 24.};
    mh->AddMixingVariable(VarManager::kPsi2B, fPsi2B.size(), fPsi2B);
  }
  if (!nameStr.compare("Psi2C1")) {
    std::vector<float> fPsi2C = {-TMath::Pi() / 2., 0.0f, TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2C, fPsi2C.size(), fPsi2C);
  }
  if (!nameStr.compare("Psi2C2")) {
    std::vector<float> fPsi2C = {-TMath::Pi() / 2., -TMath::Pi() / 4., 0.0f, TMath::Pi() / 4., TMath::Pi() / 2.};
    mh->AddMixingVariable(VarManager::kPsi2C, fPsi2C.size(), fPsi2C);
  }
  if (!nameStr.compare("Psi2C3")) {
    std::vector<float> fPsi2C = {-4 * TMath::Pi() / 8., -3 * TMath::Pi() / 8., -2 * TMath::Pi() / 8., -TMath::Pi() / 8., 0.0f, TMath::Pi() / 8., 2 * TMath::Pi() / 8., 3 * TMath::Pi() / 8., 4 * TMath::Pi() / 8.};
    mh->AddMixingVariable(VarManager::kPsi2C, fPsi2C.size(), fPsi2C);
  }
  if (!nameStr.compare("Psi2C4")) {
    std::vector<float> fPsi2C = {-8 * TMath::Pi() / 16., -7 * TMath::Pi() / 16., -6 * TMath::Pi() / 16., -5 * TMath::Pi() / 16., -4 * TMath::Pi() / 16., -3 * TMath::Pi() / 16., -2 * TMath::Pi() / 16., -TMath::Pi() / 16., 0.0f, TMath::Pi() / 16., 2 * TMath::Pi() / 16., 3 * TMath::Pi() / 16., 4 * TMath::Pi() / 16., 5 * TMath::Pi() / 16., 6 * TMath::Pi() / 16., 7 * TMath::Pi() / 16., 8 * TMath::Pi() / 16.};
    mh->AddMixingVariable(VarManager::kPsi2C, fPsi2C.size(), fPsi2C);
  }
  if (!nameStr.compare("Psi2C5")) {
    std::vector<float> fPsi2C = {-12 * TMath::Pi() / 24., -11 * TMath::Pi() / 24., -10 * TMath::Pi() / 24., -9 * TMath::Pi() / 24., -8 * TMath::Pi() / 24., -7 * TMath::Pi() / 24., -6 * TMath::Pi() / 24., -5 * TMath::Pi() / 24., -4 * TMath::Pi() / 24., -3 * TMath::Pi() / 24., -2 * TMath::Pi() / 24., -TMath::Pi() / 24., 0.0f, TMath::Pi() / 24., 2 * TMath::Pi() / 24., 3 * TMath::Pi() / 24., 4 * TMath::Pi() / 24., 5 * TMath::Pi() / 24., 6 * TMath::Pi() / 24., 7 * TMath::Pi() / 24., 8 * TMath::Pi() / 24., 9 * TMath::Pi() / 24., 10 * TMath::Pi() / 24., 11 * TMath::Pi() / 24., 12 * TMath::Pi() / 24.};
    mh->AddMixingVariable(VarManager::kPsi2C, fPsi2C.size(), fPsi2C);
  }
  if (!nameStr.compare("MedianTimeA1")) {
    std::vector<float> fMTLimsHashing = {-100.0f, -40.0f, -20.0f, 20.0f, 40.0f, 100.0f};
    mh->AddMixingVariable(VarManager::kNTPCmedianTimeLongA, fMTLimsHashing.size(), fMTLimsHashing);
  }
  if (!nameStr.compare("MedianTimeA2")) {
    std::vector<float> fMTLimsHashing = {-100.0f, -80.0f, -60.0f, -40.0f, -20.0f, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f};
    mh->AddMixingVariable(VarManager::kNTPCmedianTimeLongA, fMTLimsHashing.size(), fMTLimsHashing);
  }
  if (!nameStr.compare("MedianTimeA3")) {
    std::vector<float> fMTLimsHashing = {-100.0f, -80.0f, -60.0f, -40.0f, -30.0f, -20.0f, -10.0f, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 60.0f, 80.0f, 100.0f};
    mh->AddMixingVariable(VarManager::kNTPCmedianTimeLongA, fMTLimsHashing.size(), fMTLimsHashing);
  }
  if (!nameStr.compare("PileUpA1")) {
    std::vector<float> fPileUpLimsHashing = {0.0f, 1000.0f, 2000.0f, 6000.0f, 10000.0f, 20000.0f};
    mh->AddMixingVariable(VarManager::kNTPCcontribLongA, fPileUpLimsHashing.size(), fPileUpLimsHashing);
  }
  if (!nameStr.compare("PileUpA2")) {
    std::vector<float> fPileUpLimsHashing = {0.0f, 1000.0f, 2000.0f, 4000.0f, 6000.0f, 8000.0f, 10000.0f, 20000.0f};
    mh->AddMixingVariable(VarManager::kNTPCcontribLongA, fPileUpLimsHashing.size(), fPileUpLimsHashing);
  }
  if (!nameStr.compare("PileUpA3")) {
    std::vector<float> fPileUpLimsHashing = {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f, 5000.0f, 6000.0f, 8000.0f, 10000.0f, 20000.0f};
    mh->AddMixingVariable(VarManager::kNTPCcontribLongA, fPileUpLimsHashing.size(), fPileUpLimsHashing);
  }
  if (!nameStr.compare("PileUpA4")) {
    std::vector<float> fPileUpLimsHashing = {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f, 3500.0f, 4000.0f, 4500.0f, 5000.0f, 5500.0f, 6000.0f, 8000.0f, 10000.0f, 20000.0f};
    mh->AddMixingVariable(VarManager::kNTPCcontribLongA, fPileUpLimsHashing.size(), fPileUpLimsHashing);
  }
}
