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
}
