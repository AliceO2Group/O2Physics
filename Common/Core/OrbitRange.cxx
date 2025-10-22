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
// Container to store minimum and maximum orbit counter
//

#include "Common/Core/OrbitRange.h"

#include <TCollection.h>
#include <TIterator.h>
#include <TMathBase.h>
#include <TObject.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdio>

ClassImp(OrbitRange)

  Long64_t OrbitRange::Merge(TCollection* list)
{
  // Merge a list of OrbitRange objects
  // Stores minimum and maximum orbit among all merged orbit ranges
  // Returns the number of merged objects (including this).

  if (!list) {
    return 0;
  }

  if (list->IsEmpty()) {
    return 1;
  }

  Long64_t count = 0;
  TIterator* iter = list->MakeIterator();
  TObject* obj = nullptr;
  while ((obj = iter->Next())) {
    OrbitRange* entry = dynamic_cast<OrbitRange*>(obj);
    if (entry == nullptr) {
      continue;
    }
    if (fRunNumber != entry->GetRunNumber()) {
      printf("Warning: merging of orbit ranges for different runs is forbidden\n");
      continue;
    }
    fMinOrbit = TMath::Min(fMinOrbit, entry->GetMinOrbit());
    fMaxOrbit = TMath::Max(fMaxOrbit, entry->GetMaxOrbit());
    count++;
  }

  return count + 1;
}
