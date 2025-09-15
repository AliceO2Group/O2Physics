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

#ifndef COMMON_CORE_ORBITRANGE_H_
#define COMMON_CORE_ORBITRANGE_H_

#include <TNamed.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdint>

class TCollection;

class OrbitRange : public TNamed
{
 public:
  explicit OrbitRange(const char* name = "orbitRange") : TNamed(name, name), fRunNumber(0), fMinOrbit(0xFFFFFFFF), fMaxOrbit(0) {}
  ~OrbitRange() {}
  void SetRunNumber(uint32_t runNumber) { fRunNumber = runNumber; }
  void SetMinOrbit(uint32_t orbit) { fMinOrbit = orbit; }
  void SetMaxOrbit(uint32_t orbit) { fMaxOrbit = orbit; }
  uint32_t GetRunNumber() { return fRunNumber; }
  uint32_t GetMinOrbit() { return fMinOrbit; }
  uint32_t GetMaxOrbit() { return fMaxOrbit; }
  Long64_t Merge(TCollection* list);

 private:
  uint32_t fRunNumber;
  uint32_t fMinOrbit;
  uint32_t fMaxOrbit;
  ClassDef(OrbitRange, 1)
};

#endif // COMMON_CORE_ORBITRANGE_H_
