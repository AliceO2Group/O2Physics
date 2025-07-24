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

/// \file ProfileSubset.h/.cxx
/// \brief // Helper class to select a subrange of a TProfile
/// \author Emil Gorm Nielsen (ack. V. Vislavicius), NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_PROFILESUBSET_H_
#define PWGCF_GENERICFRAMEWORK_CORE_PROFILESUBSET_H_

#include "TProfile.h"
#include "TProfile2D.h"
#include "TError.h"

class ProfileSubset : public TProfile2D
{
 public:
  explicit ProfileSubset(TProfile2D& inpf) : TProfile2D(inpf) {}
  ~ProfileSubset() {}
  TProfile* GetSubset(bool onx, const char* name, int fb, int lb, int l_nbins = 0, double* l_binarray = 0);
  void OverrideBinContent(double x, double y, double x2, double y2, double val);
  void OverrideBinContent(double x, double y, double x2, double y2, TProfile2D* sourceProf);
  bool OverrideBinsWithZero(int xb1, int yb1, int xb2, int yb2);

  ClassDef(ProfileSubset, 2);
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_PROFILESUBSET_H_
