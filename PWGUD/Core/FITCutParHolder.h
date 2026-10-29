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
// \FIT bit thresholds
// \author Sandor Lokos, sandor.lokos@cern.ch
// \since  March 2026

#ifndef PWGUD_CORE_FITCUTPARHOLDER_H_
#define PWGUD_CORE_FITCUTPARHOLDER_H_

#include <Rtypes.h>

// object to hold customizable FIT bit thresholds
class FITCutParHolder
{
 public:
  // constructor
  FITCutParHolder(bool saveFITbitsets = true,
                  float thr1_FV0A = 8.,
                  float thr1_FT0A = 8.,
                  float thr1_FT0C = 8.,
                  float thr2_FV0A = 20.,
                  float thr2_FT0A = 20.,
                  float thr2_FT0C = 20.)
    : mSaveFITbitsets{saveFITbitsets},
      mThr1FV0A{thr1_FV0A},
      mThr1FT0A{thr1_FT0A},
      mThr1FT0C{thr1_FT0C},
      mThr2FV0A{thr2_FV0A},
      mThr2FT0A{thr2_FT0A},
      mThr2FT0C{thr2_FT0C}
  {
  }

  // setters
  void SetSaveFITbitsets(bool);
  void SetThr1FV0A(float);
  void SetThr1FT0A(float);
  void SetThr1FT0C(float);
  void SetThr2FV0A(float);
  void SetThr2FT0A(float);
  void SetThr2FT0C(float);

  // getters
  bool saveFITbitsets() const;
  float thr1_FV0A() const;
  float thr1_FT0A() const;
  float thr1_FT0C() const;
  float thr2_FV0A() const;
  float thr2_FT0A() const;
  float thr2_FT0C() const;

 private:
  bool mSaveFITbitsets;

  float mThr1FV0A;
  float mThr1FT0A;
  float mThr1FT0C;

  float mThr2FV0A;
  float mThr2FT0A;
  float mThr2FT0C;

  ClassDefNV(FITCutParHolder, 1);
};

#endif // PWGUD_CORE_FITCUTPARHOLDER_H_
