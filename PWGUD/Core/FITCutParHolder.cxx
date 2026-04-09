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

#include "FITCutParHolder.h"

// setter
void FITCutParHolder::SetSaveFITbitsets(bool saveFITbitsets)
{
  mSaveFITbitsets = saveFITbitsets;
}
void FITCutParHolder::SetThr1FV0A(float thr1_FV0A)
{
  mThr1FV0A = thr1_FV0A;
}
void FITCutParHolder::SetThr1FT0A(float thr1_FT0A)
{
  mThr1FT0A = thr1_FT0A;
}
void FITCutParHolder::SetThr1FT0C(float thr1_FT0C)
{
  mThr1FT0C = thr1_FT0C;
}
void FITCutParHolder::SetThr2FV0A(float thr2_FV0A)
{
  mThr2FV0A = thr2_FV0A;
}
void FITCutParHolder::SetThr2FT0A(float thr2_FT0A)
{
  mThr2FT0A = thr2_FT0A;
}
void FITCutParHolder::SetThr2FT0C(float thr2_FT0C)
{
  mThr2FT0C = thr2_FT0C;
}

// getter
bool FITCutParHolder::saveFITbitsets() const { return mSaveFITbitsets; }
float FITCutParHolder::thr1_FV0A() const { return mThr1FV0A; }
float FITCutParHolder::thr1_FT0A() const { return mThr1FT0A; }
float FITCutParHolder::thr1_FT0C() const { return mThr1FT0C; }
float FITCutParHolder::thr2_FV0A() const { return mThr2FV0A; }
float FITCutParHolder::thr2_FT0A() const { return mThr2FT0A; }
float FITCutParHolder::thr2_FT0C() const { return mThr2FT0C; }
