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

#include "cutHolder.h"

// setter
void cutHolder::SetNDtcoll(int ndtcoll)
{
  mNDtcoll = ndtcoll;
}
void cutHolder::SetMinNBCs(int nminbcs)
{
  mMinNBCs = nminbcs;
}
void cutHolder::SetNTracks(int MinNTracks, int MaxNTracks)
{
  mMinNTracks = MinNTracks;
  mMaxNTracks = MaxNTracks;
}
void cutHolder::SetNetCharge(int MinNetCharge, int MaxNetCharge)
{
  mMinNetCharge = MinNetCharge;
  mMaxNetCharge = MaxNetCharge;
}
void cutHolder::SetPidHypothesis(int pidHypo)
{
  mPidHypo = pidHypo;
}
void cutHolder::SetPoszRange(float MinPosz, float MaxPosz)
{
  mMinVertexPosz = MinPosz;
  mMaxVertexPosz = MaxPosz;
}

void cutHolder::SetPtRange(float minPt, float maxPt)
{
  mMinPt = minPt;
  mMaxPt = maxPt;
}
void cutHolder::SetEtaRange(float minEta, float maxEta)
{
  mMinEta = minEta;
  mMaxEta = maxEta;
}
void cutHolder::SetIVMRange(float minIVM, float maxIVM)
{
  mMinIVM = minIVM;
  mMaxIVM = maxIVM;
}
void cutHolder::SetMaxNSigmaTPC(float maxnSigma)
{
  mMaxNSigmaTPC = maxnSigma;
}
void cutHolder::SetMaxNSigmaTOF(float maxnSigma)
{
  mMaxNSigmaTOF = maxnSigma;
}

void cutHolder::SetFITAmpLimits(std::vector<float> FITAmpLimits)
{
  mFITAmpLimits = FITAmpLimits;
}

// getter
int cutHolder::NDtcoll() const { return mNDtcoll; }
int cutHolder::minNBCs() const { return mMinNBCs; }
int cutHolder::minNTracks() const { return mMinNTracks; }
int cutHolder::maxNTracks() const { return mMaxNTracks; }
int cutHolder::minNetCharge() const { return mMinNetCharge; }
int cutHolder::maxNetCharge() const { return mMaxNetCharge; }
int cutHolder::pidHypothesis() const { return mPidHypo; }
float cutHolder::minPosz() const { return mMinVertexPosz; }
float cutHolder::maxPosz() const { return mMaxVertexPosz; }
float cutHolder::minPt() const { return mMinPt; }
float cutHolder::maxPt() const { return mMaxPt; }
float cutHolder::minEta() const { return mMinEta; }
float cutHolder::maxEta() const { return mMaxEta; }
float cutHolder::minIVM() const { return mMinIVM; }
float cutHolder::maxIVM() const { return mMaxIVM; }
float cutHolder::maxNSigmaTPC() const { return mMaxNSigmaTPC; };
float cutHolder::maxNSigmaTOF() const { return mMaxNSigmaTOF; };
std::vector<float> cutHolder::FITAmpLimits() const { return mFITAmpLimits; };
