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

#include "DGCutparHolder.h"

// setter
void DGCutparHolder::SetNDtcoll(int ndtcoll)
{
  mNDtcoll = ndtcoll;
}
void DGCutparHolder::SetMinNBCs(int nminbcs)
{
  mMinNBCs = nminbcs;
}
void DGCutparHolder::SetWithFwdTracks(bool withFwdTracks)
{
  mWithFwdTracks = withFwdTracks;
}
void DGCutparHolder::SetGlobalTracksOnly(bool globalTracksOnly)
{
  mGlobalTracksOnly = globalTracksOnly;
}
void DGCutparHolder::SetMinRgtrwTOF(float rgtrwTOF)
{
  mMinRgtrwTOF = rgtrwTOF;
}
void DGCutparHolder::SetNTracks(int MinNTracks, int MaxNTracks)
{
  mMinNTracks = MinNTracks;
  mMaxNTracks = MaxNTracks;
}
void DGCutparHolder::SetNetCharges(std::vector<int> netCharges)
{
  mNetCharges = netCharges;
}
void DGCutparHolder::SetPidHypothesis(int pidHypo)
{
  mPidHypo = pidHypo;
}
void DGCutparHolder::SetPoszRange(float MinPosz, float MaxPosz)
{
  mMinVertexPosz = MinPosz;
  mMaxVertexPosz = MaxPosz;
}

void DGCutparHolder::SetPtRange(float minPt, float maxPt)
{
  mMinPt = minPt;
  mMaxPt = maxPt;
}
void DGCutparHolder::SetEtaRange(float minEta, float maxEta)
{
  mMinEta = minEta;
  mMaxEta = maxEta;
}
void DGCutparHolder::SetIVMRange(float minIVM, float maxIVM)
{
  mMinIVM = minIVM;
  mMaxIVM = maxIVM;
}
void DGCutparHolder::SetMaxNSigmaTPC(float maxnSigma)
{
  mMaxNSigmaTPC = maxnSigma;
}
void DGCutparHolder::SetMaxNSigmaTOF(float maxnSigma)
{
  mMaxNSigmaTOF = maxnSigma;
}

void DGCutparHolder::SetFITAmpLimits(std::vector<float> FITAmpLimits)
{
  mFITAmpLimits = FITAmpLimits;
}

// getter
int DGCutparHolder::NDtcoll() const { return mNDtcoll; }
int DGCutparHolder::minNBCs() const { return mMinNBCs; }
bool DGCutparHolder::withFwdTracks() const { return mWithFwdTracks; }
bool DGCutparHolder::globalTracksOnly() const { return mGlobalTracksOnly; }
float DGCutparHolder::minRgtrwTOF() const { return mMinRgtrwTOF; }
int DGCutparHolder::minNTracks() const { return mMinNTracks; }
int DGCutparHolder::maxNTracks() const { return mMaxNTracks; }
std::vector<int> DGCutparHolder::netCharges() const { return mNetCharges; }
int DGCutparHolder::pidHypothesis() const { return mPidHypo; }
float DGCutparHolder::minPosz() const { return mMinVertexPosz; }
float DGCutparHolder::maxPosz() const { return mMaxVertexPosz; }
float DGCutparHolder::minPt() const { return mMinPt; }
float DGCutparHolder::maxPt() const { return mMaxPt; }
float DGCutparHolder::minEta() const { return mMinEta; }
float DGCutparHolder::maxEta() const { return mMaxEta; }
float DGCutparHolder::minIVM() const { return mMinIVM; }
float DGCutparHolder::maxIVM() const { return mMaxIVM; }
float DGCutparHolder::maxNSigmaTPC() const { return mMaxNSigmaTPC; }
float DGCutparHolder::maxNSigmaTOF() const { return mMaxNSigmaTOF; }
std::vector<float> DGCutparHolder::FITAmpLimits() const { return mFITAmpLimits; }
