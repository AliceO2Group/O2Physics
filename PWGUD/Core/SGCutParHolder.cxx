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

#include "SGCutParHolder.h"

#include <vector>

// setter
void SGCutParHolder::SetNDtcoll(int ndtcoll)
{
  mNDtcoll = ndtcoll;
}
void SGCutParHolder::SetMinNBCs(int nminbcs)
{
  mMinNBCs = nminbcs;
}
void SGCutParHolder::SetWithFwdTracks(bool withFwdTracks)
{
  mWithFwdTracks = withFwdTracks;
}
void SGCutParHolder::SetGlobalTracksOnly(bool globalTracksOnly)
{
  mGlobalTracksOnly = globalTracksOnly;
}
void SGCutParHolder::SetITSOnlyTracks(bool ITSonlyTracks)
{
  mITSOnlyTracks = ITSonlyTracks;
}
void SGCutParHolder::SetMinRgtrwTOF(float rgtrwTOF)
{
  mMinRgtrwTOF = rgtrwTOF;
}
void SGCutParHolder::SetNTracks(int MinNTracks, int MaxNTracks)
{
  mMinNTracks = MinNTracks;
  mMaxNTracks = MaxNTracks;
}
void SGCutParHolder::SetPidHypothesis(int pidHypo)
{
  mPidHypo = pidHypo;
}
void SGCutParHolder::SetPoszRange(float MinPosz, float MaxPosz)
{
  mMinVertexPosz = MinPosz;
  mMaxVertexPosz = MaxPosz;
}

void SGCutParHolder::SetPtRange(float minPt, float maxPt)
{
  mMinPt = minPt;
  mMaxPt = maxPt;
}
void SGCutParHolder::SetEtaRange(float minEta, float maxEta)
{
  mMinEta = minEta;
  mMaxEta = maxEta;
}
void SGCutParHolder::SetMaxFITtime(float maxFITtime)
{
  mMaxFITtime = maxFITtime;
}
void SGCutParHolder::SetFITAmpLimits(std::vector<float> FITAmpLimits)
{
  mFITAmpLimits = FITAmpLimits;
}

// getter
int SGCutParHolder::NDtcoll() const { return mNDtcoll; }
int SGCutParHolder::minNBCs() const { return mMinNBCs; }
bool SGCutParHolder::withFwdTracks() const { return mWithFwdTracks; }
bool SGCutParHolder::globalTracksOnly() const { return mGlobalTracksOnly; }
bool SGCutParHolder::ITSOnlyTracks() const { return mITSOnlyTracks; }
float SGCutParHolder::minRgtrwTOF() const { return mMinRgtrwTOF; }
int SGCutParHolder::minNTracks() const { return mMinNTracks; }
int SGCutParHolder::maxNTracks() const { return mMaxNTracks; }
int SGCutParHolder::pidHypothesis() const { return mPidHypo; }
float SGCutParHolder::minPosz() const { return mMinVertexPosz; }
float SGCutParHolder::maxPosz() const { return mMaxVertexPosz; }
float SGCutParHolder::minPt() const { return mMinPt; }
float SGCutParHolder::maxPt() const { return mMaxPt; }
float SGCutParHolder::minEta() const { return mMinEta; }
float SGCutParHolder::maxEta() const { return mMaxEta; }
float SGCutParHolder::maxFITtime() const { return mMaxFITtime; }
std::vector<float> SGCutParHolder::FITAmpLimits() const { return mFITAmpLimits; }
