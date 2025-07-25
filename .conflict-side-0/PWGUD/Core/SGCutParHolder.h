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

#ifndef PWGUD_CORE_SGCUTPARHOLDER_H_
#define PWGUD_CORE_SGCUTPARHOLDER_H_

#include <Rtypes.h>

#include <vector>

// object to hold customizable cut values
class SGCutParHolder
{
 public:
  // constructor
  SGCutParHolder(int ndtcoll = 1, int nMinBCs = 2,
                 bool withFwdTracks = false,
                 bool globalTracksOnly = false,
                 bool ITSonlyTracks = true,
                 int MinNTracks = 0, int MaxNTracks = 10000,
                 int pidHypo = 211,
                 float MinPosz = -1000., float MaxPosz = 1000.,
                 float minPt = 0., float maxPt = 1000.,
                 float minEta = -1.0, float maxEta = 1.0,
                 float maxFITtime = 4,
                 float minRgtrwTOF = 0.,
                 std::vector<float> FITAmpLimits = {0., 0., 0., 0., 0.}) : mNDtcoll{ndtcoll}, mMinNBCs{nMinBCs}, mWithFwdTracks{withFwdTracks}, mGlobalTracksOnly{globalTracksOnly}, mITSOnlyTracks{ITSonlyTracks}, mMinNTracks{MinNTracks}, mMaxNTracks{MaxNTracks}, mPidHypo{pidHypo}, mMinVertexPosz{MinPosz}, mMaxVertexPosz{MaxPosz}, mMinPt{minPt}, mMaxPt{maxPt}, mMinEta{minEta}, mMaxEta{maxEta}, mMaxFITtime{maxFITtime}, mMinRgtrwTOF{minRgtrwTOF}, mFITAmpLimits{FITAmpLimits}
  {
  }

  // setter
  void SetNDtcoll(int);
  void SetMinNBCs(int);
  void SetWithFwdTracks(bool);
  void SetGlobalTracksOnly(bool);
  void SetITSOnlyTracks(bool);
  void SetNTracks(int MinNTracks, int MaxNTracks);
  void SetPidHypothesis(int pidHypo);
  void SetPoszRange(float MinPosz, float MaxPosz);
  void SetPtRange(float minPt, float maxPt);
  void SetEtaRange(float minEta, float maxEta);
  void SetMaxFITtime(float maxFITtime);
  void SetMinRgtrwTOF(float minRgtrwTOF);
  void SetFITAmpLimits(std::vector<float> FITAmpLimits);

  // getter
  int NDtcoll() const;
  int minNBCs() const;
  bool withFwdTracks() const;
  bool globalTracksOnly() const;
  bool ITSOnlyTracks() const;
  int minNTracks() const;
  int maxNTracks() const;
  int pidHypothesis() const;
  float minPosz() const;
  float maxPosz() const;
  float minPt() const;
  float maxPt() const;
  float minEta() const;
  float maxEta() const;
  float maxFITtime() const;
  float minRgtrwTOF() const;
  std::vector<float> FITAmpLimits() const;

 private:
  // number of collision time resolutions to consider
  int mNDtcoll;
  int mMinNBCs;

  // allow forward tracks
  bool mWithFwdTracks;

  // require all vertex tracks to be global tracks
  bool mGlobalTracksOnly;
  bool mITSOnlyTracks;

  // number of tracks
  int mMinNTracks, mMaxNTracks; // Number of allowed tracks

  // PID hypothesis
  int mPidHypo;

  // vertex z-position
  float mMinVertexPosz, mMaxVertexPosz; // Vertex z-position

  // kinematic cuts
  float mMinPt, mMaxPt;   // range of track pT
  float mMinEta, mMaxEta; // range of track eta

  // maximum FIT time
  float mMaxFITtime;

  float mMinRgtrwTOF;

  // lower limits for FIT signals
  std::vector<float> mFITAmpLimits;

  ClassDefNV(SGCutParHolder, 1);
};

#endif // PWGUD_CORE_SGCUTPARHOLDER_H_
