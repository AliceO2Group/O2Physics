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

#ifndef PWGUD_CORE_DGCUTPARHOLDER_H_
#define PWGUD_CORE_DGCUTPARHOLDER_H_

#include <Rtypes.h>
#include <vector>

// object to hold customizable cut values
class DGCutparHolder
{
 public:
  // constructor
  DGCutparHolder(int ndtcoll = 4, int nMinBCs = 7,
                 bool withFwdTracks = false,
                 bool globalTracksOnly = false,
                 bool ITSonlyTracks = true,
                 float minrgtrwTOF = 0.,
                 int MinNTracks = 0, int MaxNTracks = 10000,
                 std::vector<int> NetCharges = {0},
                 int pidHypo = 211,
                 float MinPosz = -1000., float MaxPosz = 1000.,
                 float minPt = 0., float maxPt = 1000.,
                 float minEta = -1.0, float maxEta = 1.0,
                 float minIVM = 0.0, float maxIVM = 1000.,
                 float maxNSigmaTPC = 1000., float maxNSigmaTOF = 1000.,
                 bool TVX = false,
                 bool TSC = false,
                 bool TCE = false,
                 bool TOR = true,
                 float maxFITtime = 4,
                 std::vector<float> FITAmpLimits = {0., 0., 0., 0., 0.},
                 std::vector<int> collisionSel = {0, 0, 0, 0, 0}) : mNDtcoll{ndtcoll}, mMinNBCs{nMinBCs}, mWithFwdTracks{withFwdTracks}, mGlobalTracksOnly{globalTracksOnly}, mITSOnlyTracks{ITSonlyTracks}, mMinRgtrwTOF{minrgtrwTOF}, mMinNTracks{MinNTracks}, mMaxNTracks{MaxNTracks}, mNetCharges{NetCharges}, mPidHypo{pidHypo}, mMinVertexPosz{MinPosz}, mMaxVertexPosz{MaxPosz}, mMinPt{minPt}, mMaxPt{maxPt}, mMinEta{minEta}, mMaxEta{maxEta}, mMinIVM{minIVM}, mMaxIVM{maxIVM}, mMaxNSigmaTPC{maxNSigmaTPC}, mMaxNSigmaTOF{maxNSigmaTOF}, mTVX{TVX}, mTSC{TSC}, mTCE{TCE}, mTOR{TOR}, mMaxFITtime{maxFITtime}, mFITAmpLimits{FITAmpLimits}, mCollisionSel{collisionSel}
  {
  }

  // setter
  void SetNDtcoll(int);
  void SetMinNBCs(int);
  void SetWithFwdTracks(bool);
  void SetGlobalTracksOnly(bool);
  void SetITSOnlyTracks(bool);
  void SetMinRgtrwTOF(float);
  void SetNTracks(int MinNTracks, int MaxNTracks);
  void SetNetCharges(std::vector<int> netCharges);
  void SetPidHypothesis(int pidHypo);
  void SetPoszRange(float MinPosz, float MaxPosz);
  void SetPtRange(float minPt, float maxPt);
  void SetEtaRange(float minEta, float maxEta);
  void SetIVMRange(float minIVM, float maxIVM);
  void SetMaxNSigmaTPC(float maxnSigma);
  void SetMaxNSigmaTOF(float maxnSigma);
  void SetTVX(bool tvx);
  void SetTSC(bool tsc);
  void SetTCE(bool tce);
  void SetTOR(bool tor);
  void SetMaxFITtime(float maxFITtime);
  void SetFITAmpLimits(std::vector<float> FITAmpLimits);
  void SetCollisionSel(std::vector<int> collisionSel);

  // getter
  int NDtcoll() const;
  int minNBCs() const;
  bool withFwdTracks() const;
  bool globalTracksOnly() const;
  bool ITSOnlyTracks() const;
  float minRgtrwTOF() const;
  int minNTracks() const;
  int maxNTracks() const;
  std::vector<int> netCharges() const;
  int pidHypothesis() const;
  float minPosz() const;
  float maxPosz() const;
  float minPt() const;
  float maxPt() const;
  float minEta() const;
  float maxEta() const;
  float minIVM() const;
  float maxIVM() const;
  float maxNSigmaTPC() const;
  float maxNSigmaTOF() const;
  bool withTVX() const;
  bool withTSC() const;
  bool withTCE() const;
  bool withTOR() const;
  float maxFITtime() const;
  std::vector<float> FITAmpLimits() const;
  std::vector<int> collisionSel() const;

 private:
  // number of collision time resolutions to consider
  int mNDtcoll;
  int mMinNBCs;

  // allow forward tracks
  bool mWithFwdTracks;

  // require all vertex tracks to be global tracks
  bool mGlobalTracksOnly;
  bool mITSOnlyTracks;

  // required minimum fraction of global tracks with TOF hit
  float mMinRgtrwTOF;

  // number of tracks
  int mMinNTracks, mMaxNTracks; // Number of allowed tracks

  // net charge of all tracks
  std::vector<int> mNetCharges;

  // PID hypothesis
  int mPidHypo;

  // vertex z-position
  float mMinVertexPosz, mMaxVertexPosz; // Vertex z-position

  // kinematic cuts
  float mMinPt, mMaxPt;   // range of track pT
  float mMinEta, mMaxEta; // range of track eta
  float mMinIVM, mMaxIVM; // range of invariant mass

  // maximum nSigma for PID
  float mMaxNSigmaTPC; // maximum nSigma TPC
  float mMaxNSigmaTOF; // maximum nSigma TOF

  // FIT vetoes
  bool mTVX;
  bool mTSC;
  bool mTCE;
  bool mTOR;

  // maximum FIT time
  float mMaxFITtime;

  // lower limits for FIT signals
  std::vector<float> mFITAmpLimits;

  // collision selections to consider from event selection task
  std::vector<int> mCollisionSel;

  ClassDefNV(DGCutparHolder, 1);
};

#endif // PWGUD_CORE_DGCUTPARHOLDER_H_
