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

#ifndef O2_ANALYSIS_DIFFCUT_HOLDER_H_
#define O2_ANALYSIS_DIFFCUT_HOLDER_H_

#include <iosfwd>
#include <Rtypes.h>
#include <TMath.h>

// object to hold customizable cut values
class cutHolder
{
 public:
  // constructor
  cutHolder(int cc = 1,
            int ndtcoll = 4,
            int MinNTracks = 0, int MaxNTracks = 10000,
            int MinNetCharge = 0, int MaxNetCharge = 0,
            int pidHypo = 211,
            float MinPosz = -1000., float MaxPosz = 1000.,
            float minPt = 0., float maxPt = 1000.,
            float minEta = -1.0, float maxEta = 1.0,
            float minIVM = 0.0, float maxIVM = 1000.,
            float maxnSigmaTPC = 1000., float maxnSigmaTOF = 1000.) : mcase{cc}, mNDtcoll{ndtcoll}, mMinNTracks{MinNTracks}, mMaxNTracks{MaxNTracks}, mMinNetCharge{MinNetCharge}, mMaxNetCharge{MaxNetCharge}, mPidHypo{pidHypo}, mMinVertexPosz{MinPosz}, mMaxVertexPosz{MaxPosz}, mMinPt{minPt}, mMaxPt{maxPt}, mMinEta{minEta}, mMaxEta{maxEta}, mMinIVM{minIVM}, mMaxIVM{maxIVM}, mMaxnSigmaTPC{maxnSigmaTPC}, mMaxnSigmaTOF{maxnSigmaTOF}
  {
  }

  // setter
  void Setcc(int);
  void SetNDtcoll(int);
  void SetNTracks(int MinNTracks, int MaxNTracks);
  void SetNetCharge(int minNetCharge, int maxNetCharge);
  void SetPidHypothesis(int pidHypo);
  void SetPoszRange(float MinPosz, float MaxPosz);
  void SetPtRange(float minPt, float maxPt);
  void SetEtaRange(float minEta, float maxEta);
  void SetIVMRange(float minIVM, float maxIVM);
  void SetMaxnSigmaTPC(float maxnSigma);
  void SetMaxnSigmaTOF(float maxnSigma);

  // getter
  int cc() const;
  int NDtcoll() const;
  int minNTracks() const;
  int maxNTracks() const;
  int minNetCharge() const;
  int maxNetCharge() const;
  int pidHypothesis() const;
  float minPosz() const;
  float maxPosz() const;
  float minPt() const;
  float maxPt() const;
  float minEta() const;
  float maxEta() const;
  float minIVM() const;
  float maxIVM() const;
  float maxnSigmaTPC() const;
  float maxnSigmaTOF() const;

 private:
  // decay channel
  //  1: 2 pion
  //  2: 4 pion
  //  3: 2 kaon
  //  4: 4 kaon
  int mcase;

  // number of collision time resolutions to consider
  int mNDtcoll;

  // number of tracks
  int mMinNTracks, mMaxNTracks; // Number of allowed tracks

  // net charge of all tracks
  int mMinNetCharge;
  int mMaxNetCharge;

  // PID hypothesis
  int mPidHypo;

  // vertex z-position
  float mMinVertexPosz, mMaxVertexPosz; // Vertex z-position

  // kinematic cuts
  float mMinPt, mMaxPt;   // range of track pT
  float mMinEta, mMaxEta; // range of track eta
  float mMinIVM, mMaxIVM; // range of invariant mass

  // maximum nSigma for PID
  float mMaxnSigmaTPC; // maximum nSigma TPC
  float mMaxnSigmaTOF; // maximum nSigma TOF

  ClassDefNV(cutHolder, 1);
};

#endif // O2_ANALYSIS_DIFFCUT_HOLDER_H_
