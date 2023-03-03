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

#ifndef PWGUD_CORE_DGPIDSELECTOR_H_
#define PWGUD_CORE_DGPIDSELECTOR_H_

#include <gandiva/projector.h>
#include <string>
#include <vector>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;

const int numDGPIDCutParameters = 9;
float particleMass(TDatabasePDG* pdg, int pid);

// -----------------------------------------------------------------------------
//  numPart:    Particle number to which the parameters apply
//  cutPID:     DPG code of particle hypothesis for nSigma calculation
//  cutDetector:detector: 1: TPC
//                        2: TOF
//  cutType:    cut type:  1: pt and nSigma within limits
//                        -1: nSigma out of limits within pt range
//                         2: pt and detector signal  within limits
//                        -2: detector signal out of limits within pt range
//  cutApply:   How to apply cut: 0: not active
//                                1: if information available
//                                2: return false if information not available
//  ptMin, ptMax, nSigmaTPCmin, nSigmaTPCmax, nSigmaTOFmin, nSigmaTOFmax: cut limits
struct DGPIDCut {

 public:
  // constructor
  DGPIDCut();
  DGPIDCut(float numPart, float cutPID, float cutDetector, float cutType, float cutApply,
           float ptMin, float ptMax, float nSigmamin, float nSigmamax);
  explicit DGPIDCut(float* cutValues);
  ~DGPIDCut();

  // setters

  // getters
  void Print();
  int nPart() { return mnumPart; }
  int cutPID() { return mcutPID; }
  int cutDetector() { return mcutDetector; }
  int cutType() { return mcutType; }
  int cutApply() { return mcutApply; }
  float cutPtMin() { return mptMin; }
  float cutPtMax() { return mptMax; }
  float cutdetValueMin() { return mdetValuemin; }
  float cutdetValueMax() { return mdetValuemax; }

 private:
  int mnumPart;
  int mcutPID;
  int mcutDetector;
  int mcutType;
  int mcutApply;
  float mptMin;
  float mptMax;
  float mdetValuemin;
  float mdetValuemax;

  // ClassDefNV(DGPIDCut, 1);
};

// =============================================================================
// Cuts in mPIDCuts are combined (&&)
struct DGPIDCuts {
 public:
  // constructor
  DGPIDCuts();
  explicit DGPIDCuts(std::vector<float> PIDCutValues);
  ~DGPIDCuts();

  // setter
  void clear()
  {
    mDGPIDCuts.clear();
  }
  void setPIDCuts(std::vector<float> PIDCutValues);

  // getter
  void Print();
  std::vector<DGPIDCut> Cuts() { return mDGPIDCuts; }

 private:
  std::vector<DGPIDCut> mDGPIDCuts;

  // ClassDefNV(DGPIDCuts, 1);
};

// =============================================================================
// object to hold customizable analysis parameters
struct DGAnaparHolder {
 public:
  // constructor
  DGAnaparHolder(int nCombine = 2, float maxDCAxy = 100., float maxDCAz = 100,
                 int dBCMin = 0, int dBCMax = 0,
                 float minptsys = 0.0, float maxptsys = 100.0,
                 float minpt = 0.0, float maxpt = 100.0,
                 float minalpha = 0.0, float maxalpha = 3.2,
                 std::vector<int> netCharges = {-2, -1, 0, 1, 2},
                 std::vector<float> DGPIDs = {211, 211},
                 std::vector<float> DGPIDCutValues = {}) : mNCombine{nCombine}, mdBCMin{dBCMin}, mdBCMax{dBCMax}, mMaxDCAxy{maxDCAxy}, mMaxDCAz{maxDCAz}, mMinpt{minpt}, mMaxpt{maxpt}, mMinptsys{minptsys}, mMaxptsys{maxptsys}, mMinAlpha{minalpha}, mMaxAlpha{maxalpha}, mNetCharges{netCharges}, mDGPIDs{DGPIDs}, mDGPIDCutValues{DGPIDCutValues}
  {
    if (mdBCMin < -16) {
      mdBCMin = -16;
    } else if (mdBCMin > 15) {
      mdBCMin = 15;
    }
    if (mdBCMax < -16) {
      mdBCMax = -16;
    } else if (mdBCMax > 15) {
      mdBCMax = 15;
    }

    makeUniquePermutations();
  }
  ~DGAnaparHolder();

  // getter
  void Print();
  int nCombine() const { return mNCombine; }
  int dBCMin() const { return mdBCMin; }
  int dBCMax() const { return mdBCMax; }
  float maxDCAxy() { return mMaxDCAxy; }
  float maxDCAz() { return mMaxDCAz; }
  float minptsys() { return mMinptsys; }
  float maxptsys() { return mMaxptsys; }
  float minpt() { return mMinpt; }
  float maxpt() { return mMaxpt; }
  float minAlpha() { return mMinAlpha; }
  float maxAlpha() { return mMaxAlpha; }
  std::vector<int> netCharges() { return mNetCharges; }
  std::vector<float> PIDs() { return mDGPIDs; }
  DGPIDCuts PIDCuts();
  std::vector<int> uniquePermutations();

 private:
  // helper functions
  void permutations(std::vector<uint>& ref, int n0, int np, std::vector<std::vector<uint>>& perms);
  int permutations(int n0, std::vector<std::vector<uint>>& perms);
  void makeUniquePermutations();

  // number of tracks to combine
  int mNCombine;

  // BBFlags
  int mdBCMin;
  int mdBCMax;

  // dca of tracks
  float mMaxDCAxy;
  float mMaxDCAz;

  // pt-range of tracks
  float mMinpt;
  float mMaxpt;

  // pt-range of system
  float mMinptsys;
  float mMaxptsys;

  // alpha-range when 2 tracks
  float mMinAlpha;
  float mMaxAlpha;

  // net charge of all tracks
  std::vector<int> mNetCharges;

  // PID information
  std::vector<float> mDGPIDs;
  std::vector<float> mDGPIDCutValues;

  // unique permutations
  std::vector<int> muniquePerms;

  // ClassDefNV(DGAnaparHolder, 1);
};

// =============================================================================
// a structure which holds the indices of tracks and their invariant mass
struct DGParticle {
 public:
  DGParticle();
  template <typename TTrack>
  DGParticle(TDatabasePDG* pdg, DGAnaparHolder anaPars, TTrack const& tracks, std::vector<uint> comb)
  {
    // compute invariant mass
    TLorentzVector lvtmp;
    auto pids = anaPars.PIDs();

    // loop over tracks and update mIVM
    mIVM = TLorentzVector(0., 0., 0., 0.);
    auto cnt = -1;
    for (auto ind : comb) {
      cnt++;
      auto track = tracks.begin() + ind;
      lvtmp.SetXYZM(track.px(), track.py(), track.pz(), particleMass(pdg, pids[cnt]));
      mIVM += lvtmp;
    }

    // set array of track indices
    mtrkinds = comb;
  }
  ~DGParticle();

  // getter
  void Print();
  std::vector<uint> trkinds() { return mtrkinds; }
  float M() { return mIVM.M(); }
  float Perp() { return mIVM.Perp(); }

 private:
  // invariant mass
  TLorentzVector mIVM;

  // indices of tracks included
  std::vector<uint> mtrkinds;

  // ClassDefNV(DGParticle, 1);
};

// =============================================================================
// A class to check PIDs of tracks and provide track combinations
struct DGPIDSelector {
 public:
  DGPIDSelector();
  ~DGPIDSelector();

  // setters
  void init(DGAnaparHolder anaPars);

  // getters
  void Print();
  template <typename TTrack>
  bool isGoodCombination(std::vector<uint> comb, TTrack const& tracks)
  {
    // compute net charge of track combination
    int netCharge = 0.;
    for (auto const& ind : comb) {
      netCharge += (tracks.begin() + ind).sign();
    }
    LOGF(debug, "Net charge %i", netCharge);

    // is this in the list of accepted net charges?
    auto netCharges = mAnaPars.netCharges();
    if (std::find(netCharges.begin(), netCharges.end(), netCharge) != netCharges.end()) {
      return true;
    }
    return false;
  }

  template <typename TTrack>
  bool isGoodTrack(TTrack track, int cnt)
  {
    // get pid of particle cnt
    auto pid = mAnaPars.PIDs()[cnt];

    // unknown PID
    auto pidhypo = pid2ind(pid);
    if (pidhypo < 0) {
      return false;
    }

    // cut on dcaXY and dcaZ
    // LOGF(debug, "mAnaPars.maxDCAxyz %f %f", mAnaPars.maxDCAxy(), mAnaPars.maxDCAz());
    // if (track.dcaXY() < -abs(mAnaPars.maxDCAxy()) || track.dcaXY() > abs(mAnaPars.maxDCAxy())) {
    //  return false;
    //}

    // if (track.dcaZ() < -abs(mAnaPars.maxDCAz()) || track.dcaZ() > abs(mAnaPars.maxDCAz())) {
    //   return false;
    // }

    // loop over all PIDCuts and apply the ones which apply to this track
    auto pidcuts = mAnaPars.PIDCuts().Cuts();
    for (auto pidcut : pidcuts) {

      // skip cut if it does not apply to this track
      LOGF(debug, "nPart %i %i, Type %i Apply %i", pidcut.nPart(), cnt, pidcut.cutType(), pidcut.cutApply());
      if (pidcut.nPart() != cnt || pidcut.cutApply() <= 0) {
        continue;
      }

      // check pt of track
      if (track.pt() < mAnaPars.minpt() || track.pt() > mAnaPars.maxpt()) {
        return false;
      }

      // check pt for pid cut
      LOGF(debug, "pT %f %f %f", track.pt(), pidcut.cutPtMin(), pidcut.cutPtMax());
      if (track.pt() < pidcut.cutPtMin() || track.pt() > pidcut.cutPtMax()) {
        continue;
      }

      // is detector information required
      LOGF(debug, "TPC %i TOF %i", track.hasTPC(), track.hasTOF());
      if (pidcut.cutApply() == 2) {
        if (pidcut.cutDetector() == 1 && !track.hasTPC()) {
          return false;
        }
        if (pidcut.cutDetector() == 2 && !track.hasTOF()) {
          return false;
        }
      }

      // get detector value
      LOGF(debug, "cutPID %i", pidcut.cutPID());
      float detValue = 0.;
      if (pidcut.cutDetector() == 1) {
        if (!track.hasTPC()) {
          continue;
        }
        switch (abs(pidcut.cutType())) {
          case 1:
            detValue = getTPCnSigma(track, pidcut.cutPID());
            break;
          case 2:
            detValue = track.tpcSignal();
        }
        LOGF(debug, "detValue TPC %f", detValue);
      } else if (abs(pidcut.cutDetector()) == 2) {
        if (!track.hasTOF()) {
          continue;
        }
        switch (abs(pidcut.cutType())) {
          case 1:
            detValue = getTOFnSigma(track, pidcut.cutPID());
            break;
          case 2:
            detValue = track.tofSignal();
        }
        LOGF(debug, "detValue TOF %f", detValue);
      } else {
        continue;
      }
      LOGF(debug, "detValue %f", detValue);

      // inclusive / exclusive
      if (pidcut.cutType() > 0 && (detValue < pidcut.cutdetValueMin() || detValue > pidcut.cutdetValueMax())) {
        return false;
      } else if (pidcut.cutType() < 0 && (detValue > pidcut.cutdetValueMin() && detValue < pidcut.cutdetValueMax())) {
        return false;
      }
    }

    return true;
  }
  template <typename TTrack>
  int computeIVMs(TTrack const& tracks)
  {
    // reset
    mIVMs.clear();

    // create combinations including permutations
    auto combs = combinations(tracks.size());

    // loop over unique combinations
    for (auto comb : combs) {
      // is combination compatible with netCharge requirements?
      if (!isGoodCombination(comb, tracks)) {
        continue;
      }
      // is tracks compatible with PID requirements?
      bool isGoodComb = true;
      auto cnt = -1;
      for (auto ind : comb) {
        cnt++;
        if (!isGoodTrack(tracks.begin() + ind, cnt)) {
          isGoodComb = false;
          break;
        }
      }

      // update list of IVMs
      if (isGoodComb) {
        DGParticle IVM(fPDG, mAnaPars, tracks, comb);
        mIVMs.push_back(IVM);
      }
    }

    return mIVMs.size();
  }

  DGAnaparHolder getAnaPars() { return mAnaPars; }
  template <typename TTrack>
  float getTPCnSigma(TTrack track, int pid)
  {
    auto hypo = pid2ind(pid);
    switch (hypo) {
      case 0:
        return track.tpcNSigmaEl();
      case 1:
        return track.tpcNSigmaPi();
      case 2:
        return track.tpcNSigmaMu();
      case 3:
        return track.tpcNSigmaKa();
      case 4:
        return track.tpcNSigmaPr();
      default:
        return 0.;
    }
  }
  template <typename TTrack>
  float getTOFnSigma(TTrack track, int pid)
  {
    auto hypo = pid2ind(pid);
    switch (hypo) {
      case 0:
        return track.tofNSigmaEl();
      case 1:
        return track.tofNSigmaPi();
      case 2:
        return track.tofNSigmaMu();
      case 3:
        return track.tofNSigmaKa();
      case 4:
        return track.tofNSigmaPr();
      default:
        return 0.;
    }
  }
  std::vector<DGParticle> IVMs() { return mIVMs; }

  int pid2ind(int pid);

 private:
  // analysis parameters
  DGAnaparHolder mAnaPars;

  // list of DGParticles
  std::vector<DGParticle> mIVMs;

  // particle properties
  TDatabasePDG* fPDG;

  // helper functions for computeIVMs
  void combinations(int n0, std::vector<uint>& pool, int np, std::vector<uint>& inds, int n,
                    std::vector<std::vector<uint>>& combs);
  int combinations(int n0, int np, std::vector<std::vector<uint>>& combs);
  std::vector<std::vector<uint>> combinations(int nPool);

  // ClassDefNV(DGPIDSelector, 1);
};

// -----------------------------------------------------------------------------
#endif // PWGUD_CORE_DGPIDSELECTOR_H_
