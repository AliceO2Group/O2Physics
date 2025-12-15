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
#include <vector>
#include <TVector3.h>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/Logger.h"

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
  DGAnaparHolder(int MinNTracks = 0, int MaxNTracks = 10000, float minrgtrwTOF = 0.,
                 float maxDCAxy = 100., float maxDCAz = 100,
                 int dBCMin = 0, int dBCMax = 0,
                 std::vector<int> FITvetoes = {0, 1, 1, 0, 0},
                 bool ITSonlyTracks = true,
                 int minNClTPC = 0, int maxNClTPC = 200,
                 float minChi2NClTPC = 0., float maxChi2NClTPC = 100.,
                 float minpt = 0.0, float maxpt = 100.0,
                 float mineta = -2.0, float maxeta = 2.0,
                 float minalpha = 0.0, float maxalpha = 3.2,
                 float minptsys = 0.0, float maxptsys = 100.0,
                 std::size_t nCombine = 2,
                 std::vector<int> netCharges = {0},
                 std::vector<int> unlikeCharges = {0},
                 std::vector<int> likeCharges = {-2, 2},
                 std::vector<int> DGPIDs = {211, 211},
                 std::vector<float> DGPIDCutValues = {}) : mMinNTracks{MinNTracks}, mMaxNTracks{MaxNTracks}, mMinRgtrwTOF{minrgtrwTOF}, mMaxDCAxy{maxDCAxy}, mMaxDCAz{maxDCAz}, mdBCMin{dBCMin}, mdBCMax{dBCMax}, mFITvetoes{FITvetoes}, mITSOnlyTracks{ITSonlyTracks}, mMinNClTPC{minNClTPC}, mMaxNClTPC{maxNClTPC}, mMinChi2NClTPC{minChi2NClTPC}, mMaxChi2NClTPC{maxChi2NClTPC}, mMinpt{minpt}, mMaxpt{maxpt}, mMineta{mineta}, mMaxeta{maxeta}, mMinAlpha{minalpha}, mMaxAlpha{maxalpha}, mMinptsys{minptsys}, mMaxptsys{maxptsys}, mNCombine{nCombine}, mNetCharges{netCharges}, mUnlikeCharges{unlikeCharges}, mLikeCharges{likeCharges}, mDGPIDs{DGPIDs}, mDGPIDCutValues{DGPIDCutValues}
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

  // helper
  void Print();

  // setter
  void SetNTracks(int, int);
  void SetMinRgtrwTOF(float);
  void SetmaxDCA(float, float);
  void SetdBC(int, int);
  void SetFITvetoes(std::vector<int>);
  void SetITSOnlyTracks(bool);
  void SetNClTPC(int, int);
  void SetChi2NClTPC(float, float);
  void Setpt(float, float);
  void Seteta(float, float);
  void SetAlpha(float, float);
  void Setptsys(float, float);
  void SetnCombine(std::size_t);
  void SetnetCharges(std::vector<int>);
  void SetunlikeCharges(std::vector<int>);
  void SetlikeCharges(std::vector<int>);
  void SetPIDs(std::vector<int>);

  // getter
  int minNTracks() const { return mMinNTracks; }
  int maxNTracks() const { return mMaxNTracks; }
  float minRgtrwTOF() const { return mMinRgtrwTOF; }
  float maxDCAxy() const { return mMaxDCAxy; }
  float maxDCAz() const { return mMaxDCAz; }
  int dBCMin() const { return mdBCMin; }
  int dBCMax() const { return mdBCMax; }
  std::vector<int> FITvetoes() { return mFITvetoes; }
  bool ITSOnlyTracks() { return mITSOnlyTracks; }
  int minNClTPC() { return mMinNClTPC; }
  int maxNClTPC() { return mMaxNClTPC; }
  float minChi2NClTPC() { return mMinChi2NClTPC; }
  float maxChi2NClTPC() { return mMaxChi2NClTPC; }
  float minpt() const { return mMinpt; }
  float maxpt() const { return mMaxpt; }
  float mineta() const { return mMineta; }
  float maxeta() const { return mMaxeta; }
  float minAlpha() const { return mMinAlpha; }
  float maxAlpha() const { return mMaxAlpha; }
  float minptsys() const { return mMinptsys; }
  float maxptsys() const { return mMaxptsys; }
  std::size_t nCombine() const { return mNCombine; }
  std::vector<int> netCharges() const { return mNetCharges; }
  std::vector<int> unlikeCharges() const { return mUnlikeCharges; }
  std::vector<int> likeCharges() const { return mLikeCharges; }
  std::vector<int> PIDs() const { return mDGPIDs; }
  DGPIDCuts PIDCuts();
  std::vector<int> uniquePermutations();

 private:
  // helper functions
  void permutations(std::vector<int>& ref, int n0, int np, std::vector<std::vector<int>>& perms);
  int permutations(int n0, std::vector<std::vector<int>>& perms);
  void makeUniquePermutations();

  // arguments
  int mMinNTracks;
  int mMaxNTracks;
  float mMinRgtrwTOF;
  float mMaxDCAxy;
  float mMaxDCAz;
  int mdBCMin;
  int mdBCMax;
  std::vector<int> mFITvetoes;
  bool mITSOnlyTracks;
  int mMinNClTPC;
  int mMaxNClTPC;
  float mMinChi2NClTPC;
  float mMaxChi2NClTPC;
  float mMinpt;
  float mMaxpt;
  float mMineta;
  float mMaxeta;
  float mMinAlpha;
  float mMaxAlpha;
  float mMinptsys;
  float mMaxptsys;
  std::size_t mNCombine;
  std::vector<int> mNetCharges;    // all PV tracks
  std::vector<int> mUnlikeCharges; // selected PV tracks
  std::vector<int> mLikeCharges;   // selected PV tracks
  std::vector<int> mDGPIDs;
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
  DGParticle(TDatabasePDG* pdg, DGAnaparHolder anaPars, TTrack const& tracks, std::vector<int> comb)
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
  std::vector<int> trkinds() { return mtrkinds; }
  float M() { return mIVM.M(); }
  float Perp() { return mIVM.Perp(); }

 private:
  // invariant mass
  TLorentzVector mIVM;

  // indices of tracks included
  std::vector<int> mtrkinds;

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
  bool isGoodCombination(std::vector<int> comb, TTrack const& tracks, std::vector<int> acceptedCharges)
  {
    // compute net charge of track combination
    int netCharge = 0.;
    for (auto const& ind : comb) {
      netCharge += (tracks.begin() + ind).sign();
    }

    // is this in the list of accepted net charges?
    if (std::find(acceptedCharges.begin(), acceptedCharges.end(), netCharge) != acceptedCharges.end()) {
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

    // check ITS only
    if (!mAnaPars.ITSOnlyTracks() && !track.hasTPC()) {
      return false;
    }

    // check ncluster TPC
    auto nTPCCL = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
    if (nTPCCL < mAnaPars.minNClTPC() || nTPCCL > mAnaPars.maxNClTPC()) {
      return false;
    }

    // check chi2 per ncluster TPC
    auto chi2NClTPC = track.tpcChi2NCl();
    if (chi2NClTPC < mAnaPars.minChi2NClTPC() || chi2NClTPC > mAnaPars.maxChi2NClTPC()) {
      return false;
    }

    // check pt of track
    if (track.pt() < mAnaPars.minpt() || track.pt() > mAnaPars.maxpt()) {
      return false;
    }

    // check eta of track
    auto v = TVector3(track.px(), track.py(), track.pz());
    if (v.Eta() < mAnaPars.mineta() || v.Eta() > mAnaPars.maxeta()) {
      return false;
    }

    // cut on dcaXY and dcaZ
    LOGF(debug, "mAnaPars.maxDCAxyz %f %f", mAnaPars.maxDCAxy(), mAnaPars.maxDCAz());
    if (track.dcaXY() < -std::abs(mAnaPars.maxDCAxy()) || track.dcaXY() > std::abs(mAnaPars.maxDCAxy())) {
      return false;
    }
    if (track.dcaZ() < -std::abs(mAnaPars.maxDCAz()) || track.dcaZ() > std::abs(mAnaPars.maxDCAz())) {
      return false;
    }

    // loop over all PIDCuts and apply the ones which apply to this track
    auto pidcuts = mAnaPars.PIDCuts().Cuts();
    for (auto pidcut : pidcuts) {

      // skip cut if it does not apply to this track
      if (pidcut.nPart() != cnt || pidcut.cutApply() <= 0) {
        continue;
      }

      // check pt for pid cut
      if (track.pt() < pidcut.cutPtMin() || track.pt() > pidcut.cutPtMax()) {
        continue;
      }

      // is detector information required
      if (pidcut.cutApply() == 2) {
        if (pidcut.cutDetector() == 1 && !track.hasTPC()) {
          return false;
        }
        if (pidcut.cutDetector() == 2 && !track.hasTOF()) {
          return false;
        }
      }

      // get detector value
      float detValue = 0.;
      if (pidcut.cutDetector() == 1) {
        if (!track.hasTPC()) {
          continue;
        }
        switch (std::abs(pidcut.cutType())) {
          case 1:
            detValue = getTPCnSigma(track, pidcut.cutPID());
            break;
          case 2:
            detValue = track.tpcSignal();
        }
        LOGF(debug, "detValue TPC %f", detValue);
      } else if (std::abs(pidcut.cutDetector()) == 2) {
        if (!track.hasTOF()) {
          continue;
        }
        switch (std::abs(pidcut.cutType())) {
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

    // the track is good if we arrive here
    return true;
  }

  template <typename TTrack>
  std::vector<int> computeIVMs(TTrack const& tracks)
  {
    // reset
    mUnlikeIVMs.clear();
    mLikeIVMs.clear();

    // create combinations including permutations
    auto combs = combinations(tracks.size());

    // loop over unique combinations
    for (auto comb : combs) {
      // are tracks compatible with PID requirements?
      bool isGoodTracks = true;
      auto cnt = -1;
      for (auto ind : comb) {
        cnt++;
        if (!isGoodTrack(tracks.begin() + ind, cnt)) {
          isGoodTracks = false;
          break;
        }
      }

      // is combination compatible with netCharge requirements?
      if (isGoodTracks) {
        DGParticle IVM(fPDG, mAnaPars, tracks, comb);
        // unlike sign
        if (isGoodCombination(comb, tracks, mAnaPars.unlikeCharges())) {
          mUnlikeIVMs.push_back(IVM);
        }
        if (isGoodCombination(comb, tracks, mAnaPars.likeCharges())) {
          mLikeIVMs.push_back(IVM);
        }
      }
    }

    return std::vector<int>{static_cast<int>(mUnlikeIVMs.size()), static_cast<int>(mLikeIVMs.size())};
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
  std::vector<DGParticle> unlikeIVMs() { return mUnlikeIVMs; }
  std::vector<DGParticle> likeIVMs() { return mLikeIVMs; }

  int pid2ind(int pid);

 private:
  // analysis parameters
  DGAnaparHolder mAnaPars;

  // list of DGParticles
  std::vector<DGParticle> mUnlikeIVMs;
  std::vector<DGParticle> mLikeIVMs;

  // particle properties
  TDatabasePDG* fPDG;

  // helper functions for computeIVMs
  void combinations(int n0, std::vector<int>& pool, int np, std::vector<int>& inds, int n,
                    std::vector<std::vector<int>>& combs);
  int combinations(int n0, int np, std::vector<std::vector<int>>& combs);
  std::vector<std::vector<int>> combinations(int nPool);

  // ClassDefNV(DGPIDSelector, 1);
};

// -----------------------------------------------------------------------------
#endif // PWGUD_CORE_DGPIDSELECTOR_H_
