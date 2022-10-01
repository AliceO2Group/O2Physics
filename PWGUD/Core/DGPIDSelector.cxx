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

#include "CommonConstants/PhysicsConstants.h"
#include "DGPIDSelector.h"

// -----------------------------------------------------------------------------
float particleMass(TDatabasePDG* pdg, int pid)
{
  auto mass = 0.;
  TParticlePDG* pdgparticle = pdg->GetParticle(pid);
  if (pdgparticle != nullptr) {
    mass = pdgparticle->Mass();
  }
  return mass;
};

// =============================================================================
// DGPIDCut
DGPIDCut::DGPIDCut()
{
}

DGPIDCut::DGPIDCut(float numPart, float cutPID, float cutDetector, float cutType, float cutApply,
                   float ptMin, float ptMax, float nSigmamin, float nSigmamax) : mnumPart{(int)numPart}, mcutPID{(int)cutPID}, mcutDetector{(int)cutDetector}, mcutType{(int)cutType}, mcutApply{(int)cutApply}, mptMin{ptMin}, mptMax{ptMax}, mdetValuemin{nSigmamin}, mdetValuemax{nSigmamax}
{
}

DGPIDCut::DGPIDCut(float* cutValues)
{
  mnumPart = (int)cutValues[0];
  mcutPID = (int)cutValues[1];
  mcutDetector = (int)cutValues[2];
  mcutType = (int)cutValues[3];
  mcutApply = (int)cutValues[4];
  mptMin = cutValues[5];
  mptMax = cutValues[6];
  mdetValuemin = cutValues[7];
  mdetValuemax = cutValues[8];
}

DGPIDCut::~DGPIDCut()
{
}

// -----------------------------------------------------------------------------
void DGPIDCut::Print()
{
  LOGF(info, "      Cut");
  LOGF(info, "        Part:        %i", mnumPart);
  LOGF(info, "        PID:         %i", mcutPID);
  LOGF(info, "        Detector:    %i", mcutDetector);
  LOGF(info, "        Type:        %i", mcutType);
  LOGF(info, "        Application: %i", mcutApply);
  LOGF(info, "        ptMin:       %f", mptMin);
  LOGF(info, "        ptMax:       %f", mptMax);
  LOGF(info, "        nSigmaMin:   %f", mdetValuemin);
  LOGF(info, "        nSigmaMax:   %f", mdetValuemax);
}

// =============================================================================
// DGPIDCuts
DGPIDCuts::DGPIDCuts()
{
  clear();
}

DGPIDCuts::DGPIDCuts(std::vector<float> PIDCutValues)
{
  setPIDCuts(PIDCutValues);
}

DGPIDCuts::~DGPIDCuts()
{
  clear();
}

// -----------------------------------------------------------------------------
void DGPIDCuts::Print()
{
  LOGF(info, "    Cuts");
  for (auto cut : mDGPIDCuts) {
    cut.Print();
  }
}

// -----------------------------------------------------------------------------
void DGPIDCuts::setPIDCuts(std::vector<float> PIDCutValues)
{
  mDGPIDCuts.clear();

  // check correct number of PIDCutValues
  if ((PIDCutValues.size() % numDGPIDCutParameters) != 0) {
    LOGF(error, "Number of PIDCutValues should be a multiple of %i, but it is %i", numDGPIDCutParameters, PIDCutValues.size());
  }

  // fill mDGPIDCuts
  auto nCuts = PIDCutValues.size() / numDGPIDCutParameters;
  for (uint ind = 0; ind < nCuts; ind++) {
    mDGPIDCuts.push_back(DGPIDCut(&PIDCutValues[ind * numDGPIDCutParameters]));
  }
}

// =============================================================================
// DGAnaparHolder
DGAnaparHolder::DGAnaparHolder()
{
}

DGAnaparHolder::DGAnaparHolder(int nCombine, std::vector<float> DGPIDs, std::vector<float> DGPIDCutValues) : mNCombine{nCombine}, mDGPIDs{DGPIDs}, mDGPIDCutValues{DGPIDCutValues}
{
  makeUniquePermutations();
}

DGAnaparHolder::~DGAnaparHolder()
{
  mDGPIDs.clear();
  mDGPIDCutValues.clear();
  muniquePerms.clear();
}

// -----------------------------------------------------------------------------
void DGAnaparHolder::Print()
{
  LOGF(info, "  DGAnaparHolder");
  LOGF(info, "    nCombine: %i", mNCombine);
  LOGF(info, "    PIDs");
  for (auto pid : mDGPIDs) {
    LOGF(info, "      %f", pid);
  }
  PIDCuts().Print();
}

// -----------------------------------------------------------------------------
DGPIDCuts DGAnaparHolder::PIDCuts()
{
  return DGPIDCuts(mDGPIDCutValues);
}

// -----------------------------------------------------------------------------
void DGAnaparHolder::makeUniquePermutations()
{
  // all permutations of mNCombine elements
  std::vector<std::vector<uint>> perms;
  permutations(mNCombine, perms);

  // compute unique permutations
  std::hash<std::string> hasher;
  std::vector<std::size_t> hashes;
  std::vector<int> perminfo(mNCombine);
  std::string hashstr;

  int cnt;
  for (auto perm : perms) {
    cnt = -1;
    for (auto ind : perm) {
      cnt++;
      perminfo[cnt] = (int)mDGPIDs[ind];
    }
    hashstr = "";
    for (auto tok : perminfo) {
      hashstr += std::to_string(tok);
    }

    // update muniquePerms
    auto hash = hasher(std::string(hashstr));
    if (std::find(hashes.begin(), hashes.end(), hash) == hashes.end()) {
      hashes.push_back(hash);
      for (auto ii = 0; ii < mNCombine; ii++) {
        muniquePerms.push_back(perm[ii]);
      }
    }
  }
}

// -----------------------------------------------------------------------------
// return unique permutations
std::vector<int> DGAnaparHolder::uniquePermutations()
{
  // create unique permutations if not done already
  if (muniquePerms.size() == 0) {
    makeUniquePermutations();
    LOGF(info, "Number of unique permutations: %i", muniquePerms.size() / mNCombine);
  }

  return muniquePerms;
}

// -----------------------------------------------------------------------------
// find all permutations of n0 elements
void DGAnaparHolder::permutations(std::vector<uint>& ref, int n0, int np, std::vector<std::vector<uint>>& perms)
{

  // create local reference
  auto ref2u = ref;

  // loop over np-1 rotations of last np elements of ref
  for (auto ii = 0; ii < np; ii++) {

    // create a new permutation
    // copy first n0-np elements from ref
    // then rotate last np elements of ref
    std::vector<uint> perm(n0, 0);
    for (auto ii = 0; ii < n0 - np; ii++) {
      perm[ii] = ref2u[ii];
    }
    for (auto ii = n0 - np + 1; ii < n0; ii++) {
      perm[ii - 1] = ref2u[ii];
    }
    perm[n0 - 1] = ref2u[n0 - np];

    // add new permutation to the list of permuutations
    if (ii < (np - 1)) {
      perms.push_back(perm);
    }

    // if np>2 then do permutation of next level
    // use the new combination as reference
    if (np > 2) {
      auto newnp = np - 1;
      permutations(perm, n0, newnp, perms);
    }

    // update reference
    ref2u = perm;
  }
}

//-----------------------------------------------------------------------------
// find all permutations of n0 elements
int DGAnaparHolder::permutations(int n0, std::vector<std::vector<uint>>& perms)
{
  // initialize with first trivial combination
  perms.clear();
  if (n0 == 0) {
    return 0;
  }

  std::vector<uint> ref(n0, 0);
  for (auto ii = 0; ii < n0; ii++) {
    ref[ii] = ii;
  }
  perms.push_back(ref);

  // iterate recursively
  permutations(ref, n0, n0, perms);

  return perms.size();
}

// =============================================================================
// DGParticle
DGParticle::DGParticle()
{
}

DGParticle::DGParticle(TDatabasePDG* pdg, DGAnaparHolder anaPars, UDTracksFull const& tracks, std::vector<uint> comb)
{
  // compute invariant mass
  TLorentzVector lvtmp;
  auto pids = anaPars.PIDs();

  // loop over tracks and update mIVM
  mIVM = TLorentzVector(0., 0., 0., 0.);
  auto cnt = -1;
  for (auto ind : comb) {
    cnt++;
    auto track = tracks.rawIteratorAt(ind);
    lvtmp.SetXYZM(track.px(), track.py(), track.pz(), particleMass(pdg, pids[cnt]));
    mIVM += lvtmp;
  }

  // set array of track indices
  mtrkinds = comb;
}

DGParticle::~DGParticle()
{
  mtrkinds.clear();
}

// -----------------------------------------------------------------------------
void DGParticle::Print()
{
  LOGF(info, "DGParticle:");
  LOGF(info, "  Number of particles: %i", mtrkinds.size());
  LOGF(info, "  Mass / pt: %f / %f", mIVM.M(), mIVM.Perp());
  LOGF(info, "");
}

// =============================================================================
// DGPIDSelector
DGPIDSelector::DGPIDSelector()
{
  fPDG = TDatabasePDG::Instance();
}

DGPIDSelector::~DGPIDSelector()
{
  mIVMs.clear();
}

// -----------------------------------------------------------------------------
void DGPIDSelector::Print()
{
  LOGF(info, "PIDSelector");
  mAnaPars.Print();
}

void DGPIDSelector::init(DGAnaparHolder anaPars)
{
  mAnaPars = anaPars;
  mIVMs.clear();
}

// -----------------------------------------------------------------------------
float DGPIDSelector::getTPCnSigma(UDTrackFull track, int pid)
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

// -----------------------------------------------------------------------------
float DGPIDSelector::getTOFnSigma(UDTrackFull track, int pid)
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

// -----------------------------------------------------------------------------
bool DGPIDSelector::isGoodTrack(UDTrackFull track, int cnt)
{
  // get pid of particle cnt
  auto pid = mAnaPars.PIDs()[cnt];

  // unknown PID
  auto pidhypo = pid2ind(pid);
  if (pidhypo < 0) {
    return false;
  }

  // check sign
  if (pid != 0) {
    auto ch = pid / abs(pid);
    if (track.sign() != ch) {
      return false;
    }
  }

  // loop over all PIDCuts and apply the ones which apply to this track
  auto pidcuts = mAnaPars.PIDCuts().Cuts();
  for (auto pidcut : pidcuts) {

    // skip cut if it does not apply to this track
    LOGF(debug, "nPart %i %i, Type %i Apply %i", pidcut.nPart(), cnt, pidcut.cutType(), pidcut.cutApply());
    if (pidcut.nPart() != cnt || pidcut.cutApply() <= 0) {
      continue;
    }

    // check pt
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
      LOGF(info, "detValue TPC %f", detValue);
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
      LOGF(info, "detValue TOF %f", detValue);
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

// -----------------------------------------------------------------------------
int DGPIDSelector::computeIVMs(UDTracksFull const& tracks)
{
  // reset
  mIVMs.clear();

  // create combinations including permutations
  auto combs = combinations(tracks.size());

  // loop over unique combinations
  for (auto comb : combs) {
    // is tracks compatible with PID requirements?
    bool isGoodComb = true;
    auto cnt = -1;
    for (auto ind : comb) {
      cnt++;
      if (!isGoodTrack(tracks.rawIteratorAt(ind), cnt)) {
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

// -----------------------------------------------------------------------------
int DGPIDSelector::pid2ind(int pid)
{
  switch (abs(pid)) {
    case 11: // electron
      return 0;
    case 211: // pion
      return 1;
    case 13: // muon
      return 2;
    case 321: // kaon
      return 3;
    case 2212: // proton
      return 4;
    default: // unknown
      return -1.;
  }
};

// -----------------------------------------------------------------------------
// find selections of np out of n0
void DGPIDSelector::combinations(int n0, std::vector<uint>& pool, int np, std::vector<uint>& inds, int n,
                                 std::vector<std::vector<uint>>& combs)
{
  // loop over pool
  for (auto ii = 0; ii < n0 - n; ii++) {

    inds[n] = pool[ii];

    // if all inds are defined then print them out
    // else get next inds
    if (np == 1) {

      std::vector<uint> comb(n + 1, 0);
      for (uint ii = 0; ii < inds.size(); ii++) {
        comb[ii] = inds[ii];
      }
      combs.push_back(comb);

    } else {

      auto n0new = n0 - ii;
      std::vector<uint> newpool(n0new, 0);
      for (auto kk = 0; kk < n0new; kk++) {
        newpool[kk] = pool[kk + ii + 1];
      }

      auto npnew = np - 1;
      auto nnew = n + 1;
      combinations(n0new, newpool, npnew, inds, nnew, combs);
    }
  }
}

// -----------------------------------------------------------------------------
// find all possible selections of np out of n0
int DGPIDSelector::combinations(int n0, int np, std::vector<std::vector<uint>>& combs)
{
  // initialisations
  combs.clear();
  if (n0 < np) {
    return 0;
  }

  std::vector<uint> pool(n0, 0);
  for (auto ii = 0; ii < n0; ii++) {
    pool[ii] = ii;
  }
  std::vector<uint> inds(np, 0);

  // iterate recursively
  combinations(n0, pool, np, inds, 0, combs);

  return combs.size();
}

// -----------------------------------------------------------------------------
std::vector<std::vector<uint>> DGPIDSelector::combinations(int nPool)
{
  // get the unique permutations
  auto uniquePerms = mAnaPars.uniquePermutations();
  auto numUniquePerms = uniquePerms.size() / mAnaPars.nCombine();

  // all selections of nCombine elements from nPool elements
  std::vector<std::vector<uint>> combs;
  combinations(nPool, mAnaPars.nCombine(), combs);

  // permute the combinations
  std::vector<std::vector<uint>> copes;
  for (auto comb : combs) {
    for (auto ii = 0u; ii < numUniquePerms; ii++) {
      std::vector<uint> cope(mAnaPars.nCombine(), 0);
      for (auto jj = 0; jj < mAnaPars.nCombine(); jj++) {
        auto ind = ii * mAnaPars.nCombine() + jj;
        cope[uniquePerms[ind]] = comb[jj];
      }
      copes.push_back(cope);
    }
  }

  // print copes
  LOGF(debug, "copes");
  for (auto cope : copes) {
    LOGF(debug, "  cope");
    for (auto ind : cope) {
      LOGF(debug, "    %i", ind);
    }
  }

  return copes;
}

// -----------------------------------------------------------------------------
