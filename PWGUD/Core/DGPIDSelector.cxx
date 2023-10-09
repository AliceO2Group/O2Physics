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
                   float ptMin, float ptMax, float nSigmamin, float nSigmamax) : mnumPart{static_cast<int>(numPart)}, mcutPID{static_cast<int>(cutPID)}, mcutDetector{static_cast<int>(cutDetector)}, mcutType{static_cast<int>(cutType)}, mcutApply{static_cast<int>(cutApply)}, mptMin{ptMin}, mptMax{ptMax}, mdetValuemin{nSigmamin}, mdetValuemax{nSigmamax}
{
}

DGPIDCut::DGPIDCut(float* cutValues)
{
  mnumPart = static_cast<int>(cutValues[0]);
  mcutPID = static_cast<int>(cutValues[1]);
  mcutDetector = static_cast<int>(cutValues[2]);
  mcutType = static_cast<int>(cutValues[3]);
  mcutApply = static_cast<int>(cutValues[4]);
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
  LOGF(info, "    min number tracks: %d", mMinNTracks);
  LOGF(info, "    max number tracks: %d", mMaxNTracks);
  LOGF(info, "    min fraction of PV contr. with TOF: %f", mMinRgtrwTOF);
  LOGF(info, "    max dcaxy:         %f", mMaxDCAxy);
  LOGF(info, "    max dcaz:          %f", mMaxDCAz);
  LOGF(info, "    min dBC:           %d", mdBCMin);
  LOGF(info, "    max dBC:           %d", mdBCMax);
  LOGF(info, "    FIT vetoes (FV0A, FT0A, FT0C, FDDA, FDDC)");
  LOGF(info, "      %d %d %d %d %d", mFITvetoes[0], mFITvetoes[1], mFITvetoes[2], mFITvetoes[3], mFITvetoes[4]);
  LOGF(info, "    min NCl TPC:       %d", mMinNClTPC);
  LOGF(info, "    max NCl TPC:       %d", mMaxNClTPC);
  LOGF(info, "    min chi^{2} TPC    %f", mMinChi2NClTPC);
  LOGF(info, "    max chi^{2} TPC    %f", mMaxChi2NClTPC);
  LOGF(info, "    min track pT:      %f", mMinpt);
  LOGF(info, "    max track pT:      %f", mMaxpt);
  LOGF(info, "    min eta:           %f", mMineta);
  LOGF(info, "    max eta:           %f", mMaxeta);
  LOGF(info, "    min alpha:         %f", mMinAlpha);
  LOGF(info, "    max alpha:         %f", mMaxAlpha);
  LOGF(info, "    min system pT:     %f", mMinptsys);
  LOGF(info, "    max system pT:     %f", mMaxptsys);
  LOGF(info, "    nCombine:          %d", mNCombine);
  LOGF(info, "    unlike charges");
  for (auto ch : mUnlikeCharges) {
    LOGF(info, "      %i", ch);
  }
  LOGF(info, "    like charges");
  for (auto ch : mLikeCharges) {
    LOGF(info, "      %i", ch);
  }
  LOGF(info, "    PIDs");
  for (auto pid : mDGPIDs) {
    LOGF(info, "      %d", pid);
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
  // reset
  muniquePerms.clear();

  // all permutations of mNCombine elements
  std::vector<std::vector<int>> perms;
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
      perminfo[cnt] = mDGPIDs[ind];
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
  if (muniquePerms.size() < mNCombine) {
    makeUniquePermutations();
  }

  return muniquePerms;
}

// -----------------------------------------------------------------------------
// find all permutations of n0 elements
void DGAnaparHolder::permutations(std::vector<int>& ref, int n0, int np, std::vector<std::vector<int>>& perms)
{
  // create local reference
  auto ref2u = ref;

  // loop over np-1 rotations of last np elements of ref
  for (auto ii = 0; ii < np; ii++) {

    // create a new permutation
    // copy first n0-np elements from ref
    // then rotate last np elements of ref
    std::vector<int> perm(n0, 0);
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
int DGAnaparHolder::permutations(int n0, std::vector<std::vector<int>>& perms)
{
  // initialize with first trivial combination
  perms.clear();
  if (n0 == 0) {
    return 0;
  }

  std::vector<int> ref(n0, 0);
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

DGParticle::~DGParticle()
{
  mtrkinds.clear();
}

// -----------------------------------------------------------------------------
void DGParticle::Print()
{
  LOGF(info, "DGParticle:");
  LOGF(info, "  Number of particles: %i", mtrkinds.size());
  LOGF(info, "  Indices:");
  for (auto ind : mtrkinds) {
    LOGF(info, "    %d", ind);
  }
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
  mUnlikeIVMs.clear();
  mLikeIVMs.clear();
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
  mUnlikeIVMs.clear();
  mLikeIVMs.clear();
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
void DGPIDSelector::combinations(int n0, std::vector<int>& pool, int np, std::vector<int>& inds, int n,
                                 std::vector<std::vector<int>>& combs)
{
  // loop over pool
  for (auto ii = 0; ii < n0 - n; ii++) {

    inds[n] = pool[ii];

    // if all inds are defined then print them out
    // else get next inds
    if (np == 1) {

      std::vector<int> comb(n + 1, 0);
      for (uint ii = 0; ii < inds.size(); ii++) {
        comb[ii] = inds[ii];
      }
      combs.push_back(comb);

    } else {

      auto n0new = n0 - ii;
      std::vector<int> newpool(n0new, 0);
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
int DGPIDSelector::combinations(int n0, int np, std::vector<std::vector<int>>& combs)
{
  // initialisations
  combs.clear();
  if (n0 < np) {
    return 0;
  }

  std::vector<int> pool(n0, 0);
  for (auto ii = 0; ii < n0; ii++) {
    pool[ii] = ii;
  }
  std::vector<int> inds(np, 0);

  // iterate recursively
  combinations(n0, pool, np, inds, 0, combs);

  return combs.size();
}

// -----------------------------------------------------------------------------
std::vector<std::vector<int>> DGPIDSelector::combinations(int nPool)
{
  // get the unique permutations
  auto uniquePerms = mAnaPars.uniquePermutations();
  auto numUniquePerms = uniquePerms.size() / mAnaPars.nCombine();

  // all selections of nCombine elements from nPool elements
  std::vector<std::vector<int>> combs;
  combinations(nPool, mAnaPars.nCombine(), combs);

  // permute the combinations
  std::vector<std::vector<int>> copes;
  for (auto comb : combs) {
    for (auto ii = 0u; ii < numUniquePerms; ii++) {
      std::vector<int> cope(mAnaPars.nCombine(), 0);
      for (auto jj = 0; jj < mAnaPars.nCombine(); jj++) {
        auto ind = ii * mAnaPars.nCombine() + jj;
        cope[uniquePerms[ind]] = comb[jj];
      }
      copes.push_back(cope);
    }
  }

  return copes;
}

// -----------------------------------------------------------------------------
