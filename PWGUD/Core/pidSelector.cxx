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

#include "TLorentzVector.h"
#include "CommonConstants/PhysicsConstants.h"
#include "pidSelector.h"

// -----------------------------------------------------------------------------
int pid2ind(int pid)
{
  switch (abs(pid)) {
    case 21: // electron
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
float particleMass(int pid)
{
  switch (abs(pid)) {
    case 21:
      return constants::physics::MassElectron;
    case 211:
      return constants::physics::MassPionCharged;
    case 13:
      return constants::physics::MassMuon;
    case 321:
      return constants::physics::MassKaonCharged;
    case 2212:
      return constants::physics::MassProton;
    default:
      return -1.;
  }
};

// -----------------------------------------------------------------------------
// find all permutations of n0 elements
void permutations(std::vector<uint>& ref, int n0, int np, std::vector<std::vector<uint>>& perms)
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
int permutations(int n0, std::vector<std::vector<uint>>& perms)
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

//-----------------------------------------------------------------------------
// find selections of np out of n0
void combinations(int n0, std::vector<uint>& pool, int np, std::vector<uint>& inds, int n,
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
int combinations(int n0, int np, std::vector<std::vector<uint>>& combs)
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
std::vector<std::vector<uint>> combinations(int nCombine, int nPool)
{
  // all permutations of nCombine elements
  std::vector<std::vector<uint>> perms;
  permutations(nCombine, perms);

  // all selections of nCombine elements from nPool elements
  std::vector<std::vector<uint>> combs;
  combinations(nPool, nCombine, combs);

  // permute the combinations
  std::vector<std::vector<uint>> copes;
  for (auto comb : combs) {
    for (auto perm : perms) {
      std::vector<uint> cope(nCombine, 0);
      for (auto ii = 0; ii < nCombine; ii++) {
        cope[perm[ii]] = comb[ii];
      }
      copes.push_back(cope);
    }
  }

  return copes;
}

// -----------------------------------------------------------------------------
DGParticle::DGParticle(anaparHolder anaPars, aod::DGTracks const& dgtracks, std::vector<uint> comb)
{
  // compute invariant mass
  TLorentzVector lvtmp, IVM;
  auto pidinfo = anaPars.TPCnSigmas();

  // loop over tracks
  auto cnt = -1;
  for (auto ind : comb) {
    cnt++;
    auto track = dgtracks.rawIteratorAt(ind);
    lvtmp.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), particleMass(pidinfo[cnt * 12]));
    IVM += lvtmp;
  }

  mM = IVM.M();
  mPerp = IVM.Perp();
  mtrkinds = comb;
}

// -----------------------------------------------------------------------------
void DGParticle::Print()
{
  LOGF(info, "DGParticle:");
  LOGF(info, "  Number of particles: %i", mtrkinds.size());
  LOGF(info, "  Mass / pt: %f / %f", mM, mPerp);
  LOGF(info, "");
}

// -----------------------------------------------------------------------------
float pidSelector::getTPCnSigma(aod::DGTrack track, int hypo)
{
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
      return 999.;
  }
}

// -----------------------------------------------------------------------------
bool pidSelector::isGoodTrack(aod::DGTrack track, int cnt)
{
  // extract PID information
  auto pidinfo = mAnaPars.TPCnSigmas();

  // get pid of particle cnt
  auto ind = cnt * 12;
  auto pidhypo = pid2ind(pidinfo[ind]);

  // check sign
  if (pidinfo[ind + 1] != 0 && track.sign() != pidinfo[ind + 1]) {
    return false;
  }

  // any PID
  if (pidhypo < 0) {
    return true;
  }

  // check nSigma
  for (auto hypo = 0; hypo < 5; hypo++) {
    auto nSigma = getTPCnSigma(track, hypo);
    ind += 2;
    if (pidinfo[ind] == 0. && pidinfo[ind + 1] == 0.) {
      continue;
    }
    if (hypo == pidhypo) {
      // inclusive limits
      if (nSigma < pidinfo[ind] || nSigma > pidinfo[ind + 1]) {
        return false;
      }
    } else {
      // exclusive limits
      if (nSigma > pidinfo[ind] && nSigma < pidinfo[ind + 1]) {
        return false;
      }
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
int pidSelector::computeIVMs(int nCombine, aod::DGTracks const& dgtracks)
{
  // reset
  mIVMs.clear();

  // create combinations including permutations
  auto combs = combinations(nCombine, dgtracks.size());

  // loop over combinations
  for (auto comb : combs) {
    // is dgtracks compatible with PID requirements?
    bool isGoodComb = true;
    auto cnt = -1;
    for (auto ind : comb) {
      cnt++;
      if (!isGoodTrack(dgtracks.rawIteratorAt(ind), cnt)) {
        isGoodComb = false;
        break;
      }
    }

    // update list of IVMs
    if (isGoodComb) {
      DGParticle IVM(mAnaPars, dgtracks, comb);
      mIVMs.push_back(IVM);
    }
  }

  return mIVMs.size();
}

// -----------------------------------------------------------------------------
