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

/// \brief cut selection and cut functions for photon candidates
/// \author marvin.hemmer@cern.ch

#include <vector>

#include "Framework/AnalysisTask.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#ifndef PWGEM_PHOTONMESON_UTILS_GAMMASELECTIONCUTS_H_
#define PWGEM_PHOTONMESON_UTILS_GAMMASELECTIONCUTS_H_

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace emccuts
{
// set the standard cuts
std::vector<float> EMC_minTime = {-20.f};
std::vector<float> EMC_maxTime = {+25.f};
std::vector<float> EMC_minM02 = {0.1f};
std::vector<float> EMC_maxM02 = {0.7f};
std::vector<float> EMC_minE = {0.7f};
std::vector<int> EMC_minNCell = {1};
std::vector<std::vector<float>> EMC_TM_Eta = {{0.01f, 4.07f, -2.5f}};
std::vector<std::vector<float>> EMC_TM_Phi = {{0.015f, 3.65f, -2.f}};
std::vector<float> EMC_Eoverp = {1.75f};
} // namespace emccuts
void gatherCutsEMC(float minTime, float maxTime, float minM02, float maxM02, float minE, int minNCell, std::vector<float> TM_Eta, std::vector<float> TM_Phi, float Eoverp)
{
  // insert the configurable values in first position
  emccuts::EMC_minTime.insert(emccuts::EMC_minTime.begin(), minTime);
  emccuts::EMC_maxTime.insert(emccuts::EMC_maxTime.begin(), maxTime);
  emccuts::EMC_minM02.insert(emccuts::EMC_minM02.begin(), minM02);
  emccuts::EMC_maxM02.insert(emccuts::EMC_maxM02.begin(), maxM02);
  emccuts::EMC_minE.insert(emccuts::EMC_minE.begin(), minE);
  emccuts::EMC_minNCell.insert(emccuts::EMC_minNCell.begin(), minNCell);
  emccuts::EMC_TM_Eta.insert(emccuts::EMC_TM_Eta.begin(), TM_Eta);
  emccuts::EMC_TM_Phi.insert(emccuts::EMC_TM_Phi.begin(), TM_Phi);
  emccuts::EMC_Eoverp.insert(emccuts::EMC_Eoverp.begin(), Eoverp);

  // fill up the rest of the vectors to size 64 to ensure no crashes can happen
  emccuts::EMC_minTime.resize(64, 0);
  emccuts::EMC_maxTime.resize(64, 0);
  emccuts::EMC_minM02.resize(64, 0);
  emccuts::EMC_maxM02.resize(64, 0);
  emccuts::EMC_minE.resize(64, 0);
  emccuts::EMC_minNCell.resize(64, 0);
  emccuts::EMC_TM_Eta.resize(64, {0, 0, 0});
  emccuts::EMC_TM_Phi.resize(64, {0, 0, 0});
  emccuts::EMC_Eoverp.resize(64, 0);
}

uint64_t doTimeCutEMC(int iCut, uint64_t cutbit, aod::SkimEMCCluster const& cluster, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) {                                                             // check if current cut should be applied
    if (cluster.time() <= emccuts::EMC_maxTime.at(iCut) && cluster.time() >= emccuts::EMC_minTime.at(iCut)) { // check cut itself
      cut_return |= (1 << iCut);                                                                              // set bit of current cut to 1 for passing the cut
    } else {
      registry.fill(HIST("hCaloCuts_EMC"), 1, iCut);
    }
  }
  return cut_return;
}

uint64_t doM02CutEMC(int iCut, uint64_t cutbit, aod::SkimEMCCluster const& cluster, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) {                                                         // check if current cut should be applied
    if (cluster.m02() <= emccuts::EMC_maxM02.at(iCut) && cluster.m02() >= emccuts::EMC_minM02.at(iCut)) { // check cut itself
      cut_return |= (1 << iCut);                                                                          // set bit of current cut to 1 for passing the cut
    } else {
      registry.fill(HIST("hCaloCuts_EMC"), 2, iCut);
    }
  }
  return cut_return;
}

uint64_t doMinECutEMC(int iCut, uint64_t cutbit, aod::SkimEMCCluster const& cluster, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) {   // check if current cut should be applied
    if (cluster.e() > emccuts::EMC_minE.at(iCut)) { // check cut itself
      cut_return |= (1 << iCut);                    // set bit of current cut to 1 for passing the cut
    } else {
      registry.fill(HIST("hCaloCuts_EMC"), 3, iCut);
    }
  }
  return cut_return;
}

uint64_t doNCellCutEMC(int iCut, uint64_t cutbit, aod::SkimEMCCluster const& cluster, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) {             // check if current cut should be applied
    if (cluster.nCells() >= emccuts::EMC_minNCell.at(iCut)) { // check cut itself
      cut_return |= (1 << iCut);                              // set bit of current cut to 1 for passing the cut
    } else {
      registry.fill(HIST("hCaloCuts_EMC"), 4, iCut);
    }
  }
  return cut_return;
}

uint64_t doTrackMatchingEMC(int iCut, uint64_t cutbit, aod::SkimEMCCluster const& cluster, aod::SkimEMCMTs const& tracks, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  double dEta, dPhi;
  if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) { // check if current cut should be applied
    bool hasMatchedTrack_EMC = false;
    for (const auto& track : tracks) {
      dEta = track.tracketa() - cluster.eta();
      dPhi = track.trackphi() - cluster.phi();
      if ((fabs(dEta) <= emccuts::EMC_TM_Eta.at(iCut).at(0) + pow(track.trackpt() + emccuts::EMC_TM_Eta.at(iCut).at(1), emccuts::EMC_TM_Eta.at(iCut).at(2))) &&
          (fabs(dPhi) <= emccuts::EMC_TM_Phi.at(iCut).at(0) + pow(track.trackpt() + emccuts::EMC_TM_Phi.at(iCut).at(1), emccuts::EMC_TM_Phi.at(iCut).at(2))) &&
          cluster.e() / track.trackp() < emccuts::EMC_Eoverp.at(iCut)) { // check cut itself
        hasMatchedTrack_EMC = true;                                      // set bit of current cut to 1 for passing the cut
      }
    }
    if (hasMatchedTrack_EMC) {
      registry.fill(HIST("hCaloCuts_EMC"), 5, iCut);
    } else {
      cut_return |= (1 << iCut);
    }
  }
  return cut_return;
}

uint64_t doPhotonCutsEMC(uint64_t cutbit, aod::SkimEMCCluster const& cluster, aod::SkimEMCMTs const& tracks, Preslice<o2::aod::SkimEMCMTs> perEMCClusterMT, HistogramRegistry& registry)
{
  uint64_t cut_return = 0;
  auto tracksMatchedEMC = tracks.sliceBy(perEMCClusterMT, cluster.globalIndex());
  for (int iCut = 0; iCut < 64; iCut++) {           // loop over max number of cut settings
    if (cutbit & ((uint64_t)1 << (uint64_t)iCut)) { // check each cut setting if it is selected
      registry.fill(HIST("hClusterEIn"), cluster.e(), iCut);
      cut_return = doTimeCutEMC(iCut, cutbit, cluster, registry);
      // use cut_return instead of cutbit from here on to only check cut settings that we want to look at
      // where the cluster did not fail in the previous cut(s)
      cut_return = doM02CutEMC(iCut, cut_return, cluster, registry);
      cut_return = doMinECutEMC(iCut, cut_return, cluster, registry);
      cut_return = doNCellCutEMC(iCut, cut_return, cluster, registry);
      cut_return = doTrackMatchingEMC(iCut, cut_return, cluster, tracksMatchedEMC, registry);

      registry.fill(HIST("hCaloCuts_EMC"), 0, iCut);
      if (cut_return & ((uint64_t)1 << (uint64_t)iCut)) {
        registry.fill(HIST("hClusterEOut"), cluster.e(), iCut);
        registry.fill(HIST("hCaloCuts_EMC"), 6, iCut);
      }
    }
  }
  return cut_return;
}

#endif // PWGEM_PHOTONMESON_UTILS_GAMMASELECTIONCUTS_H_
