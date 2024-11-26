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

/// \file utilsCorrelations.h
/// \brief Utilities for HFC analyses
/// \author Xu Wang <waxu@cern.ch>

#ifndef PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
#define PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_

#include <cmath>
#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"

namespace o2::analysis::hf_correlations
{
enum Region {
  Default = 0,
  Toward,
  Away,
  Transverse
};

template <typename T>
Region getRegion(T const deltaPhi)
{
  if (std::abs(deltaPhi) < o2::constants::math::PIThird) {
    return Toward;
  } else if (deltaPhi > 2. * o2::constants::math::PIThird && deltaPhi < 4. * o2::constants::math::PIThird) {
    return Away;
  } else {
    return Transverse;
  }
}
//apply PID selection on associated tracks
template <typename Atrack, typename SpeciesContainer, typename T1, typename T2>
bool passPIDSelection(Atrack const& track, SpeciesContainer const& mPIDspecies, 
                      T1 const maxTPC, T2 const maxTOF, double ptThreshold = 0.75) {
    size_t speciesIndex = 0;

    // Ensure size consistency
    if (mPIDspecies.size() != maxTPC.value.size() || mPIDspecies.size() != maxTOF.value.size()) {
        LOGF(error, "Size of particle species and corresponding nSigma selection array should be same");
    }

    for (auto const& pid : mPIDspecies) {
        auto nSigmaTPC = o2::aod::pidutils::tpcNSigma(pid, track);
        if (std::abs(nSigmaTPC) > maxTPC->at(speciesIndex)) {
            return false;
        }

        if (track.pt() > ptThreshold && track.hasTOF()) {
            auto nSigmaTOF = o2::aod::pidutils::tofNSigma(pid, track);
            if (std::abs(nSigmaTOF) > maxTOF->at(speciesIndex)) {
                return false;
            }
        }
        //std::cout<<"species "<<pid<<"  number i "<<speciesIndex<<"  tpcSigmaMax  "<<maxTPC->at(speciesIndex)<<"  tofSigmaMax  "<<maxTOF->at(speciesIndex)<<std::endl;
        //std::cout<<"proton nSigma TPC "<<track.tpcNSigmaPr()<<" "<<nSigmaTPC<<"  tof "<<track.tofNSigmaPr()<<std::endl;
        ++speciesIndex;
    }
    return true;
}

//template <typename Atrack, typename species, typename T1, typename T2>
//void passPIDSelction(Atrack const& track, species const mPIDspecies, T1 const maxTPC, T2  maxTOF){
//   int iSpecie =0;
//          for (auto pid : mPIDspecies){
//          auto nSigmaTPC = o2::aod::pidutils::tpcNSigma(pid, track);
//          auto nSigmaTOF = o2::aod::pidutils::tofNSigma(pid, track);
//          if(std::abs(nSigmaTPC) > maxTPC->at(iSpecie)) continue;
//          if(track.pt() > 0.75 && track.hasTOF()){
//            if(std::abs(nSigmaTOF) > maxTOF->at(iSpecie)) continue;
//          }
//          iSpecie++;
//        }
//}
// ========= Find Leading Particle ==============
template <typename TTracks, typename T1, typename T2, typename T3>
int findLeadingParticle(TTracks const& tracks, T1 const dcaXYTrackMax, T2 const dcaZTrackMax, T3 const etaTrackMax)
{
  auto leadingParticle = tracks.begin();
  for (auto const& track : tracks) {
    if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
      continue;
    }
    if (std::abs(track.eta()) > etaTrackMax) {
      continue;
    }
    if (track.pt() > leadingParticle.pt()) {
      leadingParticle = track;
    }
  }
  return leadingParticle.globalIndex();
}

// ======= Find Leading Particle for McGen ============
template <typename TMcParticles, typename T1, typename T2>
int findLeadingParticleMcGen(TMcParticles const& mcParticles, T1 const etaTrackMax, T2 const ptTrackMin)
{
  auto leadingParticle = mcParticles.begin();
  for (auto const& mcParticle : mcParticles) {
    if (std::abs(mcParticle.eta()) > etaTrackMax) {
      continue;
    }
    if (mcParticle.pt() < ptTrackMin) {
      continue;
    }
    if ((std::abs(mcParticle.pdgCode()) != kElectron) && (std::abs(mcParticle.pdgCode()) != kMuonMinus) && (std::abs(mcParticle.pdgCode()) != kPiPlus) && (std::abs(mcParticle.pdgCode()) != kKPlus) && (std::abs(mcParticle.pdgCode()) != kProton)) {
      continue;
    }
    if (mcParticle.pt() > leadingParticle.pt()) {
      leadingParticle = mcParticle;
    }
  }
  return leadingParticle.globalIndex();
}
} // namespace o2::analysis::hf_correlations
#endif // PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
