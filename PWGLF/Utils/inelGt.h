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

///
/// \file  inelGt.h
/// \since 07/05/2024
/// \brief Utilities for counting charged particles in |eta| region and define the INEL > N selection in the MC at the generated level
///        Based on cascqaanalysis.cxx
///

#ifndef PWGLF_UTILS_INELGT_H_
#define PWGLF_UTILS_INELGT_H_

#include <vector>

namespace o2
{
namespace pwglf
{

struct EtaCharge {
  float eta;
  int charge;
};

// Struct for counting charged particles in |eta| region
template <typename TMcParticles, typename pdgDatabase>
bool isINELgtNmc(TMcParticles particles, int nChToSatisfySelection, pdgDatabase pdgDB)
{
  // INEL > N (at least N+1 charged particles in |eta| < 1.0)
  EtaCharge etaCharge;
  std::vector<EtaCharge> ParticlesEtaAndCharge(particles.size());
  unsigned int nParticles = 0;
  for (const auto& particle : particles) {
    if (particle.isPhysicalPrimary() == 0) {
      continue; // consider only primaries
    }
    etaCharge = {999, 0}; // refresh init. for safety
    TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
    if (!p) {
      switch (std::to_string(particle.pdgCode()).length()) {
        case 10: // nuclei
        {
          etaCharge = {particle.eta(), static_cast<int>(particle.pdgCode() / 10000 % 1000)};
          ParticlesEtaAndCharge[nParticles++] = etaCharge;
          break;
        }
        default:
          break;
      }
    } else {
      etaCharge = {particle.eta(), static_cast<int>(p->Charge())};
      ParticlesEtaAndCharge[nParticles++] = etaCharge;
    }
  }

  ParticlesEtaAndCharge.resize(nParticles);

  auto etaChargeConditionFunc = [](EtaCharge elem) {
    return ((TMath::Abs(elem.eta) < 1.0) && (TMath::Abs(elem.charge) > 0.001));
  };

  if (std::count_if(ParticlesEtaAndCharge.begin(), ParticlesEtaAndCharge.end(), etaChargeConditionFunc) > nChToSatisfySelection) {
    return true;
  } else {
    return false;
  }
}

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_INELGT_H_
