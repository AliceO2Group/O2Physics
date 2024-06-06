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



template <typename pdgDatabase>
struct ParticleCounter {
  // Service<o2::framework::O2DatabasePDG> mPdgDatabase;
  // ParticleCounter(pdgDatabase& pdgDB) { mPdgDatabase = &pdgDB; }
  bool mSelectPrimaries = true;

  pdgDatabase* mPdgDatabase;

  template <float etamin, float etamax>
  float countMultInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0;
    for (const auto& particle : mcParticles) {

      // primary
      if (mSelectPrimaries && !particle.isPhysicalPrimary()) {
        continue;
      }
      // if ((particle.vx() * particle.vx() + particle.vy() * particle.vy()) > maxProdRadius * maxProdRadius) {
      //   return false;
      // }

      // has pdg
      TParticlePDG* p = mPdgDatabase->get()->GetParticle(particle.pdgCode());
      if (!p) {
        continue;
      }
      // is charged
      if (abs(p->Charge()) == 0) {
        continue;
      }
      // in acceptance
      if (particle.eta() > etamin && particle.eta() < etamax) {
        counter++;
      }
    }
    return counter;
  }

  template <float etamin, float etamax, bool requireNeutral = false>
  float countEnergyInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0.f;
    for (const auto& particle : mcParticles) {

      // primary
      if (mSelectPrimaries && !particle.isPhysicalPrimary()) {
        continue;
      }
      // has pdg
      TParticlePDG* p = mPdgDatabase->get()->GetParticle(particle.pdgCode());
      if (!p) {
        continue;
      }
      // is neutral
      if constexpr (requireNeutral) {
        if (abs(p->Charge()) > 1e-3)
          continue;
      } else {
        if (abs(p->Charge()) <= 1e-3)
          continue;
      }
      // in acceptance
      if (particle.eta() > etamin && particle.eta() < etamax) {
        counter += particle.e();
      }
    }
    return counter;
  }

  float countFT0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<3.5f, 4.9f>(mcParticles); }
  float countFT0C(const aod::McParticles& mcParticles) { return countMultInAcceptance<-3.3f, -2.1f>(mcParticles); }
  float countFV0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<2.2f, 5.1f>(mcParticles); }
  float countFDDA(const aod::McParticles& mcParticles) { return countMultInAcceptance<4.9f, 6.3f>(mcParticles); }
  float countFDDC(const aod::McParticles& mcParticles) { return countMultInAcceptance<-7.f, -4.9f>(mcParticles); }

  float countZNA(const aod::McParticles& mcParticles) { return countEnergyInAcceptance<8.8f, 100.f>(mcParticles); }
  float countZNC(const aod::McParticles& mcParticles) { return countEnergyInAcceptance<-100.f, -8.8f>(mcParticles); }

  float countITSIB(const aod::McParticles& mcParticles) { return countMultInAcceptance<-2.f, 2.f>(mcParticles); }
};



} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_INELGT_H_
