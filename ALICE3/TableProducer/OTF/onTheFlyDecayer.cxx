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

/// \file onTheFlyDecayer.cxx
///
/// \brief Task to decay long lived particles and to propagate the information to other tasks
///
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, UniUPO
///

#include "Framework/O2DatabasePDGPlugin.h"


using namespace o2;
using namespace o2::framework;
using std::array;

struct OnTheFlyDecayer {
      Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(o2::framework::InitContext&)
  {
  }

template <typename ParticleType>
  bool decayParticle(const auto & particle){


bool canDecay = false;

switch(particle.pdgCode()){
    case 3312:
    canDecay = true;
}
if(!canDecay){
    return false;
}
// Check that it does not have daughters
if(particle.hasDaughters()){
LOG(fatal) << "Particle has daughters";
}


      const auto& pdgInfo = pdgDB->GetParticle(particle.pdgCode());
      if (!pdgInfo) {
        LOG(fatal) << "PDG code " << particle.pdgCode() << " not found in the database";
      }
      



    const double u = rand.Uniform(0, 1);
    double xi_mass = o2::constants::physics::MassXiMinus;
    double la_mass = o2::constants::physics::MassLambda;
    double pi_mass = o2::constants::physics::MassPionCharged;
    double pr_mass = o2::constants::physics::MassProton;

    double mass = 0.;
    double tau = 0.;
// Compute channel
    switch (particle.pdgCode())
    {
        case 3312:
        mass = xi_mass;
        tau = 4.91;
        break;
    }

    const double gamma = 1 / sqrt(1 + (particle.p() * particle.p()) / (mass * mass));
    const double ctau = tau * gamma;
    const double rxyz = (-ctau * log(1.0 - u));
    // If the particle is charged, then propagate in the mag field
    o2::math_utils::CircleXYf_t circle;
    if (pdgInfo->Charge() != 0) {
    float sna, csa;
    track.getCircleParams(magneticField, circle, sna, csa);
      }
      else{ // Neutral particles

      }
    const double rxy = rxyz / sqrt(1. + track.getTgl() * track.getTgl());
    const double theta = rxy / circle.rC;
    const double newX = ((particle.vx() - circle.xC) * std::cos(theta) - (particle.vy() - circle.yC) * std::sin(theta)) + circle.xC;
    const double newY = ((particle.vy() - circle.yC) * std::cos(theta) + (particle.vx() - circle.xC) * std::sin(theta)) + circle.yC;
    const double newPx = particle.px() * std::cos(theta) - particle.py() * std::sin(theta);
    const double newPy = particle.py() * std::cos(theta) + particle.px() * std::sin(theta);
    xiDecayVertex.push_back(newX);
    xiDecayVertex.push_back(newY);
    xiDecayVertex.push_back(particle.vz() + rxyz * (particle.pz() / particle.p()));

    std::vector<double> xiDaughters = {la_mass, pi_mass};
    TLorentzVector xi(newPx, newPy, particle.pz(), particle.e());
    TGenPhaseSpace xiDecay;
    xiDecay.SetDecay(xi, 2, xiDaughters.data());
    xiDecay.Generate();
    decayDaughters.push_back(*xiDecay.GetDecay(1));
    return true;

    TLorentzVector la = *xiDecay.GetDecay(0);

    double la_gamma = 1 / sqrt(1 + (la.P() * la.P()) / (la_mass * la_mass));
    double la_ctau = 7.89 * la_gamma;
    std::vector<double> laDaughters = {pi_mass, pr_mass};
    double la_rxyz = (-la_ctau * log(1 - u));
    laDecayVertex.push_back(xiDecayVertex[0] + la_rxyz * (xiDecay.GetDecay(0)->Px() / xiDecay.GetDecay(0)->P()));
    laDecayVertex.push_back(xiDecayVertex[1] + la_rxyz * (xiDecay.GetDecay(0)->Py() / xiDecay.GetDecay(0)->P()));
    laDecayVertex.push_back(xiDecayVertex[2] + la_rxyz * (xiDecay.GetDecay(0)->Pz() / xiDecay.GetDecay(0)->P()));

    TGenPhaseSpace laDecay;
    laDecay.SetDecay(la, 2, laDaughters.data());
    laDecay.Generate();
    decayDaughters.push_back(*laDecay.GetDecay(0));
    decayDaughters.push_back(*laDecay.GetDecay(1));

}


  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
for(const auto & particle : mcParticles) {
decayParticle(particle);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
