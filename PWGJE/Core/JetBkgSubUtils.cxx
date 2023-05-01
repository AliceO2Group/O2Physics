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

// jet finder task
//
// Author: Hadi Hassan, Universiy of Jväskylä, hadi.hassan@cern.ch
#include <memory>
#include <tuple>
#include "Framework/Logger.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/FastJetUtilities.h"

JetBkgSubUtils::JetBkgSubUtils(float jetBkgR_out, float constSubAlpha_out, float constSubRMax_out, float bkgEtaMin_out, float bkgEtaMax_out, float bkgPhiMin_out, float bkgPhiMax_out, fastjet::GhostedAreaSpec ghostAreaSpec_out, int nHardReject) : jetBkgR(jetBkgR_out),
                                                                                                                                                                                                                                                      constSubAlpha(constSubAlpha_out),
                                                                                                                                                                                                                                                      constSubRMax(constSubRMax_out),
                                                                                                                                                                                                                                                      bkgEtaMin(bkgEtaMin_out),
                                                                                                                                                                                                                                                      bkgEtaMax(bkgEtaMax_out),
                                                                                                                                                                                                                                                      bkgPhiMin(bkgPhiMin_out),
                                                                                                                                                                                                                                                      bkgPhiMax(bkgPhiMax_out),
                                                                                                                                                                                                                                                      ghostAreaSpec(ghostAreaSpec_out)

{
  // Note: if you are using the PerpCone method you should jetBkgR to be the same as the anit_kt jets R, otherwise use R=0.2
  jetDefBkg = fastjet::JetDefinition(algorithmBkg, jetBkgR, recombSchemeBkg, fastjet::Best);
  areaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghostAreaSpec);
  selRho = fastjet::SelectorRapRange(bkgEtaMin, bkgEtaMax) && fastjet::SelectorPhiRange(bkgPhiMin, bkgPhiMax) && !fastjet::SelectorNHardest(nHardReject); // here we have to put rap range, to be checked!
}

std::tuple<double, double> JetBkgSubUtils::estimateRhoAreaMedian(const std::vector<fastjet::PseudoJet>& inputParticles)
{

  if (inputParticles.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  // cluster the kT jets
  fastjet::ClusterSequenceArea clusterSeq(inputParticles, jetDefBkg, areaDefBkg);

  // select jets in detector acceptance
  std::vector<fastjet::PseudoJet> alljets = selRho(clusterSeq.inclusive_jets());

  double totaljetAreaPhys(0), totalAreaCovered(0);
  std::vector<double> rhovector;
  std::vector<double> rhoMvector;

  // Fill a vector for pT/area to be used for the median
  for (auto& ijet : alljets) {

    // Physical area/ Physical jets (no ghost)
    if (!clusterSeq.is_pure_ghost(ijet)) {
      rhovector.push_back(ijet.perp() / ijet.area());
      rhoMvector.push_back(getM(ijet) / ijet.area());

      totaljetAreaPhys += ijet.area();
    }
    // Full area
    totalAreaCovered += ijet.area();
  }
  // calculate Rho as the median of the jet pT / jet area
  double rho = TMath::Median<double>(rhovector.size(), rhovector.data());
  double rhoM = TMath::Median<double>(rhoMvector.size(), rhoMvector.data());

  if (doSparseSub) {
    // calculate The ocupancy factor, which the ratio of covered area / total area
    double occupancyFactor = totalAreaCovered > 0 ? totaljetAreaPhys / totalAreaCovered : 1.;
    rho *= occupancyFactor;
    rhoM *= occupancyFactor;
  }

  return std::make_tuple(rho, rhoM);
}

std::tuple<double, double> JetBkgSubUtils::estimateRhoPerpCone(const std::vector<fastjet::PseudoJet>& inputParticles, const std::vector<fastjet::PseudoJet>& jets)
{

  if (inputParticles.size() == 0 || jets.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  double perpPtDensity1 = 0;
  double perpPtDensity2 = 0;
  double perpMtDensity1 = 0;
  double perpMtDensity2 = 0;

  fastjet::Selector selectJet = fastjet::SelectorEtaRange(bkgEtaMin, bkgEtaMax) && fastjet::SelectorPhiRange(bkgPhiMin, bkgPhiMax);

  std::vector<fastjet::PseudoJet> selectedJets = fastjet::sorted_by_pt(selectJet(jets));

  if (selectedJets.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  fastjet::PseudoJet leadingJet = selectedJets[0];

  double dPhi1 = 999.;
  double dPhi2 = 999.;
  double dEta = 999.;
  double PerpendicularConeAxisPhi1 = 999., PerpendicularConeAxisPhi2 = 999.;
  // build 2 perp cones in phi around the leading jet (right and left of the jet)
  PerpendicularConeAxisPhi1 = RecoDecay::constrainAngle<double, double>(leadingJet.phi() + (M_PI / 2.)); // This will contrain the angel between 0-2Pi
  PerpendicularConeAxisPhi2 = RecoDecay::constrainAngle<double, double>(leadingJet.phi() - (M_PI / 2.)); // This will contrain the angel between 0-2Pi

  for (auto& particle : inputParticles) {
    // sum the momentum of all paricles that fill the two cones
    dPhi1 = particle.phi() - PerpendicularConeAxisPhi1;
    dPhi1 = RecoDecay::constrainAngle<double, double>(dPhi1, -M_PI);
    dPhi2 = particle.phi() - PerpendicularConeAxisPhi2;
    dPhi2 = RecoDecay::constrainAngle<double, double>(dPhi2, -M_PI);
    dEta = leadingJet.eta() - particle.eta();
    if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= jetBkgR) {
      perpPtDensity1 += particle.perp();
      perpMtDensity1 += particle.mt();
    }

    if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= jetBkgR) {
      perpPtDensity2 += particle.perp();
      perpMtDensity2 += particle.mt();
    }
  }

  // Caculate rho as the ratio of average pT of the two cones / the cone area
  double perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * M_PI * jetBkgR * jetBkgR);
  double perpMtDensity = (perpMtDensity1 + perpMtDensity2) / (2 * M_PI * jetBkgR * jetBkgR);

  return std::make_tuple(perpPtDensity, perpMtDensity);
}

/// Sets the background subtraction object
fastjet::Subtractor JetBkgSubUtils::setSub(std::vector<fastjet::PseudoJet>& inputParticles, std::vector<fastjet::PseudoJet>& jets, float& rhoParam, float& rhoMParam, BkgSubEstimator bkgSubEst, BkgSubMode bkgSubMode)
{

  if (bkgSubEst == BkgSubEstimator::medianRhoSparse) {
    doSparseSub = true;
  }

  fastjet::Subtractor sub;
  fastjet::Selector selRemoveHFCand = !FastJetUtilities::SelectorIsHFCand();

  if (bkgSubEst == BkgSubEstimator::medianRho || bkgSubEst == BkgSubEstimator::medianRhoSparse) {

    std::tie(rhoParam, rhoMParam) = estimateRhoAreaMedian(removeHFCand ? selRemoveHFCand(inputParticles) : inputParticles);
  } else if (bkgSubEst == BkgSubEstimator::perpCone) {

    if (jets.size() == 0) {
      return sub;
    }
    std::tie(rhoParam, rhoMParam) = estimateRhoPerpCone(removeHFCand ? selRemoveHFCand(inputParticles) : inputParticles, jets);
  } else {
    rhoParam = 0.;
    rhoMParam = 0.;
    if (bkgSubEst != BkgSubEstimator::none) {
      LOGF(error, "requested estimator not implemented!");
    }
  }

  if (bkgSubMode == BkgSubMode::eventConstSub) { // eventwise subtraction

    fastjet::contrib::ConstituentSubtractor constituentSub(rhoParam, rhoMParam);
    constituentSub.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); /// deltaR=sqrt((y_i-y_j)^2+(phi_i-phi_j)^2)), longitudinal Lorentz invariant
    constituentSub.set_max_distance(constSubRMax);
    constituentSub.set_alpha(constSubAlpha);
    constituentSub.set_ghost_area(ghostAreaSpec.ghost_area());
    constituentSub.set_max_eta(bkgEtaMax);
    if (removeHFCand) {
      constituentSub.set_particle_selector(&selRemoveHFCand);
    }

    constituentSub.set_do_mass_subtraction();

    inputParticles = constituentSub.subtract_event(inputParticles);

    for (auto& track : inputParticles) {
      if (track.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
        std::cout << "HF candidate found after eventwise constituent subtraction\n";
      }
    }

  } else if (bkgSubMode == BkgSubMode::jetConstSub) { // jetwise subtraction

    if (jets.size() == 0) {
      return sub;
    }

    // FIXME There is a bug in fastjet ConstituentSubtractor at line 180, subtraction should be done on selected particles
    // The ConstituentSubtractor also needs to be modified such that it write the area information into the output jet
    fastjet::contrib::ConstituentSubtractor subtractor(rhoParam, rhoMParam);
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); /// deltaR=sqrt((y_i-y_j)^2+(phi_i-phi_j)^2)), longitudinal Lorentz invariant
    subtractor.set_max_distance(constSubRMax);
    subtractor.set_alpha(constSubAlpha);
    subtractor.set_ghost_area(ghostAreaSpec.ghost_area());
    subtractor.set_max_eta(bkgEtaMax);
    if (removeHFCand) {
      subtractor.set_particle_selector(&selRemoveHFCand);
    }

    // by default, the masses of all particles are set to zero. With this flag the jet mass will also be subtracted
    subtractor.set_do_mass_subtraction();

    // FIXME There are some problems with this method, it doesn't propagate the constituents user_info
    // Another problem is that, This method doesn't propagate the area information, since after constituent subtraction
    // the jet structure will change, so it no longer has the same area. fastjet developers said calculatig the area
    // information will difficult
    jets = subtractor(jets);

  } else if (bkgSubMode == BkgSubMode::rhoAreaSub) {
    sub = fastjet::Subtractor(rhoParam, rhoMParam);
  } else {
    rhoParam = 0.;
    rhoMParam = 0.;
    if (bkgSubMode != BkgSubMode::none) {
      LOGF(error, "requested subtraction mode not implemented!");
    }
  }

  return sub;
}

double JetBkgSubUtils::getM(fastjet::PseudoJet jet) const
{

  double sum(0);
  for (auto cons : jet.constituents()) {
    sum += TMath::Sqrt(cons.m() * cons.m() + cons.pt() * cons.pt()) - cons.pt(); // sqrt(E^2-P^2+pt^2)=sqrt(E^2-pz^2)
  }

  return sum;
}
