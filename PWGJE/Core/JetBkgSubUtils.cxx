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

JetBkgSubUtils::JetBkgSubUtils(float jetBkgR_out, float bkgEtaMin_out, float bkgEtaMax_out, float bkgPhiMin_out, float bkgPhiMax_out, float constSubAlpha_out, float constSubRMax_out, int nHardReject_out, fastjet::GhostedAreaSpec ghostAreaSpec_out) : jetBkgR(jetBkgR_out),
                                                                                                                                                                                                                                                          bkgEtaMin(bkgEtaMin_out),
                                                                                                                                                                                                                                                          bkgEtaMax(bkgEtaMax_out),
                                                                                                                                                                                                                                                          bkgPhiMin(bkgPhiMin_out),
                                                                                                                                                                                                                                                          bkgPhiMax(bkgPhiMax_out),
                                                                                                                                                                                                                                                          constSubAlpha(constSubAlpha_out),
                                                                                                                                                                                                                                                          constSubRMax(constSubRMax_out),
                                                                                                                                                                                                                                                          nHardReject(nHardReject_out),
                                                                                                                                                                                                                                                          ghostAreaSpec(ghostAreaSpec_out)

{
}

void JetBkgSubUtils::initialise()
{
  // Note: if you are using the PerpCone method you should jetBkgR to be the same as the anit_kt jets R, otherwise use R=0.2
  jetDefBkg = fastjet::JetDefinition(algorithmBkg, jetBkgR, recombSchemeBkg, fastjet::Best);
  areaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghostAreaSpec);
  selRho = fastjet::SelectorRapRange(bkgEtaMin, bkgEtaMax) && fastjet::SelectorPhiRange(bkgPhiMin, bkgPhiMax) && !fastjet::SelectorNHardest(nHardReject); // here we have to put rap range, to be checked!
}

std::tuple<double, double> JetBkgSubUtils::estimateRhoAreaMedian(const std::vector<fastjet::PseudoJet>& inputParticles, bool doSparseSub)
{
  JetBkgSubUtils::initialise();

  if (inputParticles.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  // cluster the kT jets
  fastjet::ClusterSequenceArea clusterSeq(inputParticles, jetDefBkg, areaDefBkg);

  // select jets in detector acceptance
  std::vector<fastjet::PseudoJet> alljets = selRho(clusterSeq.inclusive_jets());

  double totaljetAreaPhys(0), totalAreaCovered(0);
  std::vector<double> rhovector;
  std::vector<double> rhoMdvector;

  // Fill a vector for pT/area to be used for the median
  for (auto& ijet : alljets) {

    // Physical area/ Physical jets (no ghost)
    if (!clusterSeq.is_pure_ghost(ijet)) {
      rhovector.push_back(ijet.perp() / ijet.area());
      rhoMdvector.push_back(getMd(ijet) / ijet.area());

      totaljetAreaPhys += ijet.area();
    }
    // Full area
    totalAreaCovered += ijet.area();
  }
  // calculate Rho as the median of the jet pT / jet area

  double rho = 0.0;
  double rhoM = 0.0;
  if (rhovector.size() != 0) {
    rho = TMath::Median<double>(rhovector.size(), rhovector.data());
    rhoM = TMath::Median<double>(rhoMdvector.size(), rhoMdvector.data());
  }

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

  JetBkgSubUtils::initialise();
  if (inputParticles.size() == 0 || jets.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  double perpPtDensity1 = 0;
  double perpPtDensity2 = 0;
  double perpMdDensity1 = 0;
  double perpMdDensity2 = 0;

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
    dPhi1 = RecoDecay::constrainAngle<double, double>(dPhi1, -M_PI); // This will contrain the angel between -pi & Pi
    dPhi2 = particle.phi() - PerpendicularConeAxisPhi2;
    dPhi2 = RecoDecay::constrainAngle<double, double>(dPhi2, -M_PI); // This will contrain the angel between -pi & Pi
    dEta = leadingJet.eta() - particle.eta();                        // The perp cone eta is the same as the leading jet since the cones are perpendicular only in phi
    if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= jetBkgR) {
      perpPtDensity1 += particle.perp();
      perpMdDensity1 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }

    if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= jetBkgR) {
      perpPtDensity2 += particle.perp();
      perpMdDensity2 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }
  }

  // Caculate rho as the ratio of average pT of the two cones / the cone area
  double perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * M_PI * jetBkgR * jetBkgR);
  double perpMdDensity = (perpMdDensity1 + perpMdDensity2) / (2 * M_PI * jetBkgR * jetBkgR);

  return std::make_tuple(perpPtDensity, perpMdDensity);
}

fastjet::PseudoJet JetBkgSubUtils::doRhoAreaSub(fastjet::PseudoJet& jet, double rhoParam, double rhoMParam)
{

  fastjet::Subtractor sub = fastjet::Subtractor(rhoParam, rhoMParam);
  if (doRhoMassSub) {
    sub.set_safe_mass();
  }
  return sub(jet);
}

std::vector<fastjet::PseudoJet> JetBkgSubUtils::doEventConstSub(std::vector<fastjet::PseudoJet>& inputParticles, double rhoParam, double rhoMParam)
{
  JetBkgSubUtils::initialise();
  fastjet::contrib::ConstituentSubtractor constituentSub(rhoParam, rhoMParam);
  constituentSub.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); /// deltaR=sqrt((y_i-y_j)^2+(phi_i-phi_j)^2)), longitudinal Lorentz invariant
  constituentSub.set_max_distance(constSubRMax);
  constituentSub.set_alpha(constSubAlpha);
  constituentSub.set_ghost_area(ghostAreaSpec.ghost_area());
  constituentSub.set_max_eta(maxEtaEvent);

  // by default, the masses of all particles are set to zero. With this flag the jet mass will also be subtracted
  if (doRhoMassSub) {
    constituentSub.set_do_mass_subtraction();
  }

  return constituentSub.subtract_event(inputParticles, maxEtaEvent);
}

std::vector<fastjet::PseudoJet> JetBkgSubUtils::doJetConstSub(std::vector<fastjet::PseudoJet>& jets, double rhoParam, double rhoMParam)
{
  JetBkgSubUtils::initialise();
  if (jets.size() == 0) {
    return std::vector<fastjet::PseudoJet>();
  }

  // FIXME, this method works only if the input jets "jets" are reconstructed with area def "active_area_explicit_ghosts"
  //  because it needs the ghosts to estimate the backgeound
  fastjet::contrib::ConstituentSubtractor constituentSub(rhoParam, rhoMParam);
  constituentSub.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); /// deltaR=sqrt((y_i-y_j)^2+(phi_i-phi_j)^2)), longitudinal Lorentz invariant
  constituentSub.set_max_distance(constSubRMax);
  constituentSub.set_alpha(constSubAlpha);
  constituentSub.set_ghost_area(ghostAreaSpec.ghost_area());
  constituentSub.set_max_eta(bkgEtaMax);

  // by default, the masses of all particles are set to zero. With this flag the jet mass will also be subtracted
  if (doRhoMassSub) {
    constituentSub.set_do_mass_subtraction();
  }

  // FIXME, This method doesn't propagate the area information, since after constituent subtraction
  // the jet structure will change, so it no longer has the same area. fastjet developers said calculatig the area
  // information will difficult
  return constituentSub(jets);
}

double JetBkgSubUtils::getMd(fastjet::PseudoJet jet) const
{
  // Refere to https://arxiv.org/abs/1211.2811 for the rhoM caclulation
  double sum(0);
  for (auto constituent : jet.constituents()) {
    sum += TMath::Sqrt(constituent.m() * constituent.m() + constituent.pt() * constituent.pt()) - constituent.pt();
  }

  return sum;
}
