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
// Author: Hadi Hassan
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "Framework/Logger.h"
#include "TVector2.h"

JetBkgSubUtils::JetBkgSubUtils(float jetBkgR, float bkgPhiMin, float bkgPhiMax, float bkgEtaMin, float bkgEtaMax, float constSubAlpha, float constSubRMax, fastjet::GhostedAreaSpec ghostAreaSpec) : mJetBkgR(jetBkgR),
                                                                                                                                                                                                     mBkgPhiMin(bkgPhiMin),
                                                                                                                                                                                                     mBkgPhiMax(bkgPhiMax),
                                                                                                                                                                                                     mBkgEtaMin(bkgEtaMin),
                                                                                                                                                                                                     mBkgEtaMax(bkgEtaMax),
                                                                                                                                                                                                     mConstSubAlpha(constSubAlpha),
                                                                                                                                                                                                     mConstSubRMax(constSubRMax),
                                                                                                                                                                                                     mGhostAreaSpec(ghostAreaSpec)

{
  mJetDefBkg = fastjet::JetDefinition(mAlgorithmBkg, mJetBkgR, mRecombSchemeBkg, fastjet::Best);
  mAreaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, mGhostAreaSpec);
  mSelRho = fastjet::SelectorRapRange(mBkgEtaMin, mBkgEtaMax) && fastjet::SelectorPhiRange(mBkgPhiMin, mBkgPhiMax) && !fastjet::SelectorNHardest(2); // here we have to put rap range, to be checked!
}

double JetBkgSubUtils::estimateRhoAreaMedian(std::vector<fastjet::PseudoJet>& inputParticles)
{

  if (inputParticles.size() == 0) {
    return 0.0;
  }

  // cluster the kT jets
  fastjet::ClusterSequenceArea clusterSeq(inputParticles, mJetDefBkg, mAreaDefBkg);

  // select jets in detector acceptance
  std::vector<fastjet::PseudoJet> alljets = mSelRho(clusterSeq.inclusive_jets());

  double totaljetAreaPhys(0), totalAreaCovered(0);
  std::vector<double> rhovector;

  // Fill a vector for pT/area to be used for the median
  for (auto& ijet : alljets) {

    // Physical area/ Physical jets (no ghost)
    if (!clusterSeq.is_pure_ghost(ijet)) {
      rhovector.push_back(ijet.perp() / ijet.area());
      totaljetAreaPhys += ijet.area();
    }
    // Full area
    totalAreaCovered += ijet.area();
  }
  // calculate Rho as the median of the jet pT / jet area
  double rho = TMath::Median<double>(rhovector.size(), rhovector.data());

  if (mDoSparseSub) {
    // calculate The ocupancy factor, which the ratio of covered area / total area
    double occupancyFactor = totalAreaCovered > 0 ? totaljetAreaPhys / totalAreaCovered : 1.;
    rho *= occupancyFactor;
  }

  return rho;
}

double JetBkgSubUtils::estimateRhoPerpCone(const std::vector<fastjet::PseudoJet>& inputParticles, const std::vector<fastjet::PseudoJet>& jets)
{

  if (inputParticles.size() == 0 || jets.size() == 0) {
    return 0.0;
  }

  double perpPtDensity1 = 0;
  double perpPtDensity2 = 0;

  fastjet::Selector selectJet = fastjet::SelectorEtaRange(mBkgEtaMin, mBkgEtaMax) && fastjet::SelectorPhiRange(mBkgPhiMin, mBkgPhiMax);

  std::vector<fastjet::PseudoJet> selectedJets = fastjet::sorted_by_pt(selectJet(jets));

  if (selectedJets.size() == 0) {
    return 0.0;
  }

  fastjet::PseudoJet leadingJet = selectedJets[0];

  double dPhi1 = 999.;
  double dPhi2 = 999.;
  double dEta = 999.;
  double PerpendicularConeAxisPhi1 = 999., PerpendicularConeAxisPhi2 = 999.;
  // build 2 perp cones in phi around the leading jet (right and left of the jet)
  PerpendicularConeAxisPhi1 = TVector2::Phi_0_2pi(leadingJet.phi() + (M_PI / 2.));
  PerpendicularConeAxisPhi2 = TVector2::Phi_0_2pi(leadingJet.phi() - (M_PI / 2.));

  for (auto& particle : inputParticles) {
    // sum the momentum of all paricles that fill the two cones
    dPhi1 = particle.phi() - PerpendicularConeAxisPhi1;
    dPhi1 = TVector2::Phi_mpi_pi(dPhi1);
    dPhi2 = particle.phi() - PerpendicularConeAxisPhi2;
    dPhi2 = TVector2::Phi_mpi_pi(dPhi2);
    dEta = leadingJet.eta() - particle.eta();
    if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= mJetBkgR) {
      perpPtDensity1 += particle.perp();
    }

    if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= mJetBkgR) {
      perpPtDensity2 += particle.perp();
    }
  }

  // Caculate rho as the ratio of average pT of the two cones / the cone area
  Double_t perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * M_PI * mJetBkgR * mJetBkgR);

  return perpPtDensity;
}

/// Sets the background subtraction estimater pointer
void JetBkgSubUtils::setBkgE()
{
  mbkgE = decltype(mbkgE)(new fastjet::JetMedianBackgroundEstimator(mSelRho, mJetDefBkg, mAreaDefBkg));
}

/// Sets the background subtraction object
fastjet::Subtractor JetBkgSubUtils::setSub(std::vector<fastjet::PseudoJet>& inputParticles, float& rhoParam, BkgSubMode bkgSubMode, std::vector<fastjet::PseudoJet>& jets)
{

  if (bkgSubMode == BkgSubMode::rhoSparseSub) {
    mDoSparseSub = true;
  }

  fastjet::Subtractor sub;

  if (bkgSubMode == BkgSubMode::rhoMedianAreaSub || bkgSubMode == BkgSubMode::rhoSparseSub) {

    rhoParam = estimateRhoAreaMedian(inputParticles);
    sub = fastjet::Subtractor(rhoParam);

  } else if (bkgSubMode == BkgSubMode::rhoPerpConeSub) {

    if (jets.size() == 0) {
      return sub;
    }
    rhoParam = estimateRhoPerpCone(inputParticles, jets);
    sub = fastjet::Subtractor(rhoParam);

  } else if (bkgSubMode == BkgSubMode::rhoAreaSub) {

    setBkgE();
    mbkgE->set_particles(inputParticles);
    sub = fastjet::Subtractor(mbkgE.get());
    rhoParam = mbkgE->estimate().rho();

  } else if (bkgSubMode == BkgSubMode::constSub) { // event or jetwise

    setBkgE();
    mbkgE->set_particles(inputParticles);
    rhoParam = mbkgE->estimate().rho();
    std::unique_ptr<fastjet::contrib::ConstituentSubtractor> constituentSub;
    constituentSub = decltype(constituentSub){new fastjet::contrib::ConstituentSubtractor{mbkgE.get()}};
    constituentSub->set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    constituentSub->set_max_distance(mConstSubRMax);
    constituentSub->set_alpha(mConstSubAlpha);
    constituentSub->set_ghost_area(mGhostAreaSpec.ghost_area());
    constituentSub->set_max_eta(mBkgEtaMax);
    constituentSub->set_background_estimator(mbkgE.get());

    inputParticles = constituentSub->subtract_event(inputParticles);

  } else if (bkgSubMode == BkgSubMode::jetconstSub) {

    if (jets.size() == 0) {
      return sub;
    }

    setBkgE();
    fastjet::BackgroundJetScalarPtDensity scalarPtDensity;
    mbkgE->set_jet_density_class(&scalarPtDensity); // this changes the computation of pt of patches from vector sum to scalar sum. The scalar sum seems more reasonable.
    mbkgE->set_particles(inputParticles);

    fastjet::contrib::ConstituentSubtractor subtractor(mbkgE.get());
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    subtractor.set_max_distance(mConstSubRMax);
    subtractor.set_alpha(mConstSubAlpha);
    subtractor.set_ghost_area(mGhostAreaSpec.ghost_area());
    subtractor.set_max_eta(mBkgEtaMax);

    // by default, the masses of all particles are set to zero.
    // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
    subtractor.set_do_mass_subtraction();
    subtractor.set_common_bge_for_rho_and_rhom();

    // Unfortunately this method doesn't work, it crashes with error:
    // fastjet::Error:  One or more of this composite jet's pieces does not support area
    jets = subtractor(jets);

  } else {
    rhoParam = 0.;
    if (bkgSubMode != BkgSubMode::none) {
      LOGF(error, "requested subtraction mode not implemented!");
    }
  }

  return sub;
}
