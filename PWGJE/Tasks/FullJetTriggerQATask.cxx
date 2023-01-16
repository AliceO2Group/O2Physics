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

// jet analysis tasks (subscribing to jet finder task)
//

// Author: Gijs van Weelden
//

#include "jettriggerqa.h"
/*
void JetTriggerQA::NEF(aod::Jet const& jet, std::vector<fastjet::PseudoJet> const& jetClusterConstituents)
{
  float nef = 0.;
  if (jetClusterConstituents.size() > 0) {
    for (const auto& constituent : jetClusterConstituents) {
      nef += constituent.E();
    }
    nef /= jet.energy();
  }
  spectra.fill(HIST("hJetNEF"), jet.pt(), nef);
}

void JetTriggerQA::SoftDrop(aod::Jet const& jet,
                            std::vector<fastjet::PseudoJet> const& jetReclustered,
                            float zg,
                            float rg,
                            int nsd)
{
  bool softDropped = false;
  fastjet::PseudoJet daughterSubJet = jetReclustered[0];
  fastjet::PseudoJet parentSubJet1;
  fastjet::PseudoJet parentSubJet2;
  while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
    if (parentSubJet1.perp() < parentSubJet2.perp()) {
      std::swap(parentSubJet1, parentSubJet2);
    }
    auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
    auto r = parentSubJet1.delta_R(parentSubJet2);
    if (z >= f_SD_zCut * TMath::Power(r / f_jetR, f_SD_beta)) {
      if (!softDropped) {
        zg = z;
        rg = r;
        spectra.fill(HIST("hJetZg"), jet.pt(), zg);
        spectra.fill(HIST("hJetRg"), jet.pt(), rg);
        softDropped = true;
      }
      nsd++;
    }
    daughterSubJet = parentSubJet1;
  }
  spectra.fill(HIST("hJetnsd"), jet.pt(), nsd);
}

void JetTriggerQA::ChargeFragmentation(aod::Jet const& jet, std::vector<fastjet::PseudoJet> const& jetConstituents)
{
  float z = 0.;
  for (auto track : jetConstituents){
    z = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
    z /= (jet.p() * jet.p());
    spectra.fill(HIST("hJetChargeFrag"), jet.pt(), z);
  }
}

float JetTriggerQA::Angularity(aod::Jet const& jet,
                               std::vector<fastjet::PseudoJet> const& jetConstituents,
                               std::vector<fastjet::PseudoJet> const& jetClusterConstituents,
                               float alpha,
                               float kappa)
{
  if (alpha == 0 || kappa == 0){
    std::cout << "Warning: requested angularity for zero alpha (" << alpha << ") or kappa (" << kappa << ") values."
              << std::endl
              << "Raising a number to the zeroth power will not give 1! No fix implemented yet."
              << std::endl;
  }
  if (alpha < 0 || kappa < 0){
    std::cout << "Error: requested angularity for negative alpha (" << alpha << ") or kappa (" << kappa << ") values." << std::endl;
    return -999.;
  }
  float deltaR = 0.;
  float ang = 0.;
  for (auto const& constituent : jetConstituents){
    deltaR = TMath::Sqrt(TMath::Power(jet.phi() - constituent.phi(), 2) + TMath::Power(jet.eta() - constituent.eta(), 2));
    ang += TMath::Power(constituent.pt()/jet.pt(), kappa) * TMath::Power(deltaR/f_jetR, alpha);
  }
  for (auto const& constituent : jetClusterConstituents){
    deltaR = TMath::Sqrt(TMath::Power(jet.phi() - constituent.phi(), 2) + TMath::Power(jet.eta() - constituent.eta(), 2));
    ang += TMath::Power(constituent.pt()/jet.pt(), kappa) * TMath::Power(deltaR/f_jetR, alpha);
  }
  return ang;
}

void JetTriggerQA::ptD(aod::Jet const& jet,
                       std::vector<fastjet::PseudoJet> const& jetConstituents,
                       std::vector<fastjet::PseudoJet> const& jetClusterConstituents)
{
  float ptD = 0.;
  for (auto const& constituent : jetConstituents){
    ptD += TMath::Power(constituent.pt(), 2);
  }
  for (auto const& constituent : jetClusterConstituents) {
    ptD += TMath::Power(constituent.pt(), 2);
  }
  ptD = TMath::Sqrt(ptD) / jet.pt();
  spectra.fill(HIST("hJetPtD"), jet.pt(), ptD);
}
*/
