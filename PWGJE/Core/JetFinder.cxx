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
// Author: Jochen Klein, Nima Zardoshti
#include "PWGJE/Core/JetFinder.h"
#include "Framework/Logger.h"
#include "PWGJE/Core/JetBkgSubUtils.h"

/// Sets the jet finding parameters
void JetFinder::setParams()
{

  if (!isReclustering) {
    jetEtaMin = etaMin + jetR; // in aliphysics this was (-etaMax + 0.95*jetR)
    jetEtaMax = etaMax - jetR;
  } else {
    jetR = 5.0 * jetR;
  }

  // selGhosts =fastjet::SelectorRapRange(ghostEtaMin,ghostEtaMax) && fastjet::SelectorPhiRange(phiMin,phiMax);
  // ghostAreaSpec=fastjet::GhostedAreaSpec(selGhosts,ghostRepeatN,ghostArea,gridScatter,ktScatter,ghostktMean);
  ghostAreaSpec = fastjet::GhostedAreaSpec(ghostEtaMax, ghostRepeatN, ghostArea, gridScatter, ktScatter, ghostktMean); // the first argument is rapidity not pseudorapidity, to be checked
  jetDef = fastjet::JetDefinition(algorithm, jetR, recombScheme, strategy);
  areaDef = fastjet::AreaDefinition(areaType, ghostAreaSpec);
  selJets = fastjet::SelectorPtRange(jetPtMin, jetPtMax) && fastjet::SelectorEtaRange(jetEtaMin, jetEtaMax) && fastjet::SelectorPhiRange(jetPhiMin, jetPhiMax);

  mSubUtils = std::make_unique<JetBkgSubUtils>(jetR, bkgPhiMin, bkgPhiMax, bkgEtaMin, bkgEtaMax, constSubAlpha, constSubRMax, ghostAreaSpec);
  mSubUtils->setJetAlgorithmAndScheme(algorithmBkg, recombSchemeBkg);
}

/// Performs jet finding
/// \note the input particle and jet lists are passed by reference
/// \param inputParticles vector of input particles/tracks
/// \param jets veector of jets to be filled
/// \return ClusterSequenceArea object needed to access constituents
fastjet::ClusterSequenceArea JetFinder::findJets(std::vector<fastjet::PseudoJet>& inputParticles, std::vector<fastjet::PseudoJet>& jets) // ideally find a way of passing the cluster sequence as a reeference
{
  setParams();
  jets.clear();

  sub = mSubUtils->setSub(inputParticles, mRho, bkgSubMode);

  fastjet::ClusterSequenceArea clusterSeq(inputParticles, jetDef, areaDef);

  if (bkgSubMode == BkgSubMode::rhoPerpConeSub) {
    sub = mSubUtils->setSub(inputParticles, mRho, bkgSubMode, clusterSeq.inclusive_jets());
  }

  jets = (mRho > DBL_EPSILON && bkgSubMode != BkgSubMode::constSub) ? (sub)(clusterSeq.inclusive_jets()) : clusterSeq.inclusive_jets();
  jets = selJets(jets);
  if (isReclustering) {
    jetR = jetR / 5.0;
  }
  return clusterSeq;
}
