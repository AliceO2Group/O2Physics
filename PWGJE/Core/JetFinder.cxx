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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "PWGJE/Core/JetFinder.h"

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

#include <vector>

/// Sets the jet finding parameters
void JetFinder::setParams()
{
  if (jetEtaDefault) {

    jetEtaMin = etaMin + jetR; // in aliphysics this was (-etaMax + 0.95*jetR)
    jetEtaMax = etaMax - jetR;

    if (isReclustering || isTriggering) {
      jetEtaMin -= jetR;
      jetEtaMax += jetR;
    }
  }
  float jetRForClustering = isReclustering ? 5.0 * jetR : jetR;

  selGhosts = fastjet::SelectorEtaRange(ghostEtaMin, ghostEtaMax) && fastjet::SelectorPhiRange(phiMin, phiMax);
  // ghostAreaSpec=fastjet::GhostedAreaSpec(selGhosts,ghostRepeatN,ghostArea,gridScatter,ktScatter,ghostktMean);
  ghostAreaSpec = fastjet::GhostedAreaSpec(selGhosts, ghostRepeatN, ghostArea, gridScatter, ktScatter, ghostktMean); // the first argument is rapidity not pseudorapidity, to be checked
  jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, jetRForClustering, recombScheme, strategy);
  if (fastjetExtraParam > -98.0) {
    jetDef.set_extra_param(fastjetExtraParam);
  }
  jetDef.set_jet_algorithm(algorithm);
  areaDef = fastjet::AreaDefinition(areaType, ghostAreaSpec);
  selJets = fastjet::SelectorPtRange(jetPtMin, jetPtMax) && fastjet::SelectorEtaRange(jetEtaMin, jetEtaMax) && fastjet::SelectorPhiRange(jetPhiMin, jetPhiMax);
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
  fastjet::ClusterSequenceArea clusterSeq(inputParticles, jetDef, areaDef);
  jets = clusterSeq.inclusive_jets();
  jets = selJets(jets);
  jets = fastjet::sorted_by_pt(jets);
  return clusterSeq;
}
