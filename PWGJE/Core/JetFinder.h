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
// Authors: Nima Zardoshti, Jochen Klein
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#ifndef PWGJE_CORE_JETFINDER_H_
#define PWGJE_CORE_JETFINDER_H_

#include <memory>
#include <vector>

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TMath.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/tools/Subtractor.hh"

enum class JetType {
  full = 0,
  charged = 1,
  neutral = 2,
};

class JetFinder
{

 public:
  /// Performs jet finding
  /// \note the input particle and jet lists are passed by reference
  /// \param inputParticles vector of input particles/tracks
  /// \param jets veector of jets to be filled
  /// \return ClusterSequenceArea object needed to access constituents
  // fastjet::ClusterSequenceArea findJets(std::vector<fastjet::PseudoJet> &inputParticles, std::vector<fastjet::PseudoJet> &jets);

  static constexpr float mPion = 0.139; // TDatabasePDG::Instance()->GetParticle(211)->Mass(); //can be removed when pion mass becomes default for unidentified tracks

  float phiMin = 0.;
  float phiMax = 2. * M_PI;
  float etaMin = -.9;
  float etaMax = .9;

  float jetR = .4;
  float jetPtMin = 0.;
  float jetPtMax = 1000.;
  float jetPhiMin = 0.;
  float jetPhiMax = 2. * M_PI;
  float jetEtaMin = -99.;
  float jetEtaMax = 99.;
  bool jetEtaDefault = false;

  float ghostEtaMin = -.9;
  float ghostEtaMax = .9;
  float ghostArea = .005;
  int ghostRepeatN = 1;
  double ghostktMean = 1.e-100;
  float gridScatter = 1.;
  float ktScatter = .1;

  bool isReclustering = false;
  bool isTriggering = false;

  fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::AreaType areaType = fastjet::active_area;
  fastjet::GhostedAreaSpec ghostAreaSpec;
  fastjet::JetDefinition jetDef;
  fastjet::AreaDefinition areaDef;
  fastjet::Selector selJets;
  fastjet::Selector selGhosts;

  /// Sets the jet finding parameters
  void setParams();

  /// Performs jet finding
  /// \note the input particle and jet lists are passed by reference
  /// \param inputParticles vector of input particles/tracks
  /// \param jets veector of jets to be filled
  /// \return ClusterSequenceArea object needed to access constituents
  fastjet::ClusterSequenceArea findJets(std::vector<fastjet::PseudoJet>& inputParticles, std::vector<fastjet::PseudoJet>& jets); // ideally find a way of passing the cluster sequence as a reeference

 private:
  ClassDefNV(JetFinder, 1);
};

#endif // PWGJE_CORE_JETFINDER_H_
