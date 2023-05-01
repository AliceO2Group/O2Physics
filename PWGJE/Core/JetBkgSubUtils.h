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

/// \file JetBkgSubUtils.h
/// \brief Jet background subtraction utilities
///
/// \author Hadi Hassan <hadi.hassan@cern.ch>, JYU

#ifndef PWGJE_CORE_JETBKGSUBUTILS_H_
#define PWGJE_CORE_JETBKGSUBUTILS_H_

#include <string>
#include <memory>
#include <tuple>
#include <vector>
#include <TMath.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "Framework/Logger.h"

enum class BkgSubEstimator { none = 0,
                             medianRho = 1,
                             medianRhoSparse = 2,
                             perpCone = 3
};

enum class BkgSubMode { none = 0,
                        rhoAreaSub = 1,
                        eventConstSub = 2,
                        jetConstSub = 3
};

class JetBkgSubUtils
{

 public:
  // Default contructor
  JetBkgSubUtils() = default;

  JetBkgSubUtils(float jetBkgR_out = 0.2, float constSubAlpha_out = 1., float constSubRMax_out = 0.6, float bkgEtaMin_out = -0.8, float bkgEtaMax_out = 0.8,
                 float bkgPhiMin_out = 0., float bkgPhiMax_out = 2 * M_PI, fastjet::GhostedAreaSpec ghostAreaSpec_out = fastjet::GhostedAreaSpec(), int nHardReject = 2);

  // Default destructor
  ~JetBkgSubUtils() = default;

  /// @brief Setting the jet algorithm and the recombination scheme
  void setJetAlgorithmAndScheme(fastjet::JetAlgorithm algorithmBkg_out, fastjet::RecombinationScheme recombSchemeBkg_out)
  {
    algorithmBkg = algorithmBkg_out;
    recombSchemeBkg = recombSchemeBkg_out;
  }

  /// @brief Method for estimating the jet background density using the median method or the sparse method
  /// @param inputParticles (all particles in the event)
  /// @return Rho, RhoM the underlying event density
  std::tuple<double, double> estimateRhoAreaMedian(const std::vector<fastjet::PseudoJet>& inputParticles);

  /// @brief Background estimator using the perpendicular cone method
  /// @param inputParticles
  /// @param jets (all jets in the event)
  /// @return Rho, RhoM the underlying event density
  std::tuple<double, double> estimateRhoPerpCone(const std::vector<fastjet::PseudoJet>& inputParticles, const std::vector<fastjet::PseudoJet>& jets);

  /// @brief method that reconstruct the bkgsub object and estimates rho
  /// @param inputParticles
  /// @param jets (all jets in the event)
  /// @param rhoParam the underlying evvent density (to be set)
  /// @param bkgSubMode the background subtraction method
  /// @return sub, the background subtractor
  fastjet::Subtractor setSub(std::vector<fastjet::PseudoJet>& inputParticles, std::vector<fastjet::PseudoJet>& jets, float& rhoParam, float& rhoMParam, BkgSubEstimator bkgSubEst, BkgSubMode bkgSubMode);

  // Setters
  void setJetBkgR(float jetbkgR_out) { jetBkgR = jetbkgR_out; }
  void setPhiMinMax(float phimin_out, float phimax_out)
  {
    bkgPhiMin = phimin_out;
    bkgPhiMax = phimax_out;
  }
  void setEtaMinMax(float etamin_out, float etamax_out)
  {
    bkgEtaMin = etamin_out;
    bkgEtaMax = etamax_out;
  }
  void setConstSubAlphaRMax(float alpha_out, float rmax_out)
  {
    constSubAlpha = alpha_out;
    constSubRMax = rmax_out;
  }
  void setDoRhoSparseSub(bool dosparse_out = true) { doSparseSub = dosparse_out; }
  void setRemoveHFCandidate(bool removecandidate = true) { removeHFCand = removecandidate; }
  void setGhostAreaSpec(fastjet::GhostedAreaSpec ghostAreaSpec_out) { ghostAreaSpec = ghostAreaSpec_out; }
  void setJetDefinition(fastjet::JetDefinition jetdefbkg_out) { jetDefBkg = jetdefbkg_out; }
  void setAreaDefinition(fastjet::AreaDefinition areaDefBkg_out) { areaDefBkg = areaDefBkg_out; }
  void setRhoSelector(fastjet::Selector selRho_out) { selRho = selRho_out; }

  // Getters
  float getJetBkgR() const { return jetBkgR; }
  float getPhiMin() const { return bkgPhiMin; }
  float getPhiMax() const { return bkgPhiMax; }
  float getEtaMin() const { return bkgEtaMin; }
  float getEtaMax() const { return bkgEtaMax; }
  float getConstSubAlpha() const { return constSubAlpha; }
  float getConstSubRMax() const { return constSubRMax; }
  float getDoRhoSparseSub() const { return doSparseSub; }
  float getRemoveHFCandidate() const { return removeHFCand; }
  fastjet::GhostedAreaSpec getGhostAreaSpec() const { return ghostAreaSpec; }
  fastjet::JetDefinition getJetDefinition() const { return jetDefBkg; }
  fastjet::AreaDefinition getAreaDefinition() const { return areaDefBkg; }
  fastjet::Selector getRhoSelector() const { return selRho; }

  // Calculate the jet mass
  double getM(fastjet::PseudoJet jet) const;

 protected:
  float jetBkgR = 0.2;
  float constSubAlpha = 1.0;
  float constSubRMax = 0.6;
  float bkgEtaMin = -0.8;
  float bkgEtaMax = 0.8;
  float bkgPhiMin = 0.;
  float bkgPhiMax = 2 * M_PI;
  bool doSparseSub = false;  /// flag whether to do background subtraction for sparse systems.
  bool removeHFCand = false; /// flag whether to remove the HF candidate from the list of particles

  fastjet::GhostedAreaSpec ghostAreaSpec = fastjet::GhostedAreaSpec();
  fastjet::JetAlgorithm algorithmBkg = fastjet::kt_algorithm;
  fastjet::RecombinationScheme recombSchemeBkg = fastjet::E_scheme;
  fastjet::JetDefinition jetDefBkg = fastjet::JetDefinition(algorithmBkg, jetBkgR, recombSchemeBkg, fastjet::Best);
  fastjet::AreaDefinition areaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghostAreaSpec);
  fastjet::Selector selRho = fastjet::Selector();

}; // class JetBkgSubUtils

#endif // PWGJE_CORE_JETBKGSUBUTILS_H_
