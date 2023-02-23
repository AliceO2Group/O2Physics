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

#ifndef PWGJE_CORE_JETBKGSUBUTILS_H
#define PWGJE_CORE_JETBKGSUBUTILS_H

#include <memory>
#include <vector>
#include <optional>
#include <TMath.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "Framework/Logger.h"

enum class BkgSubMode { none = 0,
                        rhoAreaSub = 1,
                        constSub = 2,
                        rhoSparseSub = 3,
                        rhoPerpConeSub = 4,
                        rhoMedianAreaSub = 5,
                        jetconstSub = 6
};

class JetBkgSubUtils
{

 public:
  // Default contructor
  JetBkgSubUtils() = default;

  JetBkgSubUtils(float jetBkgR_out = 0.2, float bkgPhiMin_out = 0., float bkgPhiMax_out = 2 * M_PI, float bkgEtaMin_out = -0.9, float bkgEtaMax_out = 0.9,
                 float constSubAlpha_out = 1., float constSubRMax_out = 0.6, fastjet::GhostedAreaSpec ghostAreaSpec_out = fastjet::GhostedAreaSpec());

  // Default destructor
  ~JetBkgSubUtils() = default;

  /// @brief Setting the jet algorithm and the recombination scheme
  void setJetAlgorithmAndScheme(fastjet::JetAlgorithm algorithmBkg, fastjet::RecombinationScheme recombSchemeBkg)
  {
    mAlgorithmBkg = algorithmBkg;
    mRecombSchemeBkg = recombSchemeBkg;
  }

  /// @brief Method for estimating the jet background density using the median method or the sparse method
  /// @param inputParticles (all particles in the event)
  /// @return Rho, the underlying event density
  double estimateRhoAreaMedian(std::vector<fastjet::PseudoJet>& inputParticles);

  /// @brief Background estimator using the perpendicular cone method
  /// @param inputParticles
  /// @param jets (all jets in the event)
  /// @return Rho, the underlying event density
  double estimateRhoPerpCone(const std::vector<fastjet::PseudoJet>& inputParticles, const std::vector<fastjet::PseudoJet>& jets);

  /// @brief method that reconstruct the bkgsub object and estimates rho
  /// @param inputParticles
  /// @param jets (all jets in the event)
  /// @param rhoParam the underlying evvent density (to be set)
  /// @param bkgSubMode the background subtraction method
  /// @param jets signal jets (needed for the perp cone method)
  /// @return sub, the background subtractor
  fastjet::Subtractor setSub(std::vector<fastjet::PseudoJet>& inputParticles, float& rhoParam, BkgSubMode bkgSubMode, std::vector<fastjet::PseudoJet>& jets);

  /// Sets the background subtraction estimater pointer
  void setBkgE();

  // Setters
  void setJetBkgR(float jetbkgR) { mJetBkgR = jetbkgR; }
  void setPhiMinMax(float phimin, float phimax)
  {
    mBkgPhiMin = phimin;
    mBkgPhiMax = phimax;
  }
  void setEtaMinMax(float etamin, float etamax)
  {
    mBkgEtaMin = etamin;
    mBkgEtaMax = etamax;
  }
  void setConstSubAlphaRMax(float alpha, float rmax)
  {
    mConstSubAlpha = alpha;
    mConstSubRMax = rmax;
  }
  void setDoRhoSparseSub(bool dosparse = true) { mDoSparseSub = dosparse; }
  void setGhostAreaSpec(fastjet::GhostedAreaSpec ghostAreaSpec) { mGhostAreaSpec = ghostAreaSpec; }
  void setJetDefinition(fastjet::JetDefinition jetdefbkg) { mJetDefBkg = jetdefbkg; }
  void setAreaDefinition(fastjet::AreaDefinition areaDefBkg) { mAreaDefBkg = areaDefBkg; }
  void setRhoSelector(fastjet::Selector selRho) { mSelRho = selRho; }

  // Getters
  float getJetBkgR() const { return mJetBkgR; }
  float getPhiMin() const { return mBkgPhiMin; }
  float getPhiMax() const { return mBkgPhiMax; }
  float getEtaMin() const { return mBkgEtaMin; }
  float getEtaMax() const { return mBkgEtaMax; }
  float getConstSubAlpha() const { return mConstSubAlpha; }
  float getConstSubRMax() const { return mConstSubRMax; }
  float getDoRhoSparseSub() const { return mDoSparseSub; }
  fastjet::GhostedAreaSpec getGhostAreaSpec() const { return mGhostAreaSpec; }
  fastjet::JetDefinition getJetDefinition() const { return mJetDefBkg; }
  fastjet::AreaDefinition getAreaDefinition() const { return mAreaDefBkg; }
  fastjet::Selector getRhoSelector() const { return mSelRho; }

 protected:
  float mJetBkgR = 0.2;
  float mBkgPhiMin = 0.;
  float mBkgPhiMax = 2 * M_PI;
  float mBkgEtaMin = -0.9;
  float mBkgEtaMax = 0.9;
  float mConstSubAlpha = 1.0;
  float mConstSubRMax = 0.6;
  bool mDoSparseSub = false; /// flag whether to do background subtraction for sparse systems.

  std::unique_ptr<fastjet::JetMedianBackgroundEstimator> mbkgE;

  fastjet::GhostedAreaSpec mGhostAreaSpec = fastjet::GhostedAreaSpec();
  fastjet::JetAlgorithm mAlgorithmBkg = fastjet::kt_algorithm;
  fastjet::RecombinationScheme mRecombSchemeBkg = fastjet::E_scheme;
  fastjet::JetDefinition mJetDefBkg = fastjet::JetDefinition(mAlgorithmBkg, mJetBkgR, mRecombSchemeBkg, fastjet::Best);
  fastjet::AreaDefinition mAreaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, mGhostAreaSpec);
  fastjet::Selector mSelRho = fastjet::Selector();

}; // class JetBkgSubUtils

#endif // PWGJE_CORE_JETBKGSUBUTILS_H
