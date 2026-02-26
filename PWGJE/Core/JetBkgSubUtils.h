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

#include <fastjet/AreaDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

#include <tuple>
#include <vector>

#include <math.h>

enum class BkgSubEstimator { none = 0,
                             medianRho = 1,
                             medianRhoSparse = 2
                             // perpendicular cone method is in JetUtilities
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

  JetBkgSubUtils(float jetBkgR_out, float bkgEtaMin_out = -0.9, float bkgEtaMax_out = 0.9,
                 float bkgPhiMin_out = 0., float bkgPhiMax_out = 2 * M_PI, float constSubAlpha_out = 1., float constSubRMax_out = 0.6, int nHardReject_out = 2, fastjet::GhostedAreaSpec ghostAreaSpec_out = fastjet::GhostedAreaSpec());

  // Default destructor
  ~JetBkgSubUtils() = default;

  /// @brief Setting the selectors after the input values have been initialised
  void initialise();

  /// @brief Setting the jet algorithm and the recombination scheme
  void setJetAlgorithmAndScheme(fastjet::JetAlgorithm algorithmBkg_out, fastjet::RecombinationScheme recombSchemeBkg_out)
  {
    algorithmBkg = algorithmBkg_out;
    recombSchemeBkg = recombSchemeBkg_out;
  }

  /// @brief Method for estimating the jet background density using the median method or the sparse method
  /// @param inputParticles (all particles in the event)
  /// @param doSparseSub weather to do rho sparse subtraction
  /// @return Rho, RhoM the underlying event density
  std::tuple<double, double> estimateRhoAreaMedian(const std::vector<fastjet::PseudoJet>& inputParticles, bool doSparseSub);

  /// @brief method that subtracts the background from jets using the area method
  /// @param jet input jet to be background subtracted
  /// @param rhoParam the underlying evvent density vs pT (to be set)
  /// @param rhoParam the underlying evvent density vs jet mass (to be set)
  /// @return jet, background subtracted jet
  fastjet::PseudoJet doRhoAreaSub(const fastjet::PseudoJet& jet, double rhoParam, double rhoMParam);

  /// @brief method that subtracts the background from the input particles using the event-wise cosntituent subtractor
  /// @param inputParticles (all the tracks/clusters/particles in the event)
  /// @param rhoParam the underlying evvent density vs pT (to be set)
  /// @param rhoParam the underlying evvent density vs jet mass (to be set)
  /// @return inputParticles, a vector of background subtracted input particles
  std::vector<fastjet::PseudoJet> doEventConstSub(std::vector<fastjet::PseudoJet>& inputParticles, double rhoParam, double rhoMParam);

  /// @brief method that subtracts the background from jets using the jet-wise constituent subtractor
  /// @param jets (all jets in the event)
  /// @param rhoParam the underlying evvent density vs pT (to be set)
  /// @param rhoParam the underlying evvent density vs jet mass (to be set)
  /// @return jets, a vector of background subtracted jets
  std::vector<fastjet::PseudoJet> doJetConstSub(std::vector<fastjet::PseudoJet>& jets, double rhoParam, double rhoMParam);

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
  void setDoRhoMassSub(bool doMSub_out = true) { doRhoMassSub = doMSub_out; }
  void setGhostAreaSpec(fastjet::GhostedAreaSpec ghostAreaSpec_out) { ghostAreaSpec = ghostAreaSpec_out; }

  // Getters
  float getJetBkgR() const { return jetBkgR; }
  float getPhiMin() const { return bkgPhiMin; }
  float getPhiMax() const { return bkgPhiMax; }
  float getEtaMin() const { return bkgEtaMin; }
  float getEtaMax() const { return bkgEtaMax; }
  float getConstSubAlpha() const { return constSubAlpha; }
  float getConstSubRMax() const { return constSubRMax; }
  float getDoRhoMassSub() const { return doRhoMassSub; }
  fastjet::GhostedAreaSpec getGhostAreaSpec() const { return ghostAreaSpec; }
  fastjet::JetDefinition getJetDefinition() const { return jetDefBkg; }
  fastjet::AreaDefinition getAreaDefinition() const { return areaDefBkg; }
  fastjet::Selector getRhoSelector() const { return selRho; }

  // Calculate the jet mass
  double getMd(fastjet::PseudoJet jet) const;

 protected:
  float jetBkgR = 0.2;
  float bkgEtaMin = -0.9;
  float bkgEtaMax = 0.9;
  float bkgPhiMin = 0.0;
  float bkgPhiMax = 2.0 * M_PI;
  float constSubAlpha = 1.0;
  float constSubRMax = 0.24;
  int nHardReject = 2;
  bool doRhoMassSub = false; /// flag whether to do jet mass subtraction with the const sub

  fastjet::GhostedAreaSpec ghostAreaSpec = fastjet::GhostedAreaSpec();
  fastjet::JetAlgorithm algorithmBkg = fastjet::kt_algorithm;
  fastjet::RecombinationScheme recombSchemeBkg = fastjet::E_scheme;
  fastjet::JetDefinition jetDefBkg = fastjet::JetDefinition(algorithmBkg, jetBkgR, recombSchemeBkg, fastjet::Best);
  fastjet::AreaDefinition areaDefBkg = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghostAreaSpec);
  fastjet::Selector selRho = fastjet::Selector();

}; // class JetBkgSubUtils

#endif // PWGJE_CORE_JETBKGSUBUTILS_H_
