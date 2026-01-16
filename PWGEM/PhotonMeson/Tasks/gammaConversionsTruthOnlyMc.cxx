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

/// \file gammaConversionsTruthOnlyMc.cxx
/// \brief extract relevant mc truth information that allows to compute efficiency, purity and more quantities for the photon conversion analysis.
/// \author stephan.friedrich.stiefelmaier@cern.ch
/// dependencies: none
// SFS todo: update description

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <TVector3.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct gammaConversionsTruthOnlyMc {

  Configurable<bool> fPhysicalPrimaryOnly{"fPhysicalPrimaryOnly", true, "fPhysicalPrimaryOnly"};
  Configurable<float> fEtaMax{"fEtaMax", 0.8, "aMaximum photon eta"};
  Configurable<float> fV0RMin{"fV0RMin", 0., "minimum conversion radius of the V0s"};
  Configurable<float> fV0RMax{"fV0RMax", 180., "maximum conversion radius of the V0s"};
  Configurable<float> LineCutZ0{"fLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
  Configurable<float> LineCutZRSlope{"LineCutZRSlope", std::tan(2.f * std::atan(std::exp(-fEtaMax.value))), "The slope for the line cut"};

  HistogramRegistry registry{
    "registry",
    {
      gHistoSpec_hCollisionZ_all_MCTrue,
      gHistoSpec_hCollisionZ_MCTrue,
      gHistoSpec_hCollisionZ_MCRec,

      {"hGammaProdAfterCutsP_MCTrue", "hGammaProdAfterCutsP_MCTrue;p (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},
      {"hGammaProdAfterCutsPt_MCTrue", "hGammaProdAfterCutsPt_MCTrue;p_T (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},

      {"hGammaConvertedP_Rsel_MCTrue", "hGammaConvertedP_Rsel_MCTrue;p (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},
      {"hGammaConvertedPt_Rsel_MCTrue", "hGammaConvertedPt_Rsel_MCTrue;p_T (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},
      {"hGammaConvertedR_MCTrue", "hGammaConvertedR_MCTrue;conversion radius (cm);counts", {HistType::kTH1F, {{2000, 0.f, 500.f}}}},

      {"hGammaConvertedEtaP_MCTrue", "hGammaConvertedEtaP_MCTrue;#eta;p (GeV/c)", {HistType::kTH2F, {gAxis_eta2d, gAxis_pT2d}}},
      {"hGammaConvertedEtaR_MCTrue", "hGammaConvertedEtaR_MCTrue;#eta;conversion radius (cm)", {HistType::kTH2F, {gAxis_eta2d, gAxis_r2d}}},
      {"hGammaConvertedEtaZ_MCTrue", "hGammaConvertedEtaZ_MCTrue;#eta;conversion z (cm)", {HistType::kTH2F, {gAxis_eta2d, gAxis_z2d}}},
      {"hGammaConvertedRP_MCTrue", "hGammaConvertedRP_MCTrue;conversion radius (cm);conversion z (cm)", {HistType::kTH2F, {gAxis_r2d, gAxis_pT2d}}},
      {"hGammaConvertedRZ_MCTrue", "hGammaConvertedRZ_MCTrue;conversion radius (cm);conversion z (cm)", {HistType::kTH2F, {gAxis_r2d, gAxis_z2d}}},
      {"hGammaConvertedRPt_MCTrue", "hGammaConvertedRPt_MCTrue;conversion radius (cm);p_T (GeV/c)", {HistType::kTH2F, {gAxis_r2d, gAxis_pT2d}}},
      {"hGammaConvertedXY_MCTrue", "hGammaConvertedXY_MCTrue;conversion x (cm);conversion y (cm)", {HistType::kTH2F, {gAxis_z2d, gAxis_z2d}}},
      {"hGammaConvertedZP_MCTrue", "hGammaConvertedZP_MCTrue;conversion z (cm);p (GeV/c)", {HistType::kTH2F, {gAxis_z2d, gAxis_pT2d}}},
      {"hGammaConvertedpeDivpGamma", "hpeDivpGamma;p (GeV/c);p_{e}/p_{#gamma};counts", {HistType::kTH2F, {gAxis_pT, {220, 0.f, 1.1f}}}},

      // debugging histograms
      {"hNDaughters_MCTrue", "hNDaughters_MCTrue;nDaughters;counts", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
    },
  };

  template <typename MCGAMMA>
  void fillConversionHistograms(MCGAMMA const& theMcConvGamma)
  {
    float lConversionRadius = theMcConvGamma.v0Radius();
    // 1d histos
    registry.fill(HIST("hGammaConvertedR_MCTrue"), lConversionRadius);
    // 2d histos
    registry.fill(HIST("hGammaConvertedEtaP_MCTrue"), theMcConvGamma.eta(), theMcConvGamma.p());
    registry.fill(HIST("hGammaConvertedEtaR_MCTrue"), theMcConvGamma.eta(), lConversionRadius);
    registry.fill(HIST("hGammaConvertedEtaZ_MCTrue"), theMcConvGamma.eta(), theMcConvGamma.conversionZ());
    registry.fill(HIST("hGammaConvertedRP_MCTrue"), lConversionRadius, theMcConvGamma.p());
    registry.fill(HIST("hGammaConvertedRPt_MCTrue"), lConversionRadius, theMcConvGamma.pt());
    registry.fill(HIST("hGammaConvertedRZ_MCTrue"), lConversionRadius, theMcConvGamma.conversionZ());
    registry.fill(HIST("hGammaConvertedXY_MCTrue"), theMcConvGamma.conversionX(), theMcConvGamma.conversionY());
    registry.fill(HIST("hGammaConvertedZP_MCTrue"), theMcConvGamma.conversionZ(), theMcConvGamma.p());

    if (lConversionRadius > fV0RMin && lConversionRadius < fV0RMax) {
      registry.fill(HIST("hGammaConvertedP_Rsel_MCTrue"), theMcConvGamma.p());
      registry.fill(HIST("hGammaConvertedPt_Rsel_MCTrue"), theMcConvGamma.pt());
    }
  }
  template <typename MCGAMMA, typename MCDAUONE, typename MCDAUTWO>
  void fillAsymmetryHistograms(MCGAMMA const& theMcConvGamma, MCDAUONE const theFirstDaughter, MCDAUTWO theSecondDaughter)
  {
    float lConversionRadius = theMcConvGamma.v0Radius();
    float lGammaMomentum = theMcConvGamma.p();

    if (lConversionRadius > fV0RMin && lConversionRadius < fV0RMax) {
      registry.fill(HIST("hGammaConvertedpeDivpGamma"), lGammaMomentum, theFirstDaughter.p() / lGammaMomentum);
      registry.fill(HIST("hGammaConvertedpeDivpGamma"), lGammaMomentum, theSecondDaughter.p() / lGammaMomentum);
    }
  }

  template <typename MCGAMMA>
  bool photonPassesCuts(MCGAMMA const& theMcGamma)
  {
    if (fPhysicalPrimaryOnly && !theMcGamma.isPhysicalPrimary()) {
      // fill histo
      return false;
    }

    if (std::abs(theMcGamma.eta()) > fEtaMax) {
      // fill histo
      return false;
    }

    if (std::abs(theMcGamma.conversionZ()) > LineCutZ0 + theMcGamma.v0Radius() * LineCutZRSlope) {
      return false;
    }
    return true;
  }

  // loop over MC truth McCollisions
  void process(aod::McCollision const& theMcCollision,
               soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                          aod::Collisions>> const& theCollisions,
               aod::McGammasTrue const& theMcGammas,
               aod::McDaughterTrue const&)
  {
    registry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());
    if (theCollisions.size() == 0) {
      return;
    }
    registry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());

    for (const auto& lCollision : theCollisions) {
      registry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());
    }

    for (const auto& lMcGamma : theMcGammas) {

      if (!photonPassesCuts(lMcGamma)) {
        continue;
      }

      registry.fill(HIST("hGammaProdAfterCutsP_MCTrue"), lMcGamma.p());
      registry.fill(HIST("hGammaProdAfterCutsPt_MCTrue"), lMcGamma.pt());

      int const lNDaughters = lMcGamma.nDaughters();
      registry.fill(HIST("hNDaughters_MCTrue"), 0.5 + lNDaughters);
      if (lNDaughters == 2) {

        fillConversionHistograms(lMcGamma);

        if (lMcGamma.has_mcDaughterTrueOne() && lMcGamma.has_mcDaughterTrueTwo()) {
          auto lMcDaughterOne = lMcGamma.mcDaughterTrueOne();
          auto lMcDaughterTwo = lMcGamma.mcDaughterTrueTwo();

          fillAsymmetryHistograms(lMcGamma, lMcDaughterOne, lMcDaughterTwo);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsTruthOnlyMc>(cfgc)};
}
