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

// SFS todo: update description
/// \brief extract relevant mc truth information that allows to compute efficiency, purity and more quantities for the photon conversion analysis.
/// dependencies: none
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"

#include "TVector3.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct gammaConversionsTruthOnlyMc {

  Configurable<bool> fPhysicalPrimaryOnly{"fPhysicalPrimaryOnly", true, "fPhysicalPrimaryOnly"};
  Configurable<float> fEtaMax{"fEtaMax", 0.8, "aMaximum photon eta"};

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
      {"hGammaConvertedR_MCTrue", "hGammaConvertedR_MCTrue;conversion radius (cm);counts", {HistType::kTH1F, {{1600, 0.f, 500.f}}}},

      {"hGammaConvertedEtaP_MCTrue", "hGammaConvertedEtaP_MCTrue;#eta;p (GeV/c)", {HistType::kTH2F, {gAxis_eta2d, gAxis_pT2d}}},
      {"hGammaConvertedEtaR_MCTrue", "hGammaConvertedEtaR_MCTrue;#eta;conversion radius (cm)", {HistType::kTH2F, {gAxis_eta2d, gAxis_r2d}}},
      {"hGammaConvertedEtaZ_MCTrue", "hGammaConvertedEtaZ_MCTrue;#eta;conversion z (cm)", {HistType::kTH2F, {gAxis_eta2d, gAxis_z2d}}},
      {"hGammaConvertedRP_MCTrue", "hGammaConvertedRP_MCTrue;conversion radius (cm);conversion z (cm)", {HistType::kTH2F, {gAxis_r2d, gAxis_pT2d}}},
      {"hGammaConvertedRZ_MCTrue", "hGammaConvertedRZ_MCTrue;conversion radius (cm);conversion z (cm)", {HistType::kTH2F, {gAxis_r2d, gAxis_z2d}}},
      {"hGammaConvertedRPt_MCTrue", "hGammaConvertedRPt_MCTrue;conversion radius (cm);p_T (GeV/c)", {HistType::kTH2F, {gAxis_r2d, gAxis_pT2d}}},
      {"hGammaConvertedXY_MCTrue", "hGammaConvertedXY_MCTrue;conversion x (cm);conversion y (cm)", {HistType::kTH2F, {gAxis_z2d, gAxis_z2d}}},
      {"hGammaConvertedZP_MCTrue", "hGammaConvertedZP_MCTrue;conversion z (cm);p (GeV/c)", {HistType::kTH2F, {gAxis_z2d, gAxis_pT2d}}},

      // debugging histograms
      {"hPeculiarOccurences_MCTrue", "hPeculiarOccurences_MCTrue", {HistType::kTH1F, {{50, -25.f, 25.f}}}},
      {"hNElectrons_MCTrue", "hNElectrons_MCTrue", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      {"hNDaughters_MCTrue", "hNDaughters_MCTrue;nDaughters;counts", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      {"hPdgCodeDaughters_MCTrue", "hPdgCodeDaughters_MCTrue;pdg code;counts", {HistType::kTH1F, {{2000, -1000.f, 1000.f}}}},
    },
  };

  template <typename MCGAMMA>
  void fillConversionHistograms(MCGAMMA const& theMcConvGamma)
  {
    // this produces another convpoint r distribution than the uncommend version below. why?
    // access first daughter to get conversion point

    /*auto const &lDaughter0 = lDaughters.begin();
    float lConversionRadius = std::sqrt(std::pow(lDaughter0.vx(), 2) + std::pow(lDaughter0.vy(), 2));
    registry.fill(HIST("hGammaConvertedEtaP"), lMcGamma.eta(), lMcGamma.p());
    registry.fill(HIST("hGammaConvertedEtaR"), lMcGamma.eta(), lConversionRadius);
    registry.fill(HIST("hGammaConvertedEtaZ"), lMcGamma.eta(), lDaughter0.vz());
    registry.fill(HIST("hGammaConvertedR"), lConversionRadius);
    registry.fill(HIST("hGammaConvertedRP"), lConversionRadius, lMcGamma.p());
    registry.fill(HIST("hGammaConvertedRZ"), lConversionRadius, lDaughter0.vz());
    registry.fill(HIST("hGammaConvertedZP"), lDaughter0.vz(), lMcGamma.p());
    registry.fill(HIST("hGammaConvertedRPt"), lConversionRadius, lMcGamma.pt());
    TVector3 lDaughter0Vtx(lDaughter0.vx(),lDaughter0.vy(), lDaughter0.vz());
    float_t lEtaDiff = lDaughter0Vtx.Eta() - lMcGamma.eta();
    registry.fill(HIST("hEtaDiff"), lEtaDiff);*/

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

    if (lConversionRadius > 5. && lConversionRadius < 180.) {
      registry.fill(HIST("hGammaConvertedP_Rsel_MCTrue"), theMcConvGamma.p());
      registry.fill(HIST("hGammaConvertedPt_Rsel_MCTrue"), theMcConvGamma.pt());
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
    return true;
  }

  // loop over MC truth McCollisions
  void process(aod::McCollision const& theMcCollision,
               soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                          aod::Collisions>> const& theCollisions,
               aod::McGammasTrue const& theMcGammas,
               aod::McGammaDaughtersTrue const& theMcGammaDaughters)
  {
    registry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());
    if (theCollisions.size() == 0) {
      return;
    }
    registry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());

    for (auto& lCollision : theCollisions) {
      registry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());
    }

    for (auto& lMcGamma : theMcGammas) {

      if (!photonPassesCuts(lMcGamma)) {
        continue;
      }

      registry.fill(HIST("hGammaProdAfterCutsP_MCTrue"), lMcGamma.p());
      registry.fill(HIST("hGammaProdAfterCutsPt_MCTrue"), lMcGamma.pt());

      int const lNDaughters = lMcGamma.nDaughters();
      registry.fill(HIST("hNDaughters_MCTrue"), 0.5 + lNDaughters);
      if (lNDaughters == 2) {

        fillConversionHistograms(lMcGamma);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsTruthOnlyMc>(cfgc)};
}
