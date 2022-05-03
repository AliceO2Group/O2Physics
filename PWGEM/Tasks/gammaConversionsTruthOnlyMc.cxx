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

#include "gammaTables.h"

#include "TVector3.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct gammaConversionsTruthOnlyMc {

  Configurable<bool> fOnlyPrimary{"fOnlyPrimary", true, "fOnlyPrimary"};
  Configurable<float> fEtaMax{"fEtaMax", 0.8, "aMaximum photon eta"};

  HistogramRegistry registry{
    "registry",
    {
      {"hMcCollisionZ_MCTrue", "hMcCollisionZ_MCTrue", {HistType::kTH1F, {{400, -50.f, 50.f}}}},
      {"hEtaDiff_MCTrue", "hEtaDiff_MCTrue", {HistType::kTH1F, {{400, -2.f, 2.f}}}},

      {"hNDaughters_MCTrue", "hNDaughters_MCTrue", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      {"hPdgCodeDaughters_MCTrue", "hPdgCodeDaughters_MCTrue", {HistType::kTH1F, {{2000, -1000.f, 1000.f}}}},
      {"hNElectrons_MCTrue", "hNElectrons_MCTrue", {HistType::kTH1F, {{50, 0.f, 50.f}}}},

      {"hGammaProdAfterCutsP_MCTrue", "hGammaProdAfterCutsP_MCTrue", {HistType::kTH1F, {{800, 0.f, 25.f}}}},

      {"hGammaConvertedR_MCTrue", "hGammaConvertedR_MCTrue", {HistType::kTH1F, {{1600, 0.f, 500.f}}}},
      {"hGammaConvertedRselP_MCTrue", "hGammaConvertedRselP_MCTrue", {HistType::kTH1F, {{800, 0.f, 25.f}}}},

      {"hGammaConvertedEtaP_MCTrue", "hGammaConvertedEtaP_MCTrue", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedEtaR_MCTrue", "hGammaConvertedEtaR_MCTrue", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, 0.f, 250.f}}}},
      {"hGammaConvertedEtaZ_MCTrue", "hGammaConvertedEtaZ_MCTrue", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, -250.f, 250.f}}}},

      {"hGammaConvertedRP_MCTrue", "hGammaConvertedRP_MCTrue", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedRZ_MCTrue", "hGammaConvertedRZ_MCTrue", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, -250.f, 250.f}}}},

      {"hGammaConvertedZP_MCTrue", "hGammaConvertedZP_MCTrue", {HistType::kTH2F, {{400, -250.f, 250.f}, {400, 0.f, 25.f}}}},

      {"hGammaProdAfterCutsPt_MCTrue", "hGammaProdAfterCutsPt_MCTrue", {HistType::kTH1F, {{800, 0.f, 25.f}}}},

      {"hGammaConvertedRPt_MCTrue", "hGammaConvertedRPt_MCTrue", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedRselPt_MCTrue", "hGammaConvertedRselPt_MCTrue", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hPeculiarOccurences_MCTrue", "hPeculiarOccurences_MCTrue", {HistType::kTH1F, {{50, -25.f, 25.f}}}},
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
    registry.fill(HIST("hGammaConvertedEtaP_MCTrue"), theMcConvGamma.eta(), theMcConvGamma.p());
    registry.fill(HIST("hGammaConvertedEtaR_MCTrue"), theMcConvGamma.eta(), lConversionRadius);
    registry.fill(HIST("hGammaConvertedEtaZ_MCTrue"), theMcConvGamma.eta(), theMcConvGamma.conversionZ());
    registry.fill(HIST("hGammaConvertedR_MCTrue"), lConversionRadius);
    registry.fill(HIST("hGammaConvertedRP_MCTrue"), lConversionRadius, theMcConvGamma.p());
    registry.fill(HIST("hGammaConvertedRPt_MCTrue"), lConversionRadius, theMcConvGamma.pt());
    registry.fill(HIST("hGammaConvertedRZ_MCTrue"), lConversionRadius, theMcConvGamma.conversionZ());
    registry.fill(HIST("hGammaConvertedZP_MCTrue"), theMcConvGamma.conversionZ(), theMcConvGamma.p());

    TVector3 lDaughter0Vtx(theMcConvGamma.conversionX(), theMcConvGamma.conversionY(), theMcConvGamma.conversionZ());
    float_t lEtaDiff = lDaughter0Vtx.Eta() - theMcConvGamma.eta();
    registry.fill(HIST("hEtaDiff_MCTrue"), lEtaDiff);

    if (lConversionRadius > 5. && lConversionRadius < 180.) {
      registry.fill(HIST("hGammaConvertedRselP_MCTrue"), theMcConvGamma.p());
      registry.fill(HIST("hGammaConvertedRselPt_MCTrue"), theMcConvGamma.pt());
    }
  }

  template <typename MCGAMMA>
  bool photonPassesCuts(MCGAMMA const& theMcGamma)
  {
    if (fOnlyPrimary && !theMcGamma.isPhysicalPrimary()) {
      // fill histo
      return false;
    }

    if (std::abs(theMcGamma.eta()) >= fEtaMax) { // SFS todo: track v0 eta??
      // fill histo
      return false;
    }
    return true;
  }

  // loop over MC truth McCollisions
  void process(aod::McCollision const& theMcCollision,
               aod::McGammasTrue const& theMcGammas,
               aod::McGammaDaughtersTrue const& theMcGammaDaughters)
  {
    registry.fill(HIST("hMcCollisionZ_MCTrue"), theMcCollision.posZ());

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
