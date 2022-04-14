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

/// \brief extract relevant mc truth information that allows to compute efficiency, purity and more quantities for the photon conversion analysis.
/// dependencies: none
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using tracksAndMcLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;

struct gammaConversionsRubenDirect {

  HistogramRegistry registry{
    "registry",
    {
      {"hGammaReconctructedPtTrue", "hGammaReconctructedPtTrue", {HistType::kTH1F, {{800, -0.f, 25.f}}}},
    },
  };

  void process(aod::Collisions::iterator const& theCollision,
               aod::V0s const& theV0s,
               tracksAndMcLabels const& theTracks,
               aod::McParticles const& theMcParticles)
  {
    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndMcLabels>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndMcLabels>(); // negative daughter

      auto lMcPos = lTrackPos.template mcParticle_as<aod::McParticles_001>();
      auto lMcNeg = lTrackNeg.template mcParticle_as<aod::McParticles_001>();

      //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
      std::vector<int> lMothers;
      for (auto& mP : lMcPos.template mothers_as<aod::McParticles_001>()) {
        LOGF(info, "   mother index mP: %d", mP.globalIndex());
        lMothers.push_back(mP.globalIndex());
      }

      if (lMothers.size() > 0) {
        for (auto& mN : lMcNeg.template mothers_as<aod::McParticles_001>()) {
          LOGF(info, "   mother index mN: %d", mN.globalIndex());
          lMothers.push_back(mN.globalIndex());
        }
      }

      // SFS verify theyre all the same category and remove
      int lSame = 0;
      {
        if (lMothers.size() == 2) {
          if (lMothers[0] == lMothers[1]) {
            LOGF(info, "size2: 01");
            lSame = 1;
          }
        }

        if (lMothers.size() == 3) {
          if (lMothers[0] == lMothers[1]) {
            LOGF(info, "size3: 01");
            lSame = 2;
          }
          if (lMothers[0] == lMothers[2]) {
            LOGF(info, "size2: 02");
            lSame = 3;
          }
          if (lMothers[1] == lMothers[2]) {
            LOGF(info, "size2: 12");
            lSame = 4;
          }
        }

        if (lMothers.size() == 4) {
          if (lMothers[0] == lMothers[2]) {
            LOGF(info, "size4 02");
            lSame = 4;
          }
          if (lMothers[1] == lMothers[3]) {
            LOGF(info, "size4 13");
            lSame = 5;
          }
          if (lMothers[0] == lMothers[3]) {
            LOGF(info, "size4 03");
            lSame = 6;
          }
          if (lMothers[1] == lMothers[2]) {
            LOGF(info, "size4 12");
            lSame = 7;
          }
          if (lMothers[0] == lMothers[1]) {
            LOGF(info, "size4 01");
            lSame = 8;
          }
          if (lMothers[2] == lMothers[3]) {
            LOGF(info, "size4 23");
            lSame = 9;
          }
        }
      }

      if (lSame) {

        // SFS todo: actually no loop required here, for this
        for (auto& lMother : lMcNeg.template mothers_as<aod::McParticles_001>()) {

          if (lMother.pdgCode() == 22) {
            registry.fill(HIST("hGammaReconctructedPtTrue"), lMother.pt());
          }
        }
      }
    }
  }
};

struct gammaConversionsMcTruthOnly {

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionZ", "hCollisionZ", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hGammaProdInAcc", "hGammaProdInAcc", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaMoreThanTwoDaughtersPt", "hGammaMoreThanTwoDaughtersPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaConvertedPt", "hGammaConvertedPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaConvertedR", "hGammaConvertedR", {HistType::kTH1F, {{1600, 0.f, 500.f}}}},
      {"hGammaConvertedRPt", "hGammaConvertedRPt", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
    },
  };

  // loop over MC truth McCollisions
  void process(aod::McCollision const& mcCollision,
               soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                          aod::Collisions>> const& collisions,
               aod::McParticles const& theMcParticles)
  {
    for (auto& collision : collisions) {

      registry.fill(HIST("hCollisionZ"), collision.posZ());
      for (auto& lMcParticle : theMcParticles) {

        if ((lMcParticle.pdgCode() == 22) &&
            (std::abs(lMcParticle.eta()) < 0.8)) { // SFS todo: track v0 eta??

          registry.fill(HIST("hGammaProdInAcc"), lMcParticle.pt());

          size_t lBothElectrons = 0;
          if (lMcParticle.has_daughters()) {
            for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
              lBothElectrons += std::abs(lDaughter.pdgCode()) == 11;
            }
            if (lBothElectrons == 2) {
              for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
                float lConversionRadius = std::sqrt(std::pow(lDaughter.vx(), 2) + std::pow(lDaughter.vy(), 2));
                registry.fill(HIST("hGammaConvertedR"), lConversionRadius);
                registry.fill(HIST("hGammaConvertedRPt"), lConversionRadius, lMcParticle.pt());

                if (lConversionRadius > 5. && lConversionRadius < 180.) {
                  registry.fill(HIST("hGammaConvertedPt"), lMcParticle.pt());
                }
                break;
              }
            }
            if (lBothElectrons > 2) {
              registry.fill(HIST("hGammaMoreThanTwoDaughtersPt"), lMcParticle.pt());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsMcTruthOnly>(cfgc),
                      adaptAnalysisTask<gammaConversionsRubenDirect>(cfgc)};
}
