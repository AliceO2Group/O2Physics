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

/// \file taskMultiplicityEstimatorCorrelation.cxx
/// \brief Task for correlating the multiplicity estimator with generated dN/deta
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskMultiplicityEstimatorCorrelation {
  HistogramRegistry registry{"registry", {}};
  static constexpr int8_t NEstimators = 8;
  static constexpr std::array<std::string_view, NEstimators> EstimatorsNames = {"FV0A", "FT0A", "FT0C", "FT0M", "FDDA", "FDDC", "FDDM", "NTPV"};

  std::vector<unsigned> consideredParticles = {
    kElectron,
    kMuonMinus,
    kPiPlus,
    kKPlus,
    kProton};

  ConfigurableAxis axisFV0A = {"axisFV0A", {100, 0., 20000.}, "axis for FV0A estimator"};
  ConfigurableAxis axisFT0A = {"axisFT0A", {100, 0., 10000.}, "axis for FT0A estimator"};
  ConfigurableAxis axisFT0C = {"axisFT0C", {100, 0., 5000.}, "axis for FT0C estimator"};
  ConfigurableAxis axisFT0M = {"axisFT0M", {100, 0., 10000.}, "axis for FT0M estimator"};
  ConfigurableAxis axisFDDA = {"axisFDDA", {100, 0., 20000.}, "axis for FDDA estimator"};
  ConfigurableAxis axisFDDC = {"axisFDDC", {100, 0., 5000.}, "axis for FDDC estimator"};
  ConfigurableAxis axisFDDM = {"axisFDDM", {100, 0., 20000.}, "axis for FDDM estimator"};
  ConfigurableAxis axisNTPV = {"axisNTPV", {100, 0., 100.}, "axis for NTPV estimator"};
  ConfigurableAxis axisdNdEta = {"axisdNdEta", {100, 0., 100.}, "axis for dN/deta"};

  std::vector<ConfigurableAxis*> estimatorsAxes = {&axisFV0A, &axisFT0A, &axisFT0C, &axisFT0M, &axisFDDA, &axisFDDC, &axisFDDM, &axisNTPV};

  Preslice<aod::McParticles> particlesPerCollision = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  using CollisionsWithMult = soa::Join<aod::Collisions, aod::PVMults, aod::MultZeqs, aod::EvSels, aod::McCollisionLabels>;

  void init(InitContext&)
  {
    for (int8_t i = 0; i < NEstimators; i++) {
      registry.add<TH2>(("etaPFive/" + std::string(EstimatorsNames[i]) + "VsdNdeta").c_str(), (std::string(EstimatorsNames[i]) + "VsdNdeta;" + std::string(EstimatorsNames[i]) + ";<dN_{ch}/d#eta>").c_str(), HistType::kTH2F, {*(estimatorsAxes[i]), axisdNdEta});
      registry.add<TH2>(("etaOne/" + std::string(EstimatorsNames[i]) + "VsdNdeta").c_str(), (std::string(EstimatorsNames[i]) + "VsdNdeta;" + std::string(EstimatorsNames[i]) + ";<dN_{ch}/d#eta>").c_str(), HistType::kTH2F, {*(estimatorsAxes[i]), axisdNdEta});
    }
  }

  void process(CollisionsWithMult const& collisions,
               aod::McCollisions const& mcCollisions,
               aod::McParticles const& particles,
               soa::Join<aod::BCs, aod::Timestamps> const&)
  {
    for (auto const& collision : mcCollisions) {

      // Get multiplicity for the reconstructed collision with the highest number of contributors
      unsigned maxNumContrib = 0;
      CollisionsWithMult::iterator collisionMaxNumContrib;
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, collision.globalIndex());
      for (const auto& recCol : recoCollsPerMcColl) {
        if (recCol.numContrib() > maxNumContrib) {
          maxNumContrib = recCol.numContrib();
          collisionMaxNumContrib = recCol;
        }
      }
      std::vector<float> multiplicity = {
        collisionMaxNumContrib.multZeqFV0A(),
        collisionMaxNumContrib.multZeqFT0A(),
        collisionMaxNumContrib.multZeqFT0C(),
        collisionMaxNumContrib.multZeqFT0A() + collisionMaxNumContrib.multZeqFT0C(),
        collisionMaxNumContrib.multZeqFDDA(),
        collisionMaxNumContrib.multZeqFDDC(),
        collisionMaxNumContrib.multZeqFDDA() + collisionMaxNumContrib.multZeqFDDC(),
        collisionMaxNumContrib.multZeqNTracksPV()};

      // Get the dN/deta for the generated collision
      unsigned nChargedInEtaFive = 0;
      unsigned nChargedInEtaOne = 0;
      const auto& particlesPerMcColl = particles.sliceBy(particlesPerCollision, collision.globalIndex());
      for (auto const& particle : particlesPerMcColl) {
        if (particle.isPhysicalPrimary()) {
          bool isCharged = false;
          for (auto const& consideredParticle : consideredParticles) {
            if (static_cast<unsigned int>(std::abs(particle.pdgCode())) == consideredParticle) {
              isCharged = true;
              break;
            }
          }
          if (!isCharged) {
            continue;
          }
          if (std::abs(particle.eta()) < 0.5) {
            nChargedInEtaFive++;
          }
          if (std::abs(particle.eta()) < 1.0) {
            nChargedInEtaOne++;
          }
        }
      }

      float dNdetaFive = nChargedInEtaFive;
      float dNdetaOne = nChargedInEtaOne / 2.0;
      for (int i = 0; i < NEstimators; i++) {
        static_for<0, NEstimators - 1>([&](auto j) {
          constexpr int Index = j.value;
          registry.fill(HIST("etaPFive/") + HIST(EstimatorsNames[Index]) + HIST("VsdNdeta"), multiplicity[Index], dNdetaFive);
          registry.fill(HIST("etaOne/") + HIST(EstimatorsNames[Index]) + HIST("VsdNdeta"), multiplicity[Index], dNdetaOne);
        });
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMultiplicityEstimatorCorrelation>(cfgc)};
}
