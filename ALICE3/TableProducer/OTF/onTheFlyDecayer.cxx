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

///
/// \file onTheFlyDecayer.cxx
/// \brief pre-processing for on-the-fly analysis
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>
///

#include "ALICE3/Core/Decayer.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFMCParticle.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <DetectorsVertexing/PVertexerHelpers.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <SimulationDataFormat/InteractionSampler.h>

#include <TPDGCode.h>

#include <string>
#include <gsl/span>
#include <span>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr int kNumDecays = 7;
static constexpr int kNumParameters = 1;
static constexpr int defaultParameters[kNumDecays][kNumParameters]{{1}, {1}, {1}, {1}, {1}, {1}, {1}};
static const std::vector<std::string> parameterNames{"enable"};
static const std::vector<std::string> particleNames{"K0s",
                                                    "Lambda",
                                                    "Anti-Lambda",
                                                    "Xi",
                                                    "Anti-Xi",
                                                    "Omega",
                                                    "Anti-Omega"};
static const std::vector<int> pdgCodes{kK0Short,
                                       kLambda0,
                                       kLambda0Bar,
                                       kXiMinus,
                                       kXiPlusBar,
                                       kOmegaMinus,
                                       kOmegaPlusBar};

struct OnTheFlyDecayer {
  Produces<aod::McParticlesWithDau> tableMcParticlesWithDau;
  o2::upgrade::Decayer decayer;
  Service<o2::framework::O2DatabasePDG> pdgDB;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<int64_t, std::vector<o2::upgrade::OTFParticle>> mDecayDaughters;

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
    {defaultParameters[0], kNumDecays, kNumParameters, particleNames, parameterNames},
    "Enable option for particle to be decayed: 0 - no, 1 - yes"};


  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};


  std::vector<int> mEnabledDecays;
  void init(o2::framework::InitContext&)
  {
	  decayer.setSeed(seed);
    decayer.setBField(magneticField);
    for (int i = 0; i < kNumDecays; ++i) {
      if (enabledDecays->get(particleNames[i].c_str(), "enable")) {
        mEnabledDecays.push_back(pdgCodes[i]);
      }
    }

    histos.add("K0S/hGenK0S", "hGenK0S", kTH1D, {axisPt});
    histos.add("Lambda/hGenLambda", "hGenLambda", kTH1D, {axisPt});
    histos.add("AntiLambda/hGenAntiLambda", "hGenAntiLambda", kTH1D, {axisPt});
    histos.add("Xi/hGenXi", "hGenXi", kTH1D, {axisPt});
    histos.add("AntiXi/hGenAntiXi", "hGenAntiXi", kTH1D, {axisPt});
    histos.add("Omega/hGenOmega", "hGenOmega", kTH1D, {axisPt});
    histos.add("AntiOmega/hGenAntiOmega", "hGenAntiOmega", kTH1D, {axisPt});
  }

  bool canDecay(const int pdgCode)
  {
    return std::find(mEnabledDecays.begin(), mEnabledDecays.end(), pdgCode) != mEnabledDecays.end();
  } 

  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    mDecayDaughters.clear();
    u_int64_t nStoredDaughters = 0;
    for (int64_t index{0}; index < mcParticles.size(); ++index) {
      const auto& particle = mcParticles.iteratorAt(index);
      std::vector<o2::upgrade::OTFParticle> decayDaughters;
      static constexpr int kMaxNestedDecays = 10;
      int nDecays = 0;
      if (canDecay(particle.pdgCode())) {
        o2::track::TrackParCov o2track;
        o2::upgrade::convertMCParticleToO2Track(particle, o2track, pdgDB);
        decayDaughters = decayer.decayParticle(pdgDB, o2track, particle.pdgCode());
        for (size_t idau{0}; idau < decayDaughters.size(); ++idau) {
          const auto& dau = decayDaughters[idau];
          o2::track::TrackParCov dauTrack;
          o2::upgrade::convertOTFParticleToO2Track(dau, dauTrack, pdgDB);
          if (canDecay(dau.pdgCode())) {
            std::vector<o2::upgrade::OTFParticle> cascadingDaughers = decayer.decayParticle(pdgDB, dauTrack, dau.pdgCode());
            for (const auto& daudau : cascadingDaughers) {
              decayDaughters.push_back(daudau);
              if (kMaxNestedDecays < ++nDecays) {
                LOG(error) << "Seemingly stuck trying to perpetually decay products from pdg: " << particle.pdgCode();
              }
            }
          }
        }

        if (decayDaughters.empty()) {
          LOG(error) << "Attempted to decay " << particle.pdgCode() << " but resulting vector of daugthers were empty";
          continue;
        }
      }
      
      
      if (particle.pdgCode() == kK0Short) {
        histos.fill(HIST("K0S/hGenK0S"), particle.pt());
      } else if (particle.pdgCode() == kLambda0) {
        histos.fill(HIST("Lambda/hGenLambda"), particle.pt());
      } else if (particle.pdgCode() == kLambda0Bar) {
        histos.fill(HIST("AntiLambda/hGenAntiLambda"), particle.pt());
      } else if (particle.pdgCode() == kXiMinus) {
        histos.fill(HIST("Xi/hGenXi"), particle.pt());
      } else if (particle.pdgCode() == kXiPlusBar) {
        histos.fill(HIST("AntiXi/hGenAntiXi"), particle.pt());
      } else if (particle.pdgCode() == kOmegaMinus) {
        histos.fill(HIST("Omega/hGenOmega"), particle.pt());
      } else if (particle.pdgCode() == kOmegaPlusBar) {
        histos.fill(HIST("AntiOmega/hGenAntiOmega"), particle.pt());
      }

      int daughtersIdSlice[2];
      if (canDecay(particle.pdgCode())) {
        daughtersIdSlice[0] = static_cast<int>(mcParticles.size() + nStoredDaughters);
        daughtersIdSlice[1] = static_cast<int>(mcParticles.size() + nStoredDaughters + decayDaughters.size());
      } else {
        daughtersIdSlice[0] = static_cast<int>(particle.daughtersIds()[0]);
        daughtersIdSlice[1] = static_cast<int>(particle.daughtersIds()[1]);
      }

      std::span<const int> motherSpan(particle.mothersIds().data(), particle.mothersIds().size());
      mDecayDaughters.emplace(index, decayDaughters);
      nStoredDaughters += decayDaughters.size();

      float phi = o2::constants::math::PI + std::atan2(-1.0f * particle.py(), -1.0f * particle.px());
      float eta; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
      float pt = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py());
      float p = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py() + particle.pz() * particle.pz());
      float y; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

      if ((p - particle.pz()) < 1e-7f) {
        eta = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        eta = 0.5f * std::log((p + particle.pz()) / (p - particle.pz()));
      }

      if ((particle.e() - particle.pz()) < 1e-7f) {
        y = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        y = 0.5f * std::log((particle.e() + particle.pz()) / (particle.e() - particle.pz()));
      }

      // TODO: Particle status code
      // TODO: Expression columns
      tableMcParticlesWithDau(particle.mcCollisionId(), particle.pdgCode(), particle.statusCode(), 
                              particle.flags(), motherSpan, daughtersIdSlice, particle.weight(),
                              particle.px(), particle.py(), particle.pz(), particle.e(), 
                              particle.vx(), particle.vy(), particle.vz(), particle.vt(),
                              phi, eta, pt, p, y);
    }

    int daughtersIdSlice[2] = {-1, -1};
    for (const auto& [index, decayDaughters] : mDecayDaughters) {
      for (const auto& dau : decayDaughters) {
        if (index >= mcParticles.size()) {
          LOG(warn) << "--- Index " << index << " out of bounds for mcParticles table of size " << mcParticles.size();
          continue;
        }

        auto mother = mcParticles.iteratorAt(index);
        std::vector<int> motherIds = { static_cast<int>(index) };
        std::span<const int> motherSpan(motherIds.data(), motherIds.size());

        float phi = o2::constants::math::PI + std::atan2(-1.0f * dau.py(), -1.0f * dau.px());
        float eta; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
        float pt = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py());
        float p = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py() + dau.pz() * dau.pz());
        float y; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

        if ((p - dau.pz()) < 1e-7f) {
          eta = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          eta = 0.5f * std::log((p + dau.pz()) / (p - dau.pz()));
        }

        if ((dau.e() - dau.pz()) < 1e-7f) {
          y = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          y = 0.5f * std::log((dau.e() + dau.pz()) / (dau.e() - dau.pz()));
        }


        // TODO: Particle status code
        // TODO: Expression columns
        tableMcParticlesWithDau(mother.mcCollisionId(), dau.pdgCode(), 1,
                                mother.flags(), motherSpan, daughtersIdSlice, mother.weight(), 
                                dau.px(), dau.py(), dau.pz(), dau.e(),
                                dau.vx(), dau.vy(), dau.vz(), mother.vt(),
                                phi, eta, pt, p, y);
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}