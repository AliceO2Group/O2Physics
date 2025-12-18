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
  Produces<aod::OTFMCParticles> tableMcParticles;
  o2::upgrade::Decayer decayer;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
    {defaultParameters[0], kNumDecays, kNumParameters, particleNames, parameterNames},
    "Enable option for particle to be decayed: 0 - no, 1 - yes"};


  ConfigurableAxis axisRadius{"axisRadius", {1000, 0, 100}, "Radius"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};


  std::vector<int> mEnabledDecays;
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);
	  decayer.setSeed(seed);
    decayer.setBField(magneticField);
    for (int i = 0; i < kNumDecays; ++i) {
      if (enabledDecays->get(particleNames[i].c_str(), "enable")) {
        mEnabledDecays.push_back(pdgCodes[i]);
      }
    }

    histos.add("K0S/hGenK0S", "hGenK0S", kTH2D, {axisRadius, axisPt});
    histos.add("Lambda/hGenLambda", "hGenLambda", kTH2D, {axisRadius, axisPt});
    histos.add("AntiLambda/hGenAntiLambda", "hGenAntiLambda", kTH2D, {axisRadius, axisPt});
    histos.add("Xi/hGenXi", "hGenXi", kTH2D, {axisRadius, axisPt});
    histos.add("AntiXi/hGenAntiXi", "hGenAntiXi", kTH2D, {axisRadius, axisPt});
    histos.add("Omega/hGenOmega", "hGenOmega", kTH2D, {axisRadius, axisPt});
    histos.add("AntiOmega/hGenAntiOmega", "hGenAntiOmega", kTH2D, {axisRadius, axisPt});
  }

  bool canDecay(const int pdgCode)
  {
    return std::find(mEnabledDecays.begin(), mEnabledDecays.end(), pdgCode) != mEnabledDecays.end();
  } 

  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    for (const auto& particle : mcParticles) {
      if (!canDecay(particle.pdgCode())) {
        continue;
      }
      
      o2::track::TrackParCov o2track;
      o2::upgrade::convertMCParticleToO2Track(particle, o2track, pdgDB);
      std::vector<o2::upgrade::OTFParticle> decayDaughters = decayer.decayParticle(pdgDB, o2track, particle.pdgCode());
      for (const auto& dau : decayDaughters) {
        o2::track::TrackParCov dauTrack;
        o2::upgrade::convertOTFParticleToO2Track(dau, dauTrack, pdgDB);
        if (canDecay(dau.pdgCode)) {
          std::vector<o2::upgrade::OTFParticle> cascadingDaughers = decayer.decayParticle(pdgDB, dauTrack, dau.pdgCode);
          for (const auto& daudau : cascadingDaughers) {
            decayDaughters.push_back(daudau);
          }
        }
      }

      if (decayDaughters.empty()) {
        LOG(error) << "Attempted to decay " << particle.pdgCode() << " but resulting vector of daugthers were empty";
        continue;
      }
      
      
      const float decayRadius2D = std::hypot(decayDaughters[0].vx, decayDaughters[0].vy);
      std::cout << particle.pdgCode() << ": " << decayRadius2D << std::endl;
      if (particle.pdgCode() == kK0Short) {
        histos.fill(HIST("K0S/hGenK0S"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kLambda0) {
        histos.fill(HIST("Lambda/hGenLambda"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kLambda0Bar) {
        histos.fill(HIST("AntiLambda/hGenAntiLambda"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kXiMinus) {
        histos.fill(HIST("Xi/hGenXi"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kXiPlusBar) {
        histos.fill(HIST("AntiXi/hGenAntiXi"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kOmegaMinus) {
        histos.fill(HIST("Omega/hGenOmega"), decayRadius2D, particle.pt());
      } else if (particle.pdgCode() == kOmegaPlusBar) {
        histos.fill(HIST("AntiOmega/hGenAntiOmega"), decayRadius2D, particle.pt());
      }
      
      for (size_t i = 0; i < decayDaughters.size(); ++i) {
        const o2::upgrade::OTFParticle& dau = decayDaughters[i];
        const bool isAlive = !canDecay(dau.pdgCode);
        tableMcParticles(particle.globalIndex(),
          isAlive, dau.pdgCode,
          dau.vx, dau.vy, dau.vz,
          dau.px, dau.py, dau.pz);
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
