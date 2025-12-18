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
  o2::upgrade::Decayer decayer;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
    {defaultParameters[0], kNumDecays, kNumParameters, particleNames, parameterNames},
    "Enable option for particle to be decayed: 0 - no, 1 - yes"};
    
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
      std::cout << particle.pdgCode() << " decayed into: ";
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
      std::vector<bool> isFinalState;
      for (size_t i = 0; i < decayDaughters.size(); ++i) {
        isFinalState.push_back(!canDecay(decayDaughters[i].pdgCode));
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
