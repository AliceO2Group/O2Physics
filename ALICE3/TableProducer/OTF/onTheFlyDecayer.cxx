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

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <ostream>
#include <span>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr int NumDecays = 7;
static constexpr int NumParameters = 1;
static constexpr int DefaultParameters[NumDecays][NumParameters]{{1}, {1}, {1}, {1}, {1}, {1}, {1}};
static const std::vector<std::string> ParameterNames{"enable"};
static const std::vector<std::string> ParticleNames{"K0s",
                                                    "Lambda",
                                                    "Anti-Lambda",
                                                    "Xi",
                                                    "Anti-Xi",
                                                    "Omega",
                                                    "Anti-Omega"};

static const std::vector<int> pdgCodes{PDG_t::kK0Short,
                                       PDG_t::kLambda0,
                                       PDG_t::kLambda0Bar,
                                       PDG_t::kXiMinus,
                                       PDG_t::kXiPlusBar,
                                       PDG_t::kOmegaMinus,
                                       PDG_t::kOmegaPlusBar};

struct OnTheFlyDecayer {
  Produces<aod::McPartWithDaus> tableMcParticlesWithDau;

  o2::upgrade::Decayer decayer;
  Service<o2::framework::O2DatabasePDG> pdgDB;
  std::map<int, std::vector<o2::upgrade::OTFParticle>> mDecayDaughters;

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
                                                {DefaultParameters[0], NumDecays, NumParameters, ParticleNames, ParameterNames},
                                                "Enable option for particle to be decayed: 0 - no, 1 - yes"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<int> mEnabledDecays;
  void init(o2::framework::InitContext&)
  {
    LOG(info) << " --- Initializing on-the-fly-decayer.";
    LOG(info) << " --- Using seed: " << seed;
    LOG(info) << " --- Using magnetic field: " << magneticField;
    decayer.setSeed(seed);
    decayer.setBField(magneticField);
    for (int i = 0; i < NumDecays; ++i) {
      if (enabledDecays->get(ParticleNames[i].c_str(), "enable")) {
        LOG(info) << " --- Decay enabled: " << pdgCodes[i];
        mEnabledDecays.push_back(pdgCodes[i]);
      }
    }

    auto hNaNBookkeeping = histos.add<TH1>("hNaNBookkeeping", "hNaNBookkeeping", kTH1D, {{2, -0.5, 1.5}});
    hNaNBookkeeping->GetXaxis()->SetBinLabel(1, "OK");
    hNaNBookkeeping->GetXaxis()->SetBinLabel(2, "NaN");
  }

  bool canDecay(const o2::upgrade::OTFParticle& particle)
  {
    if (particle.hasDaughters()) {
      return false;
    }

    return std::find(mEnabledDecays.begin(), mEnabledDecays.end(), particle.pdgCode()) != mEnabledDecays.end();
  }

  std::vector<o2::upgrade::OTFParticle> allParticles;
  void decayParticles(const int start, const int stop)
  {
    int ndau = 0;
    for (int i = start; i < stop; i++) {
      o2::upgrade::OTFParticle& particle = allParticles[i];
      if (particle.isFromMcParticles()) {
        particle.setIsPrimary(true);
        particle.setIsAlive(true);
      }

      if (!canDecay(particle)) {
        continue;
      }

      particle.setIsAlive(false);
      std::vector<o2::upgrade::OTFParticle> decayStack = decayer.decayParticle(pdgDB, particle);
      particle.setIndicesDaughter(allParticles.size(), allParticles.size() + (decayStack.size() - 1));
      for (o2::upgrade::OTFParticle daughter : decayStack) {
        daughter.setIndicesMother(i, i);
        daughter.setCollisionId(particle.collisionId());
        daughter.setIsAlive(true);
        daughter.setIsPrimary(false);
        allParticles.push_back(daughter);
        ndau++;
      }
    }

    if (start >= stop) {
      return;
    }

    decayParticles(stop, stop + ndau);
  }

  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    allParticles.clear();

    // First we copy the particles from the table into a vector that is extendable
    for (int index{0}; index < static_cast<int>(mcParticles.size()); ++index) {
      const auto& mcParticle = mcParticles.rawIteratorAt(index);
      allParticles.push_back(o2::upgrade::OTFParticle{mcParticle});
    }

    // Do all decays
    decayParticles(0, allParticles.size());

    // Fill output table
    for (int index{0}; index < static_cast<int>(allParticles.size()); ++index) {
      const auto& otfParticle = allParticles[index];

      if (otfParticle.hasNaN()) {
        histos.fill(HIST("hNaNBookkeeping"), 1);
      } else {
        histos.fill(HIST("hNaNBookkeeping"), 0);
      }

      // todo: status codes and vt
      tableMcParticlesWithDau(otfParticle.collisionId(), otfParticle.pdgCode(), otfParticle.statusCode(),
                              otfParticle.flags(), otfParticle.getMotherSpan(), otfParticle.getDaughters().data(), otfParticle.weight(),
                              otfParticle.px(), otfParticle.py(), otfParticle.pz(), otfParticle.e(),
                              otfParticle.vx(), otfParticle.vy(), otfParticle.vz(), otfParticle.vt(),
                              otfParticle.phi(), otfParticle.eta(), otfParticle.pt(), otfParticle.p(), otfParticle.y(),
                              otfParticle.isAlive(), otfParticle.isPrimary());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
