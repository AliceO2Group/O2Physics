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

#include <gsl/span>
#include <map>
#include <span>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr int NumDecays = 7;
static constexpr int NumParameters = 1;
static constexpr int DefaultParameters[NumDecays][NumParameters]{{1}, {1}, {1}, {1}, {1}, {1}, {1}};
static constexpr float Tolerance = 1e-7f;
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
  std::map<int, std::vector<o2::upgrade::OTFParticle>> mDecayDaughters;

  Configurable<int> seed{"seed", 0, "Set seed for particle decayer"};
  Configurable<float> magneticField{"magneticField", 20., "Magnetic field (kG)"};
  Configurable<LabeledArray<int>> enabledDecays{"enabledDecays",
                                                {DefaultParameters[0], NumDecays, NumParameters, particleNames, parameterNames},
                                                "Enable option for particle to be decayed: 0 - no, 1 - yes"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  struct McParticleAlice3 {
    McParticleAlice3() = default;
    ~McParticleAlice3() = default;
    McParticleAlice3(const McParticleAlice3& src) = default;
    McParticleAlice3(int collisionId,
                     int pdgCode,
                     int statusCode,
                     int flags,
                     int mother0,
                     int mother1,
                     int daughter0,
                     int daughter1,
                     float weight,
                     float px, float py, float pz, float e,
                     float vx, float vy, float vz, float vt,
                     float phi, float eta, float pt, float p, float y,
                     bool isAlive, bool isPrimary) : collisionId(collisionId),
                                                     pdgCode(pdgCode),
                                                     statusCode(statusCode),
                                                     flags(flags),
                                                     mothersIds{mother0, mother1},
                                                     daughtersIdSlice{daughter0, daughter1},
                                                     weight(weight),
                                                     px(px),
                                                     py(py),
                                                     pz(pz),
                                                     e(e),
                                                     vx(vx),
                                                     vy(vy),
                                                     vz(vz),
                                                     vt(vt),
                                                     phi(phi),
                                                     eta(eta),
                                                     pt(pt),
                                                     p(p),
                                                     y(y),
                                                     isAlive(isAlive),
                                                     isPrimary(isPrimary) {}
    int collisionId;
    int pdgCode;
    int statusCode;
    int flags;
    int mothersIds[2];
    int daughtersIdSlice[2];
    float weight;
    float px, py, pz, e;
    float vx, vy, vz, vt;
    float phi, eta, pt, p, y;
    bool isAlive;
    bool isPrimary;
  };

  std::vector<int> mEnabledDecays;
  void init(o2::framework::InitContext&)
  {
    LOG(info) << "Initializing on-the-fly-decayer.";
    decayer.setSeed(seed);
    decayer.setBField(magneticField);
    for (int i = 0; i < NumDecays; ++i) {
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

    histos.add("GeneratedElectron/hGenEl", "hGenEl", kTH1D, {axisPt});
    histos.add("GeneratedMuon/hGenMu", "hGenMu", kTH1D, {axisPt});
    histos.add("GeneratedPion/hGenPi", "hGenPi", kTH1D, {axisPt});
    histos.add("GeneratedKaon/hGenKa", "hGenKa", kTH1D, {axisPt});
    histos.add("GeneratedProton/hGenPr", "hGenPr", kTH1D, {axisPt});
  }

  bool canDecay(const int pdgCode)
  {
    return std::find(mEnabledDecays.begin(), mEnabledDecays.end(), pdgCode) != mEnabledDecays.end();
  }

  std::vector<McParticleAlice3> mcParticlesAlice3;
  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    mDecayDaughters.clear();
    mcParticlesAlice3.clear();
    u_int64_t nStoredDaughters = 0;
    for (int index{0}; index < static_cast<int>(mcParticles.size()); ++index) {
      const auto& particle = mcParticles.iteratorAt(index);
      std::vector<o2::upgrade::OTFParticle> decayDaughters;
      static constexpr int MaxNestedDecays = 10;
      int nDecays = 0;
      if (canDecay(particle.pdgCode())) {
        o2::track::TrackParCov o2track;
        o2::upgrade::convertMCParticleToO2Track(particle, o2track, pdgDB);
        decayDaughters = decayer.decayParticle(pdgDB, o2track, particle.pdgCode());
        for (size_t idau{0}; idau < decayDaughters.size(); ++idau) {
          o2::upgrade::OTFParticle dau = decayDaughters[idau];
          o2::track::TrackParCov dauTrack;
          o2::upgrade::convertOTFParticleToO2Track(dau, dauTrack, pdgDB);
          if (canDecay(dau.pdgCode())) {
            dau.setIsAlive(false);
            std::vector<o2::upgrade::OTFParticle> cascadingDaughers = decayer.decayParticle(pdgDB, dauTrack, dau.pdgCode());
            for (size_t idaudau{0}; idaudau < cascadingDaughers.size(); ++idaudau) {
              o2::upgrade::OTFParticle daudau = cascadingDaughers[idaudau];
              decayDaughters.push_back(daudau);
              if (MaxNestedDecays < ++nDecays) {
                LOG(error) << "Seemingly stuck trying to perpetually decay products from pdg: " << particle.pdgCode();
              }
            }
          } else {
            dau.setIsAlive(true);
          }
        }

        if (decayDaughters.empty()) {
          LOG(error) << "Attempted to decay " << particle.pdgCode() << " but resulting vector of daugthers were empty";
          continue;
        }

        switch (particle.pdgCode()) {
          case kK0Short:
            histos.fill(HIST("K0S/hGenK0S"), particle.pt());
            break;

          case kLambda0:
            histos.fill(HIST("Lambda/hGenLambda"), particle.pt());
            break;

          case kLambda0Bar:
            histos.fill(HIST("AntiLambda/hGenAntiLambda"), particle.pt());
            break;

          case kXiMinus:
            histos.fill(HIST("Xi/hGenXi"), particle.pt());
            break;

          case kXiPlusBar:
            histos.fill(HIST("AntiXi/hGenAntiXi"), particle.pt());
            break;

          case kOmegaMinus:
            histos.fill(HIST("Omega/hGenOmega"), particle.pt());
            break;

          case kOmegaPlusBar:
            histos.fill(HIST("AntiOmega/hGenAntiOmega"), particle.pt());
            break;

          default:
            break;
        }
      }

      int daughtersIdSlice[2];
      if (canDecay(particle.pdgCode())) {
        daughtersIdSlice[0] = static_cast<int>(mcParticles.size() + nStoredDaughters);
        daughtersIdSlice[1] = static_cast<int>(mcParticles.size() + nStoredDaughters + decayDaughters.size());
      } else {
        daughtersIdSlice[0] = static_cast<int>(particle.daughtersIds()[0]);
        daughtersIdSlice[1] = static_cast<int>(particle.daughtersIds()[1]);
      }

      mDecayDaughters.emplace(index, decayDaughters);
      nStoredDaughters += decayDaughters.size();

      const float phi = o2::constants::math::PI + std::atan2(-1.0f * particle.py(), -1.0f * particle.px());
      float eta; // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
      const float pt = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py());
      const float p = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py() + particle.pz() * particle.pz());
      float y; // As https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

      if ((p - particle.pz()) < Tolerance) {
        eta = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        eta = 0.5f * std::log((p + particle.pz()) / (p - particle.pz()));
      }

      if ((particle.e() - particle.pz()) < Tolerance) {
        y = (particle.pz() < 0.0f) ? -100.0f : 100.0f;
      } else {
        y = 0.5f * std::log((particle.e() + particle.pz()) / (particle.e() - particle.pz()));
      }

      // TODO: Particle status code
      // TODO: Expression columns
      auto mothers = particle.mothersIds();
      int mother0 = mothers.size() > 0 ? mothers[0] : -1;
      int mother1 = mothers.size() > 1 ? mothers[1] : mother0;
      mcParticlesAlice3.push_back(McParticleAlice3{particle.mcCollisionId(), particle.pdgCode(), particle.statusCode(),
                                                   particle.flags(), mother0, mother1,
                                                   daughtersIdSlice[0], daughtersIdSlice[1], particle.weight(),
                                                   particle.px(), particle.py(), particle.pz(), particle.e(),
                                                   particle.vx(), particle.vy(), particle.vz(), particle.vt(),
                                                   phi, eta, pt, p, y, !canDecay(particle.pdgCode()), true});
    }

    int daughtersIdSlice[2] = {-1, -1};
    for (const auto& [index, decayDaughters] : mDecayDaughters) {
      for (const auto& dau : decayDaughters) {
        if (index >= mcParticles.size()) {
          LOG(error) << "--- Index " << index << " out of bounds for mcParticles table of size " << mcParticles.size() << std::endl;
          continue;
        }

        const float phi = o2::constants::math::PI + std::atan2(-1.0f * dau.py(), -1.0f * dau.px());
        float eta; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1922
        const float pt = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py());
        const float p = std::sqrt(dau.px() * dau.px() + dau.py() * dau.py() + dau.pz() * dau.pz());
        float y; // Conditional as https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L1943

        if ((p - dau.pz()) < Tolerance) {
          eta = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          eta = 0.5f * std::log((p + dau.pz()) / (p - dau.pz()));
        }

        if ((dau.e() - dau.pz()) < Tolerance) {
          y = (dau.pz() < 0.0f) ? -100.0f : 100.0f;
        } else {
          y = 0.5f * std::log((dau.e() + dau.pz()) / (dau.e() - dau.pz()));
        }

        switch (dau.pdgCode()) {
          case kElectron:
            histos.fill(HIST("GeneratedElectron/hGenEl"), pt);
            break;

          case kMuonMinus:
            histos.fill(HIST("GeneratedMuon/hGenMu"), pt);
            break;

          case kPiPlus:
            histos.fill(HIST("GeneratedPion/hGenPi"), pt);
            break;

          case kKPlus:
            histos.fill(HIST("GeneratedKaon/hGenKa"), pt);
            break;

          case kProton:
            histos.fill(HIST("GeneratedProton/hGenPr"), pt);
            break;

          default:
            break;
        }

        // TODO: Particle status code
        // TODO: Expression columns
        // TODO: vt
        auto mother = mcParticles.iteratorAt(index);
        mcParticlesAlice3.push_back(McParticleAlice3{mother.mcCollisionId(), dau.pdgCode(), 1,
                                                     -1, index, index, daughtersIdSlice[0], daughtersIdSlice[1], mother.weight(),
                                                     dau.px(), dau.py(), dau.pz(), dau.e(),
                                                     dau.vx(), dau.vy(), dau.vz(), mother.vt(),
                                                     phi, eta, pt, p, y, dau.isAlive(), false});
      }
    }

    for (const auto& particle : mcParticlesAlice3) {
      std::span<const int> motherSpan(particle.mothersIds, 2);

      tableMcParticlesWithDau(particle.collisionId, particle.pdgCode, particle.statusCode,
                              particle.flags, motherSpan, particle.daughtersIdSlice, particle.weight,
                              particle.px, particle.py, particle.pz, particle.e,
                              particle.vx, particle.vy, particle.vz, particle.vt,
                              particle.phi, particle.eta, particle.pt, particle.p, particle.y,
                              particle.isAlive, particle.isPrimary);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyDecayer>(cfgc)};
}
