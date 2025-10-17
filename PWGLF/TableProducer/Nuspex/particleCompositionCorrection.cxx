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

/// \file particleCompositionCorrection.cxx
/// \brief Task to generate a table of dNdEta, pt and PID dependent weigths for MC particles to reflect the measured particle abundances
/// \author Mario Krüger <mario.kruger@cern.ch>

#include <PWGLF/DataModel/particleCompositionCorrectionTable.h>

#include "Tools/ML/model.h"

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TMCProcess.h>
#include <TPDGCode.h>

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::ml;

struct ParticleCompositionCorrection {
  // produces table of (dNdEta, pt, PID)-dependent weights for adjusting the particle composition in MC to better match the measured particle chemistry
  // applying these weights in the analysis allows mitigating on event-by-event basis the corresponding biases in efficiency and contamination of unidentified particles
  // weights are determined using measured pi,K,p,labda(as basis for sigma) spectra togehter with their counterpart from pythia simulations
  // they are interpolated in dNdEta and pt dimension using DNNs, which then provide particle fractions in data and MC and are stored as .onnx in the CCDB
  // weigths are assigned to the primary generated particle as well as its daughter particles (both from decay and material interactions)
  // weights are calculated only for particles within the configured kineamatic range and only for a distinct set of mother particles (see code)
  // assumes neutral particles require the same scaling as their charged counterparts (e.g. pi0 is scaled the same as pi+)
  // multi-strange baryons are assigned scaling factors of the sigma, which should be better than no scaling at all

  /*
   backlog:
   - store final neural networks in ccdb and implement accessing them
   - support collision systems beyond pp
   - add QA task illustrating improved mc/data matching of DCA distributions after scaling of secondaries
   - extend PCC weight table by columns with systematic variations (up/down)
   - investigate why particles can be their own mothers
   - investigate why histogram registry writing does not respect subfolder structure
   */

  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> etaCut{"etaCut", 0.8f, "eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.15f, "pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 10.f, "pt max cut"};
  Configurable<bool> skipSec{"skipSec", false, "dont calculate weights for secondaries"};
  Configurable<bool> enableQAHistos{"enableQAHistos", true, "enable qa histograms showing the effect of the PCC"};

  Configurable<std::string> modelNameData{"modelNameData", "ParticleFractions_pp_data.onnx", "Path to the .onnx file containing the particle fractions in data"};
  Configurable<std::string> modelNameMC{"modelNameMC", "ParticleFractions_pp_pythia.onnx", "Path to the .onnx file containing the particle fractions in MC"};

  OnnxModel particleFractionsData;
  OnnxModel particleFractionsMC;

  Produces<o2::aod::ParticleCompositionCorrection> pccTable;
  void init(InitContext const& cfgc);
  void process(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& particles);

  HistogramRegistry histos;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ParticleCompositionCorrection>(cfgc)};
}

void ParticleCompositionCorrection::init(InitContext const&)
{
  particleFractionsData.initModel(modelNameData.value, true);
  particleFractionsMC.initModel(modelNameMC.value, true);

  if (enableQAHistos) {
    std::vector<double> ptBinEdges = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                      0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                      2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                      6.0, 6.5, 7.0, 8.0, 9.0, 10.0};
    const AxisSpec ptAxis{ptBinEdges, "#it{p}_{T} (GeV/#it{c})", "pt"};

    histos.add("frac/data/pion", "", kTProfile, {ptAxis});
    histos.add("frac/data/kaon", "", kTProfile, {ptAxis});
    histos.add("frac/data/proton", "", kTProfile, {ptAxis});
    histos.add("frac/data/sigma", "", kTProfile, {ptAxis});
    histos.addClone("frac/data/", "frac/mc/"); // FIXME: bug in writer/registry moves one histogram out of subfolder...

    histos.add("weight/pion", "", kTProfile, {ptAxis});
    histos.add("weight/kaon", "", kTProfile, {ptAxis});
    histos.add("weight/proton", "", kTProfile, {ptAxis});
    histos.add("weight/sigma", "", kTProfile, {ptAxis});

    histos.add("weight/secDec", "", kTProfile, {ptAxis});
    histos.add("weight/secMat", "", kTProfile, {ptAxis});
  }
}

void ParticleCompositionCorrection::process(aod::McCollisions::iterator const&, aod::McParticles const& particles)
{
  // determine dNdEta of the collision
  float dNdEta = 0.f;
  for (const auto& particle : particles) {
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (std::abs(particle.eta()) >= 0.5) { // o2-linter: disable=magic-number (particle density at mid-rapidity)
      continue;
    }
    ++dNdEta;
  }

  // fill particle composition table with corresponding weights
  for (const auto& particle : particles) {
    float weight = 1.f;

    // calculate weights only inside the configured kinematic range
    if (std::abs(particle.eta()) < etaCut && particle.pt() > ptMinCut && particle.pt() < ptMaxCut) {
      auto refParticleID = particle.index();

      // find initial particle of secondaries from decays or interactions with material
      while (!particles.iteratorAt(refParticleID).producedByGenerator() && particles.iteratorAt(refParticleID).has_mothers()) {
        if (skipSec) {
          break;
        }
        const auto& mother = particle.mothers_first_as<aod::McParticles>();
        auto motherID = mother.globalIndex() - particles.offset();
        if (motherID == static_cast<uint64_t>(refParticleID)) {
          // FIXME: this needs to be investigated!
          // LOGP(error, "  - !!! particle is its own mother: {} (pdg: {}) status: {} ; process: {}", motherID, mother.pdgCode(), mother.getGenStatusCode(), mother.getProcess());
          break;
        }
        refParticleID = motherID;
      }

      if (const auto& refParticle = particles.iteratorAt(refParticleID); refParticle.producedByGenerator()) {

        auto absPDGCode = std::abs(refParticle.pdgCode());
        // translate abs PDG code to PID variable of neural networks (0: pion, 1: kaon, 2: proton, 3: sigma)
        static const std::map<int, float> mapPID = {
          {PDG_t::kPiPlus, 0.f},
          {PDG_t::kPi0, 0.f},
          {PDG_t::kKPlus, 1.f},
          {PDG_t::kK0Short, 1.f},
          {PDG_t::kK0Long, 1.f},
          {PDG_t::kProton, 2.f},
          {PDG_t::kNeutron, 2.f},
          {PDG_t::kSigmaPlus, 3.f},
          {PDG_t::kSigmaMinus, 3.f},
          {PDG_t::kLambda0, 3.f},
          {PDG_t::kSigma0, 3.f},
          {PDG_t::kXiMinus, 3.f},
          // TODO: potentially extend by xi0/eta/omega/rho/phi/Delta...
        };

        if (auto iterMapPID = mapPID.find(absPDGCode); iterMapPID != mapPID.end()) {
          // LOGP(info, "scaling a {} with status code {} from process {}", particle.pdgCode(), particle.getGenStatusCode(), particle.getProcess());
          float pt = refParticle.pt();

          // calculate particle fractions and corresponding weight for given reference particle
          std::vector<std::vector<float>> input = {{dNdEta}, {pt}, {iterMapPID->second}};
          float fracData = particleFractionsData.evalModel(input)[0];
          float fracMC = particleFractionsMC.evalModel(input)[0];
          if (fracMC) {
            weight = fracData / fracMC;
          }

          if (enableQAHistos && std::abs(particle.eta()) < 0.8) { // o2-linter: disable=magic-number (usual range of charged-partilce measurements)
            if (particle.isPhysicalPrimary()) {
              if (iterMapPID->first == PDG_t::kPiPlus) {
                histos.fill(HIST("frac/data/pion"), pt, fracData);
                histos.fill(HIST("frac/mc/pion"), pt, fracMC);
                histos.fill(HIST("weight/pion"), pt, weight);
              }
              if (iterMapPID->first == PDG_t::kKPlus) {
                histos.fill(HIST("frac/data/kaon"), pt, fracData);
                histos.fill(HIST("frac/mc/kaon"), pt, fracMC);
                histos.fill(HIST("weight/kaon"), pt, weight);
              }
              if (iterMapPID->first == PDG_t::kProton) {
                histos.fill(HIST("frac/data/proton"), pt, fracData);
                histos.fill(HIST("frac/mc/proton"), pt, fracMC);
                histos.fill(HIST("weight/proton"), pt, weight);
              }
              if (iterMapPID->first == PDG_t::kSigmaPlus || iterMapPID->first == PDG_t::kSigmaMinus) {
                histos.fill(HIST("frac/data/sigma"), pt, fracData);
                histos.fill(HIST("frac/mc/sigma"), pt, fracMC);
                histos.fill(HIST("weight/sigma"), pt, weight);
              }
            } else {
              auto pdgParticle = pdg->GetParticle(particle.pdgCode());
              if (pdgParticle && pdgParticle->Charge() != 0.) {
                // fill the decay daughter info only (usually there is max. one parent generation so no relevant double counting here)
                if (particle.getProcess() == TMCProcess::kPDecay) {
                  histos.fill(HIST("weight/secDec"), pt, weight);
                } else if (particle.getProcess() == TMCProcess::kPHInhelastic) {
                  histos.fill(HIST("weight/secMat"), pt, weight);
                }
              }
            }
          }
        }
      }
    }
    pccTable(weight, weight, weight); // TODO: put here systematic variations of the weights
  }
}
