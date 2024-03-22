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

#include <random>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

#include "ReducedTables.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
static constexpr float defparams[1][7] = {{0.142664, 1.40302, 2.61158, 1.20139, 5.24992, 2.16384, 0.871307}};

struct Reducer {
  using BCs = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos, aod::MultsExtraMC>;
  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::FT0Mults, aod::FDDMults, aod::ZDCMults, aod::PVMults>;
  using Particles = aod::McParticles;
  using Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>;

  Produces<aod::StoredRBCs> rbcs;
  Produces<aod::StoredRMCCollisions> rmcc;
  Produces<aod::StoredRMCColLabels> rmcl;
  Produces<aod::StoredRCollisions> rc;
  Produces<aod::StoredRHepMCinfos> rhepmci;

  Configurable<float> reductionFactor{"reduction-factor", 1.e-4, "Reduction factor"};
  Configurable<LabeledArray<float>> params{"params", {defparams[0], 1, 7, {"pars"}, {"0bin", "l", "a", "n1", "p1", "n2", "p2"}}, "Multiplicity distribution parameterization"};

  Preslice<aod::Collisions> cperBC = aod::collision::bcId;
  Preslice<aod::McCollisions> mccperBC = aod::mccollision::bcId;
  Preslice<aod::McParticles> perMCc = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perC = aod::track::collisionId;

  std::mt19937 randomgen;

  std::vector<int64_t> usedMCCs;
  std::vector<int64_t> usedLabels;
  std::vector<float> weights;

  double NormalizedDoubleNBD(double x)
  {
    // <n> = (p[2,4]^2)
    // k = (1+p[3,5]^2), k >= 1 -> smaller standar deviation, peak is pronounced
    // p[6] - normalization
    // alpha = p[1]^2 / ( 1 + p[1]^2), relative weight, 0 < alpha < 1

    return params->get((int)0, 6) *
           (params->get((int)0, 1) * params->get((int)0, 1) / (1. + params->get((int)0, 1) * params->get((int)0, 1)) * // alpha
                                                                                                                       // v1 +
              1. / (x * TMath::Beta(x, (1. + params->get((int)0, 3) * params->get((int)0, 3)))) * TMath::Power((params->get((int)0, 2) * params->get((int)0, 2)) / ((1. + params->get((int)0, 3) * params->get((int)0, 3)) + (params->get((int)0, 2) * params->get((int)0, 2))), x) * TMath::Power((1. + params->get((int)0, 3) * params->get((int)0, 3)) / ((1. + params->get((int)0, 3) * params->get((int)0, 3)) + (params->get((int)0, 2) * params->get((int)0, 2))), (1. + params->get((int)0, 3) * params->get((int)0, 3))) +
            1. / (1 + params->get((int)0, 1) * params->get((int)0, 1)) * // 1 - alpha
                                                                         // v2 );
              1. / (x * TMath::Beta(x, (1. + params->get((int)0, 5) * params->get((int)0, 5)))) * TMath::Power((params->get((int)0, 4) * params->get((int)0, 4)) / ((1. + params->get((int)0, 5) * params->get((int)0, 5)) + (params->get((int)0, 4) * params->get((int)0, 4))), x) * TMath::Power((1. + params->get((int)0, 5) * params->get((int)0, 5)) / ((1. + params->get((int)0, 5) * params->get((int)0, 5)) + (params->get((int)0, 4) * params->get((int)0, 4))), (1. + params->get((int)0, 5) * params->get((int)0, 5))));
  }

  void init(InitContext const&)
  {
    std::random_device randomdevice;
    randomgen.seed(randomdevice());
    LOGP(debug, ">>> Starting with params: {}, {}, {}, {}, {}, {}, {}", params->get((int)0, int(0)), params->get((int)0, 1), params->get((int)0, 2),
         params->get((int)0, 3), params->get((int)0, 4), params->get((int)0, 5), params->get((int)0, 6));
  }

  template <typename MCC>
  void checkSampling(MCC const& mccollisions)
  {
    weights.clear();
    weights.resize(mccollisions.size());
    std::fill(weights.begin(), weights.end(), (float)-1.);
    auto i = 0;
    LOGP(debug, ">>> {} MC collisions for BC", mccollisions.size());
    for (auto& mcc : mccollisions) {
      auto value = (mcc.multMCNParticlesEta10() == 0) ? params->get((int)0, (int)0) : NormalizedDoubleNBD((double)mcc.multMCNParticlesEta10());
      LOGP(debug, ">>> {} value for reduction (threshold {})", value, (float)reductionFactor);
      if (value < reductionFactor) {
        // if the fraction of events at this multiplicity is less than reduction factor, keep the event as is
        weights[i] = 1.f;
        LOGP(debug, ">>> keeping with weight {}", weights[i]);
      } else {
        // if the fraction of events at this multiplicity is greater than the reduction factor, randomly discard it or add with the weight
        std::uniform_real_distribution<float> d(0, value);
        auto r = d(randomgen);
        LOGP(debug, ">>> die roll: {}", r);
        if (r > reductionFactor) {
          // dicard
          LOGP(debug, ">>> discarding");
        } else {
          // keep
          weights[i] = value / reductionFactor;
          LOGP(debug, ">>> keeping with weight {}", weights[i]);
        }
      }
      ++i;
    }
  }

  void process(BCs::iterator const& bc,
               MCCollisions const& mccollisions,
               Collisions const& collisions)
  {
    usedMCCs.clear();
    usedLabels.clear();
    // if BC has no collisions, skip it
    if (mccollisions.size() == 0) {
      return;
    }
    checkSampling(mccollisions);
    // if all events are discarded, skip the BC
    if (std::all_of(weights.begin(), weights.end(), [](auto const& x) { return x < 0; })) {
      return;
    }
    // keep the BC
    rbcs(bc.runNumber());
    auto bcId = rbcs.lastIndex();
    auto i = 0;
    // check MC events
    for (auto& mcc : mccollisions) {
      // skip discarded events
      if (weights[i] < 0) {
        ++i;
        continue;
      }
      rmcc(bcId, weights[i], mcc.posX(), mcc.posY(), mcc.posZ(), mcc.impactParameter(), mcc.multMCFT0A(), mcc.multMCFT0C(), mcc.multMCNParticlesEta05(), mcc.multMCNParticlesEta10());
      rhepmci(rmcc.lastIndex(), mcc.xsectGen(), mcc.ptHard(), mcc.nMPI(), mcc.processId(), mcc.id1(), mcc.id2(), mcc.pdfId1(), mcc.pdfId2(), mcc.x1(), mcc.x2(), mcc.scalePdf(), mcc.pdf1(), mcc.pdf2());
      // remember used events so that the index relation can be preserved
      usedMCCs.push_back(mcc.globalIndex());
      usedLabels.push_back(rmcc.lastIndex());
      ++i;
    }
    i = 0;
    // check Reco events
    for (auto& c : collisions) {
      // discard fake events
      if (!c.has_mcCollision()) {
        ++i;
        continue;
      }
      // discard events for which MC event was discarded
      auto pos = std::find(usedMCCs.begin(), usedMCCs.end(), c.mcCollisionId());
      if (pos == usedMCCs.end()) {
        ++i;
        continue;
      }
      rc(bcId, c.posX(), c.posY(), c.posZ(), c.collisionTimeRes(), c.multFT0A(), c.multFT0C(), c.multFDDA(), c.multFDDC(), c.multZNA(), c.multZNC(), c.multNTracksPV(), c.multNTracksPVeta1(), c.multNTracksPVetaHalf());
      rmcl(usedLabels[std::distance(usedMCCs.begin(), pos)]);
      ++i;
    }
  }
};

struct ReducerTest {
  HistogramRegistry r{
    "Common",
    {
      {"ReconstructedMultiplicity", " ; N_{trk}", {HistType::kTH1F, {{301, -0.5, 300.5}}}},  //
      {"GeneratedMultiplicity", " ; N_{particles}", {HistType::kTH1F, {{301, -0.5, 300.5}}}} //
    }                                                                                        //
  };

  void process(aod::StoredRMCCollisions const& mccollisions, soa::Join<aod::StoredRCollisions, aod::StoredRMCColLabels> const& collisions)
  {
    for (auto& c : collisions) {
      r.fill(HIST("ReconstructedMultiplicity"), c.multNTracksPVeta1(), c.rmccollision_as<aod::StoredRMCCollisions>().weight());
    }
    for (auto& c : mccollisions) {
      r.fill(HIST("GeneratedMultiplicity"), c.multMCNParticlesEta10(), c.weight());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<Reducer>(cfgc), adaptAnalysisTask<ReducerTest>(cfgc)};
}
