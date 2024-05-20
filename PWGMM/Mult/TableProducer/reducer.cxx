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

#include "Selections.h"
#include "ReducedTables.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pwgmm::mult;
static constexpr float defparams[1][7] = {{0.142664, 1.40302, 2.61158, 1.20139, 5.24992, 2.16384, 0.871307}};
static constexpr double PI = o2::constants::math::PI;

struct Reducer {
  SliceCache cache;

  using BCs = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos, aod::MultsExtraMC>;
  using MCCollisionsNoHepMC = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::FT0Mults, aod::FDDMults, aod::ZDCMults, aod::PVMults>;
  using Particles = aod::McParticles;
  using Tracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  Produces<aod::StoredRBCs> rbcs;
  Produces<aod::StoredRMCCollisions> rmcc;
  Produces<aod::StoredRMCColLabels> rmcl;
  Produces<aod::StoredRCollisions> rc;
  Produces<aod::StoredRHepMCinfos> rhepmci;

  Configurable<float> reductionFactor{"reduction-factor", 1.e-4, "Reduction factor"};
  Configurable<LabeledArray<float>> params{"params", {defparams[0], 1, 7, {"pars"}, {"0bin", "l", "a", "n1", "p1", "n2", "p2"}}, "Multiplicity distribution parameterization"};

  Configurable<std::vector<double>> etaBins{"eta", {-1.5, -0.5, 0.5, 1.5}, "eta binning"};
  Configurable<std::vector<double>> phiBins{"phi", {0., PI / 2., PI, 3. * PI / 2., 2. * PI}, "phi binning"};

  Preslice<aod::Collisions> cperBC = aod::collision::bcId;
  Preslice<aod::McCollisions> mccperBC = aod::mccollision::bcId;
  Preslice<aod::McParticles> perMCc = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perC = aod::track::collisionId;

  std::random_device rd;
  std::mt19937 randomgen;
  std::uniform_real_distribution<float> dist;

  std::vector<int64_t> usedMCCs;
  std::vector<int64_t> usedLabels;
  std::vector<float> weights;

  std::vector<double> etabins;
  std::vector<double> phibins;

  std::vector<int> binned;

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
    randomgen.seed(rd());
    LOGP(debug, ">>> Starting with params: {}, {}, {}, {}, {}, {}, {}", params->get((int)0, int(0)), params->get((int)0, 1), params->get((int)0, 2),
         params->get((int)0, 3), params->get((int)0, 4), params->get((int)0, 5), params->get((int)0, 6));
    etabins = static_cast<std::vector<double>>(etaBins);
    phibins = static_cast<std::vector<double>>(phiBins);
    binned.resize((etabins.size() - 1) * (phibins.size() - 1));
  }

  int findBin(float eta, float phi)
  {
    // locate a bin in a linear array - the stride for phi is etabins.size() - 1
    if (std::abs(eta) < etabins[0]) { // underflow
      return -1;
    }
    auto e = std::lower_bound(etabins.begin(), etabins.end(), eta);
    if (e == etabins.end()) { // overflow
      return -1;
    }
    auto p = std::lower_bound(phibins.begin(), phibins.end(), phi);
    return std::distance(etabins.begin(), e) - 1 /* eta pos */ + (etabins.size() - 1) /* stride */ * (std::distance(phibins.begin(), p) - 1) /* phi pos */;
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
        auto threshold = reductionFactor / value;
        // if the fraction of events at this multiplicity is greater than the reduction factor, randomly discard it or add with the weight
        auto r = dist(randomgen);
        LOGP(debug, ">>> die roll: {}", r);
        if (r > threshold) {
          // dicard
          LOGP(debug, ">>> discarding");
        } else {
          // keep
          weights[i] = 1. / threshold;
          LOGP(debug, ">>> keeping with weight {}", weights[i]);
        }
      }
      ++i;
    }
  }

  expressions::Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                           ncheckbit(aod::track::trackCutFlag, trackSelectionITS);
  expressions::Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true);
  expressions::Filter fTrackSelectionDCA = ncheckbit(aod::track::trackCutFlag, trackSelectionDCA);

  void processFull(BCs::iterator const& bc,
                   MCCollisions const& mccollisions,
                   Collisions const& collisions,
                   Tracks const& tracks)
  {
    processGeneric(bc, mccollisions, collisions, tracks);
  }

  PROCESS_SWITCH(Reducer, processFull, "Full process with HepMC", false);

  void processLite(BCs::iterator const& bc,
                   MCCollisionsNoHepMC const& mccollisions,
                   Collisions const& collisions,
                   Tracks const& tracks)
  {
    processGeneric(bc, mccollisions, collisions, tracks);
  }

  PROCESS_SWITCH(Reducer, processLite, "Process without HepMC", true);

  template <typename TBCI, typename TMCC, typename TC, typename TT>
  void processGeneric(TBCI const& bc,
                      TMCC const& mccollisions,
                      TC const& collisions,
                      TT const& tracks)
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
      if constexpr (requires {mcc.processId(); mcc.pdf1(); }) {
        rhepmci(rmcc.lastIndex(), mcc.xsectGen(), mcc.ptHard(), mcc.nMPI(), mcc.processId(), mcc.id1(), mcc.id2(), mcc.pdfId1(), mcc.pdfId2(), mcc.x1(), mcc.x2(), mcc.scalePdf(), mcc.pdf1(), mcc.pdf2());
      }
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
      std::fill(binned.begin(), binned.end(), 0);
      auto stracks = tracks.sliceBy(perC, c.globalIndex());
      for (auto& track : stracks) {
        auto bin = findBin(track.eta(), track.phi());
        if (bin >= 0) {
          binned[bin] += 1;
        }
      }
      rc(bcId, c.posX(), c.posY(), c.posZ(), c.collisionTimeRes(), c.multFT0A(), c.multFT0C(), c.multFDDA(), c.multFDDC(), c.multZNA(), c.multZNC(), c.multNTracksPV(), c.multNTracksPVeta1(), c.multNTracksPVetaHalf(), binned);
      rmcl(usedLabels[std::distance(usedMCCs.begin(), pos)]);
      ++i;
    }
  }
};

struct ReducerTest {
  Configurable<int> maxMult{"maxMult", 300, "Max multiplicity bin"};
  HistogramRegistry r{"Common", {}};

  void init(InitContext const&)
  {
    r.add({"ReconstructedMultiplicity", " ; N_{trk}", {HistType::kTH1F, {{maxMult + 1, -0.5, maxMult + 0.5}}}});
    r.add({"GeneratedMultiplicity", " ; N_{particles}", {HistType::kTH1F, {{maxMult + 1, -0.5, maxMult + 0.5}}}});
    r.add({"ReconstructedMultiplicityUnweighted", " ; N_{trk}", {HistType::kTH1F, {{maxMult + 1, -0.5, maxMult + 0.5}}}});
    r.add({"GeneratedMultiplicityUnweighted", " ; N_{particles}", {HistType::kTH1F, {{maxMult + 1, -0.5, maxMult + 0.5}}}});
  }

  void process(aod::StoredRMCCollisions const& mccollisions, soa::Join<aod::StoredRCollisions, aod::StoredRMCColLabels> const& collisions)
  {
    for (auto& c : collisions) {
      r.fill(HIST("ReconstructedMultiplicity"), c.multNTracksPVeta1(), c.rmccollision_as<aod::StoredRMCCollisions>().weight());
      r.fill(HIST("ReconstructedMultiplicityUnweighted"), c.multNTracksPVeta1());
    }
    for (auto& c : mccollisions) {
      r.fill(HIST("GeneratedMultiplicity"), c.multMCNParticlesEta10(), c.weight());
      r.fill(HIST("GeneratedMultiplicityUnweighted"), c.multMCNParticlesEta10());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<Reducer>(cfgc), adaptAnalysisTask<ReducerTest>(cfgc)};
}
