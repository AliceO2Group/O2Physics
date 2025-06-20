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

/// \author Jasper Parkkila <jasper.parkkila@cern.ch>

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/detail/TypeTruncation.h"

#include <TFormula.h>

#include <experimental/type_traits>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

// #define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct Filter2Prong {
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 0, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgYMax, float, -1.0f, "Maximum candidate rapidity")
  O2_DEFINE_CONFIGURABLE(cfgImPart1Mass, float, o2::constants::physics::MassKPlus, "Daughter particle 1 mass in GeV")
  O2_DEFINE_CONFIGURABLE(cfgImPart2Mass, float, o2::constants::physics::MassKMinus, "Daughter particle 2 mass in GeV")
  O2_DEFINE_CONFIGURABLE(cfgImPart1PID, int, o2::track::PID::Kaon, "PID of daughter particle 1 (O2 PID ID)")
  O2_DEFINE_CONFIGURABLE(cfgImPart2PID, int, o2::track::PID::Kaon, "PID of daughter particle 2 (O2 PID ID)")
  O2_DEFINE_CONFIGURABLE(cfgImCutPt, float, 0.2f, "Minimal pT for candidates")
  O2_DEFINE_CONFIGURABLE(cfgImMinInvMass, float, 0.95f, "Minimum invariant mass (GeV)")
  O2_DEFINE_CONFIGURABLE(cfgImMaxInvMass, float, 1.07f, "Maximum invariant mass (GeV)")
  O2_DEFINE_CONFIGURABLE(cfgImSigmaFormula, std::string, "(z < 0.5 && x < 3.0) || (z >= 0.5 && x < 2.5 && y < 3.0)", "pT dependent daughter track sigma pass condition (x = TPC sigma, y = TOF sigma, z = pT)")

  HfHelper hfHelper;
  Produces<aod::CF2ProngTracks> output2ProngTracks;
  Produces<aod::CF2ProngTrackmls> output2ProngTrackmls;

  Produces<aod::CF2ProngMcParts> output2ProngMcParts;

  std::vector<float> mlvecd{};
  std::vector<float> mlvecdbar{};

  using HFCandidates = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using HFCandidatesML = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;

  template <class T>
  using HasMLProb = decltype(std::declval<T&>().mlProbD0());

  std::unique_ptr<TFormula> sigmaFormula;

  void init(InitContext&)
  {
    if (doprocessDataInvMass)
      sigmaFormula = std::make_unique<TFormula>("sigmaFormula", cfgImSigmaFormula.value.c_str());
  }

  template <class HFCandidatesType>
  void processDataT(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesType const& candidates)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision
    if (cfgVerbosity > 0 && candidates.size() > 0)
      LOGF(info, "Candidates for collision: %lu, cfcollisions: %lu, CFTracks: %lu", candidates.size(), cfcollisions.size(), cftracks.size());
    for (const auto& c : candidates) {
      int prongCFId[2] = {-1, -1};
      for (const auto& cftrack : cftracks) {
        if (c.prong0Id() == cftrack.trackId()) {
          prongCFId[0] = cftrack.globalIndex();
          break;
        }
      }
      for (const auto& cftrack : cftracks) {
        if (c.prong1Id() == cftrack.trackId()) {
          prongCFId[1] = cftrack.globalIndex();
          break;
        }
      }
      // look-up the collision id
      if ((c.hfflag() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) == 0)
        continue;
      if (cfgYMax >= 0.0f && std::abs(hfHelper.yD0(c)) > cfgYMax)
        continue;

      if (c.isSelD0() > 0) {
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           prongCFId[0], prongCFId[1], c.pt(), c.eta(), c.phi(), hfHelper.invMassD0ToPiK(c), aod::cf2prongtrack::D0ToPiK);
        if constexpr (std::experimental::is_detected<HasMLProb, typename HFCandidatesType::iterator>::value) {
          mlvecd.clear();
          for (const float val : c.mlProbD0())
            mlvecd.push_back(val);
          mlvecdbar.clear();
          for (const float val : c.mlProbD0bar())
            mlvecdbar.push_back(val);
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }

      if (c.isSelD0bar() > 0) {
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           prongCFId[0], prongCFId[1], c.pt(), c.eta(), c.phi(), hfHelper.invMassD0barToKPi(c), aod::cf2prongtrack::D0barToKPi);
        if constexpr (std::experimental::is_detected<HasMLProb, typename HFCandidatesType::iterator>::value) {
          mlvecd.clear();
          for (const float val : c.mlProbD0())
            mlvecd.push_back(val);
          mlvecdbar.clear();
          for (const float val : c.mlProbD0bar())
            mlvecdbar.push_back(val);
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }
    }
  }

  void processDataML(aod::Collisions::iterator const& col, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesML const& candidates)
  {
    processDataT(col, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processDataML, "Process data D0 candidates with ML", false);

  void processData(aod::Collisions::iterator const& col, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidates const& candidates)
  {
    processDataT(col, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processData, "Process data D0 candidates", true);

  using HFMCTrack = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  void processMC(aod::McCollisions::iterator const&, aod::CFMcParticleRefs const& cfmcparticles, [[maybe_unused]] HFMCTrack const& mcparticles)
  {
    // The main filter outputs the primary MC particles. Here we just resolve the daughter indices that are needed for the efficiency matching.
    for (const auto& r : cfmcparticles) {
      const auto& mcParticle = r.mcParticle_as<HFMCTrack>();
      if ((mcParticle.flagMcMatchGen() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) == 0 || mcParticle.daughtersIds().size() != 2) {
        output2ProngMcParts(-1, -1, aod::cf2prongtrack::Generic2Prong);
        continue;
      }
      int prongCFId[2] = {-1, -1};
      for (uint i = 0; i < 2; ++i) {
        for (const auto& cfmcpart : cfmcparticles) {
          if (mcParticle.daughtersIds()[i] == cfmcpart.mcParticleId()) {
            prongCFId[i] = cfmcpart.globalIndex();
            break;
          }
        }
      }
      output2ProngMcParts(prongCFId[0], prongCFId[1],
                          (mcParticle.pdgCode() >= 0 ? aod::cf2prongtrack::D0ToPiK : aod::cf2prongtrack::D0barToKPi) | ((mcParticle.originMcGen() & RecoDecay::OriginType::Prompt) ? aod::cf2prongmcpart::Prompt : 0));
    }
  }
  PROCESS_SWITCH(Filter2Prong, processMC, "Process MC 2-prong daughters", false);

  // Generic 2-prong invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  using PIDTrack = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  void processDataInvMass(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, PIDTrack const& tracks)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision
    for (const auto& cftrack1 : cftracks) {
      const auto& p1 = tracks.iteratorAt(cftrack1.trackId() - tracks.begin().globalIndex());
      if (p1.sign() != 1)
        continue;
      if (sigmaFormula->Eval(o2::aod::pidutils::tpcNSigma(cfgImPart1PID, p1), o2::aod::pidutils::tofNSigma(cfgImPart1PID, p1)) <= 0.0f)
        continue;
      for (const auto& cftrack2 : cftracks) {
        if (cftrack2.globalIndex() == cftrack1.globalIndex())
          continue;
        const auto& p2 = tracks.iteratorAt(cftrack2.trackId() - tracks.begin().globalIndex());
        if (p2.sign() != -1)
          continue;
        if (sigmaFormula->Eval(o2::aod::pidutils::tpcNSigma(cfgImPart2PID, p2), o2::aod::pidutils::tofNSigma(cfgImPart2PID, p2)) <= 0.0f)
          continue;
        ROOT::Math::PtEtaPhiMVector vec1(p1.pt(), p1.eta(), p1.phi(), cfgImPart1Mass);
        ROOT::Math::PtEtaPhiMVector vec2(p2.pt(), p2.eta(), p2.phi(), cfgImPart2Mass);
        ROOT::Math::PtEtaPhiMVector s = vec1 + vec2;
        if (s.pt() < cfgImCutPt || s.M() < cfgImMinInvMass || s.M() > cfgImMaxInvMass)
          continue;

        float phi = RecoDecay::constrainAngle(s.Phi(), 0.0f);
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           cftrack1.globalIndex(), cftrack2.globalIndex(), s.pt(), s.eta(), phi, s.M(), aod::cf2prongtrack::Generic2Prong);
      }
    }
  }
  PROCESS_SWITCH(Filter2Prong, processDataInvMass, "Process data generic 2-prong candidates with invariant mass method", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Filter2Prong>(cfgc, TaskName{"filter-2prong"})};
}
