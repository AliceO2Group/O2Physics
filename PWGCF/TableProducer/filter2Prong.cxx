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
#include <experimental/type_traits>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "MathUtils/detail/TypeTruncation.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

// #define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct Filter2Prong {
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 0, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgYMax, float, -1.0f, "Maximum candidate rapidity")

  HfHelper hfHelper;
  Produces<aod::CF2ProngTracks> output2ProngTracks;
  Produces<aod::CF2ProngTrackmls> output2ProngTrackmls;

  std::vector<float> mlvecd{};
  std::vector<float> mlvecdbar{};

  using HFCandidates = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using HFCandidatesML = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;

  template <class T>
  using HasMLProb = decltype(std::declval<T&>().mlProbD0());

  template <class HFCandidatesType>
  void processDataT(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesType const& candidates)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision
    if (cfgVerbosity > 0 && candidates.size() > 0)
      LOGF(info, "Candidates for collision: %lu, cfcollisions: %lu, CFTracks: %lu", candidates.size(), cfcollisions.size(), cftracks.size());
    for (auto& c : candidates) {
      int prongCFId[2] = {-1, -1};
      for (auto& cftrack : cftracks) {
        if (c.prong0Id() == cftrack.trackId()) {
          prongCFId[0] = cftrack.globalIndex();
          break;
        }
      }
      for (auto& cftrack : cftracks) {
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
          for (float val : c.mlProbD0()) {
            mlvecd.push_back(val);
          }
          mlvecdbar.clear();
          for (float val : c.mlProbD0bar()) {
            mlvecdbar.push_back(val);
          }
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }

      if (c.isSelD0bar() > 0) {
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           prongCFId[0], prongCFId[1], c.pt(), c.eta(), c.phi(), hfHelper.invMassD0barToKPi(c), aod::cf2prongtrack::D0barToKPi);
        if constexpr (std::experimental::is_detected<HasMLProb, typename HFCandidatesType::iterator>::value) {
          mlvecd.clear();
          for (float val : c.mlProbD0()) {
            mlvecd.push_back(val);
          }
          mlvecdbar.clear();
          for (float val : c.mlProbD0bar()) {
            mlvecdbar.push_back(val);
          }
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }
    }
  }

  void processDataML(aod::Collisions::iterator const& cols, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesML const& candidates)
  {
    processDataT(cols, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processDataML, "Process data D0 candidates with ML", false);

  void processData(aod::Collisions::iterator const& cols, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidates const& candidates)
  {
    processDataT(cols, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processData, "Process data D0 candidates", true);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Filter2Prong>(cfgc, TaskName{"filter-2prong"})};
}
