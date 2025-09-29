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

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2;

// Output-table definition
namespace o2::aod
{
namespace trd::col
{
DECLARE_SOA_COLUMN(Q00, q00, float);
DECLARE_SOA_COLUMN(Q01, q01, float);
DECLARE_SOA_COLUMN(Q02, q02, float);
DECLARE_SOA_COLUMN(Q03, q03, float);
DECLARE_SOA_COLUMN(Q04, q04, float);
DECLARE_SOA_COLUMN(Q05, q05, float);
DECLARE_SOA_COLUMN(Q10, q10, float);
DECLARE_SOA_COLUMN(Q11, q11, float);
DECLARE_SOA_COLUMN(Q12, q12, float);
DECLARE_SOA_COLUMN(Q13, q13, float);
DECLARE_SOA_COLUMN(Q14, q14, float);
DECLARE_SOA_COLUMN(Q15, q15, float);
DECLARE_SOA_COLUMN(Q20, q20, float);
DECLARE_SOA_COLUMN(Q21, q21, float);
DECLARE_SOA_COLUMN(Q22, q22, float);
DECLARE_SOA_COLUMN(Q23, q23, float);
DECLARE_SOA_COLUMN(Q24, q24, float);
DECLARE_SOA_COLUMN(Q25, q25, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
} // namespace trd::col
DECLARE_SOA_TABLE(TRDPID, "AOD", "TRDPID", o2::soa::Index<>,
                  trd::col::Q00, trd::col::Q01, trd::col::Q02, trd::col::Q03, trd::col::Q04, trd::col::Q05,
                  trd::col::Q10, trd::col::Q11, trd::col::Q12, trd::col::Q13, trd::col::Q14, trd::col::Q15,
                  trd::col::Q20, trd::col::Q21, trd::col::Q22, trd::col::Q23, trd::col::Q24, trd::col::Q25,
                  trd::col::Pt);
} // namespace o2::aod

struct TRDPIDStudy {
  Produces<aod::TRDPID> pidTable;
  HistogramRegistry mRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Config
  Configurable<bool> mFilterNeighbors{"filterNeighbors", true, "Filter tracklets with neighbors"};
  Configurable<bool> mFilterCrossings{"filterCrossings", true, "Filter tracklets crossing pads"};
  Configurable<size_t> mMinTracklets{"minTracklets", 4, "Minimum number of tracklets"};
  Configurable<float> mMinPt{"minPt", 0.1, "Minimum track-pt required"};
  Configurable<float> mMaxPt{"maxPt", 20, "Maximum track-pt allowed"};

  using Tracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;

  void init(InitContext& /*ctx*/)
  {
    mRegistry.add("Q0", "Q0", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("Q1", "Q1", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("Q2", "Q2", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("QTot", "QTot", HistType::kTH1F, {{300, 0, 300}});
    mRegistry.add("Q0Cor", "Q0Cor", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("Q1Cor", "Q1Cor", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("Q2Cor", "Q2Cor", HistType::kTH1F, {{127, 0, 127}});
    mRegistry.add("QTotCor", "QTotCor", HistType::kTH1F, {{300, 0, 300}});
  }

  void process(Tracks const& tracks, aod::TRDsExtra const& trdExtra)
  {
    for (const auto& trd : trdExtra) {
      const auto& track = tracks.rawIteratorAt(trd.trackId());
      std::bitset<6> good;
      if (!filterTrack(track, trd, good)) {
        continue;
      }

      const auto& q0s = trd.trdQ0s();
      const auto& q1s = trd.trdQ1s();
      const auto& q2s = trd.trdQ2s();
      const auto& q0sCor = trd.trdQ0sCorrected();
      const auto& q1sCor = trd.trdQ1sCorrected();
      const auto& q2sCor = trd.trdQ2sCorrected();

      for (int i{0}; i < 6; ++i) {
        if (!good[i]) {
          continue;
        }

        mRegistry.fill(HIST("Q0"), q0s[i]);
        mRegistry.fill(HIST("Q1"), q1s[i]);
        mRegistry.fill(HIST("Q2"), q2s[i]);
        mRegistry.fill(HIST("QTot"), q0s[i] + q1s[i] + q2s[i]);
        mRegistry.fill(HIST("Q0Cor"), q0sCor[i]);
        mRegistry.fill(HIST("Q1Cor"), q1sCor[i]);
        mRegistry.fill(HIST("Q2Cor"), q2sCor[i]);
        mRegistry.fill(HIST("QTotCor"), q0sCor[i] + q1sCor[i] + q2sCor[i]);
      }

      pidTable(
        q0sCor[0], q0sCor[1], q0sCor[2], q0sCor[3], q0sCor[4], q0sCor[5],
        q1sCor[0], q1sCor[1], q1sCor[2], q1sCor[3], q1sCor[4], q1sCor[5],
        q2sCor[0], q2sCor[1], q2sCor[2], q2sCor[3], q2sCor[4], q2sCor[5],
        track.pt());
    }
  }

  template <typename Track, typename TRD>
  bool filterTrack(Track const& trk, TRD const& trd, std::bitset<6>& good)
  {
    if (trk.trdNTracklets() < mMinTracklets) {
      return false;
    }
    if (trk.trdHasNeighbor()) {
      return false;
    }
    if (trk.trdHasCrossing()) {
      return false;
    }
    const auto& q0s = trd.trdQ0s();
    const auto& q1s = trd.trdQ1s();
    const auto& q2s = trd.trdQ2s();
    for (int i{0}; i < 6; ++i) {
      if ((trk.trdPattern() & (1 << i)) == 0) {
        continue;
      }
      if (q2s[i] >= 62 || q2s[i] < 6) {
        continue;
      }
      if (q1s[i] >= 127 || q1s[i] < 6) {
        continue;
      }
      if (q0s[i] >= 127 || q0s[i] < 6) {
        continue;
      }
      good.set(i);
    }

    return good.count() >= mMinTracklets;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TRDPIDStudy>(cfgc),
  };
}
