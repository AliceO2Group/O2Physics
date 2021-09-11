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
/// \file   pidTOFFull.cxx
/// \author Nicolo' Jacazio
/// \brief  Task to produce PID tables for TOF split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

struct trackTime {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  Configurable<bool> cut{"cut", false, "Apply track cuts"};
  TrackSelection globalTracks; // Track with cut for primaries

  void init(o2::framework::InitContext&)
  {
    globalTracks = getGlobalTrackSelection();

    const AxisSpec pAxis{1000, 0, 10, "#it{p} (GeV/#it{c})"};
    const AxisSpec diffAxis{1000, -100, 100, "diff (ps)"};
    const AxisSpec diffrelAxis{10000, -10, 10, "diffrel"};
    const AxisSpec signalAxis{1000, -10000, 10000, "t-t_{0}-t_{exp}"};
    const AxisSpec tofAxis{1000, 0, 100000, "tof"};
    const AxisSpec trkAxis{1000, 0, 100000, "trk"};

    histos.add("vs/all", "all", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/el", "el", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/mu", "mu", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/pi", "pi", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/ka", "ka", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/pr", "pr", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/de", "de", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/tr", "tr", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/he", "he", kTH2F, {tofAxis, trkAxis});
    histos.add("vs/al", "al", kTH2F, {tofAxis, trkAxis});

    histos.add("diff/all", "all", kTH2F, {pAxis, diffAxis});
    histos.add("diff/el", "el", kTH2F, {pAxis, diffAxis});
    histos.add("diff/mu", "mu", kTH2F, {pAxis, diffAxis});
    histos.add("diff/pi", "pi", kTH2F, {pAxis, diffAxis});
    histos.add("diff/ka", "ka", kTH2F, {pAxis, diffAxis});
    histos.add("diff/pr", "pr", kTH2F, {pAxis, diffAxis});
    histos.add("diff/de", "de", kTH2F, {pAxis, diffAxis});
    histos.add("diff/tr", "tr", kTH2F, {pAxis, diffAxis});
    histos.add("diff/he", "he", kTH2F, {pAxis, diffAxis});
    histos.add("diff/al", "al", kTH2F, {pAxis, diffAxis});

    histos.add("difftof/all", "all", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/el", "el", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/mu", "mu", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/pi", "pi", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/ka", "ka", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/pr", "pr", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/de", "de", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/tr", "tr", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/he", "he", kTH2F, {tofAxis, diffAxis});
    histos.add("difftof/al", "al", kTH2F, {tofAxis, diffAxis});

    histos.add("difftofrel/all", "all", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/el", "el", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/mu", "mu", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/pi", "pi", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/ka", "ka", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/pr", "pr", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/de", "de", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/tr", "tr", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/he", "he", kTH2F, {tofAxis, diffrelAxis});
    histos.add("difftofrel/al", "al", kTH2F, {tofAxis, diffrelAxis});

    histos.add("difftrk/all", "all", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/el", "el", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/mu", "mu", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/pi", "pi", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/ka", "ka", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/pr", "pr", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/de", "de", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/tr", "tr", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/he", "he", kTH2F, {trkAxis, diffAxis});
    histos.add("difftrk/al", "al", kTH2F, {trkAxis, diffAxis});

    histos.add("difftrkrel/all", "all", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/el", "el", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/mu", "mu", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/pi", "pi", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/ka", "ka", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/pr", "pr", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/de", "de", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/tr", "tr", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/he", "he", kTH2F, {trkAxis, diffrelAxis});
    histos.add("difftrkrel/al", "al", kTH2F, {trkAxis, diffrelAxis});

    histos.add("tof/el", "el", kTH2F, {pAxis, signalAxis});
    histos.add("tof/mu", "mu", kTH2F, {pAxis, signalAxis});
    histos.add("tof/pi", "pi", kTH2F, {pAxis, signalAxis});
    histos.add("tof/ka", "ka", kTH2F, {pAxis, signalAxis});
    histos.add("tof/pr", "pr", kTH2F, {pAxis, signalAxis});
    histos.add("tof/de", "de", kTH2F, {pAxis, signalAxis});
    histos.add("tof/tr", "tr", kTH2F, {pAxis, signalAxis});
    histos.add("tof/he", "he", kTH2F, {pAxis, signalAxis});
    histos.add("tof/al", "al", kTH2F, {pAxis, signalAxis});

    histos.add("trk/el", "el", kTH2F, {pAxis, signalAxis});
    histos.add("trk/mu", "mu", kTH2F, {pAxis, signalAxis});
    histos.add("trk/pi", "pi", kTH2F, {pAxis, signalAxis});
    histos.add("trk/ka", "ka", kTH2F, {pAxis, signalAxis});
    histos.add("trk/pr", "pr", kTH2F, {pAxis, signalAxis});
    histos.add("trk/de", "de", kTH2F, {pAxis, signalAxis});
    histos.add("trk/tr", "tr", kTH2F, {pAxis, signalAxis});
    histos.add("trk/he", "he", kTH2F, {pAxis, signalAxis});
    histos.add("trk/al", "al", kTH2F, {pAxis, signalAxis});
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended>;
  using Trks2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov,
                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                          aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                          aod::TrackSelection>;
  using Coll = aod::Collisions;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = tof::ExpTimes<Trks::iterator, pid>;
  void process(aod::Collision const& collision, Trks const& tracks)
  {

    constexpr auto responseEl = ResponseImplementation<PID::Electron>();
    constexpr auto responseMu = ResponseImplementation<PID::Muon>();
    constexpr auto responsePi = ResponseImplementation<PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<PID::Proton>();
    constexpr auto responseDe = ResponseImplementation<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementation<PID::Triton>();
    constexpr auto responseHe = ResponseImplementation<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementation<PID::Alpha>();

    for (auto& t : tracks) {
      if (cut && !globalTracks.IsSelected(t)) {
        continue;
      }

      float diff = 9999.f;
      float ttrk = 9999.f;
      float texp = 9999.f;
      if (!t.hasTOF()) {
        continue;
      }

#define CasePID(resp, part)                                                   \
  texp = resp.GetExpectedSignal(t);                                           \
  ttrk = t.trackTime() * 1e+3 + texp;                                         \
  diff = t.tofSignal() - ttrk;                                                \
  histos.fill(HIST("vs/" part), t.tofSignal(), ttrk);                         \
  histos.fill(HIST("diff/" part), t.p(), diff);                               \
  histos.fill(HIST("difftof/" part), t.tofSignal(), diff);                    \
  histos.fill(HIST("difftrk/" part), ttrk, diff);                             \
  histos.fill(HIST("difftofrel/" part), t.tofSignal(), diff / t.tofSignal()); \
  histos.fill(HIST("difftrkrel/" part), ttrk, diff / ttrk);

      const float evtime = t.collision().collisionTime() * 1000.f;
      const float Lc = t.length() / 0.029979246f;
      switch (t.pidForTracking()) {
        case 0:
          CasePID(responseEl, "el");
          break;
        case 1:
          CasePID(responseMu, "mu");
          break;
        case 2:
          CasePID(responsePi, "pi");
          break;
        case 3:
          CasePID(responseKa, "ka");
          break;
        case 4:
          CasePID(responsePr, "pr");
          break;
        case 5:
          CasePID(responseDe, "de");
          break;
        case 6:
          CasePID(responseTr, "tr");
          break;
        case 7:
          CasePID(responseHe, "he");
          break;
        case 8:
          CasePID(responseAl, "al");
          break;
        default:
          break;
      }
#undef CasePID

      histos.fill(HIST("vs/all"), t.tofSignal(), ttrk);
      histos.fill(HIST("diff/all"), t.p(), diff);
      histos.fill(HIST("difftof/all"), t.tofSignal(), diff);
      histos.fill(HIST("difftrk/all"), ttrk, diff);
      histos.fill(HIST("difftofrel/all"), t.tofSignal(), diff / t.tofSignal());
      histos.fill(HIST("difftrkrel/all"), ttrk, diff / ttrk);

      histos.fill(HIST("tof/el"), t.p(), t.tofSignal() - evtime - responseEl.GetExpectedSignal(t));
      histos.fill(HIST("trk/el"), t.p(), ttrk - evtime - responseEl.GetExpectedSignal(t));

      histos.fill(HIST("tof/mu"), t.p(), t.tofSignal() - evtime - responseMu.GetExpectedSignal(t));
      histos.fill(HIST("trk/mu"), t.p(), ttrk - evtime - responseMu.GetExpectedSignal(t));
      histos.fill(HIST("tof/pi"), t.p(), t.tofSignal() - evtime - responsePi.GetExpectedSignal(t));
      histos.fill(HIST("trk/pi"), t.p(), ttrk - evtime - responsePi.GetExpectedSignal(t));
      histos.fill(HIST("tof/ka"), t.p(), t.tofSignal() - evtime - responseKa.GetExpectedSignal(t));
      histos.fill(HIST("trk/ka"), t.p(), ttrk - evtime - responseKa.GetExpectedSignal(t));
      histos.fill(HIST("tof/pr"), t.p(), t.tofSignal() - evtime - responsePr.GetExpectedSignal(t));
      histos.fill(HIST("trk/pr"), t.p(), ttrk - evtime - responsePr.GetExpectedSignal(t));
      histos.fill(HIST("tof/de"), t.p(), t.tofSignal() - evtime - responseDe.GetExpectedSignal(t));
      histos.fill(HIST("trk/de"), t.p(), ttrk - evtime - responseDe.GetExpectedSignal(t));
      histos.fill(HIST("tof/tr"), t.p(), t.tofSignal() - evtime - responseTr.GetExpectedSignal(t));
      histos.fill(HIST("trk/tr"), t.p(), ttrk - evtime - responseTr.GetExpectedSignal(t));
      histos.fill(HIST("tof/he"), t.p(), t.tofSignal() - evtime - responseHe.GetExpectedSignal(t));
      histos.fill(HIST("trk/he"), t.p(), ttrk - evtime - responseHe.GetExpectedSignal(t));
      histos.fill(HIST("tof/al"), t.p(), t.tofSignal() - evtime - responseAl.GetExpectedSignal(t));
      histos.fill(HIST("trk/al"), t.p(), ttrk - evtime - responseAl.GetExpectedSignal(t));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<trackTime>(cfgc)};
  return workflow;
}
