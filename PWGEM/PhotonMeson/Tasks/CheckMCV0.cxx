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

/// \brief check the v0 phase-space
/// dependencies: o2-analysis-lf-lambdakzeromcfinder
/// \author daiki.sekihata@cern.ch felix.schlepper@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <TMath.h>
#include <TVector2.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct CheckMCV0 {
  // Track selection
  Configurable<bool> ptLogAxis{"ptLogAxis", false, "Flag to use a log momentum axis"};
  Configurable<float> minpt{"minpt", 0.001, "min pt for track in GeV/c"};
  Configurable<float> maxpt{"maxpt", 999.f, "max pt for track in GeV/c"};
  Configurable<float> maxZ{"maxZ", 200.f, "max z for track"};
  Configurable<float> maxeta{"maxeta", 999.f, "eta acceptance"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<float> minX{"minX", 0.0, "minimum X (starting point of track X)"};
  Configurable<float> maxX{"maxX", 200.0, "maximum X (starting point of track X)"};
  Configurable<float> minY{"minY", 0.0, "minimum Y (starting point of track Y)"};
  Configurable<float> maxY{"maxY", 200.0, "maximum Y (starting point of track Y)"};
  Configurable<int> minTPCNCls{"minTPCNCls", 10, "minimum number TPC clusters"};

  // Filters
  using Tracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU>;
  using FilteredTracks = soa::Filtered<Tracks>;
  Filter trackPt = aod::track::pt > minpt&& aod::track::pt < maxpt;
  Filter trackZ = nabs(aod::track::z) < 200.f;
  Filter trackX = nabs(aod::track::x) > minX&& nabs(aod::track::x) < maxX;
  Filter trackY = nabs(aod::track::y) > minY&& nabs(aod::track::y) < maxY;
  Filter trackEta = nabs(aod::track::eta) < maxeta;
  Filter trackDCA = nabs(aod::track::dcaXY) > dcamin&& nabs(aod::track::dcaXY) < dcamax;

  // Histogram Parameters
  Configurable<int> tglNBins{"tglNBins", 500, "nBins for tgl"};
  Configurable<int> zNBins{"zNBins", 1000, "nBins for z"};
  Configurable<int> xNBins{"xNBins", 1000, "nBins for x"};
  Configurable<int> yNBins{"yNBins", 1000, "nBins for y"};
  Configurable<int> ptNBins{"ptNBins", 1000, "nBins for pt"};

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr std::array<std::string_view, 6> v0types = {"ITSTPC_ITSTPC", "TPConly_TPConly", "ITSonly_ITSonly", "ITSTPC_TPConly", "ITSTPC_ITSonly", "TPConly_ITSonly"};

  void init(InitContext const& /*unused*/)
  {
    const AxisSpec tglPos{tglNBins, -5.f, +5.f, "tan(#lambda) of e^{+}"};
    const AxisSpec tglEle{tglNBins, -5.f, +5.f, "tan(#lambda) of e^{-}"};
    const AxisSpec zEle{zNBins, -maxZ, maxZ, "z of e^{-}"};
    const AxisSpec zPos{zNBins, -maxZ, maxZ, "z of e^{+}"};
    const AxisSpec zPosEle{zNBins, -maxZ, maxZ, "z of e^{+} - z of e^{-}"};
    const AxisSpec xEle{xNBins, -maxX, maxX, "x of e^{-}"};
    const AxisSpec xPos{xNBins, -maxX, maxX, "x of e^{+}"};
    const AxisSpec yPos{yNBins, -maxY, maxY, "y of e^{+}"};
    const AxisSpec yPosEle{yNBins, -maxY, maxY, "y of e^{+} - y of e^{-}"};
    AxisSpec ptPos{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{+}"};
    AxisSpec ptEle{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{-}"};
    if (ptLogAxis) {
      ptPos.makeLogarithmic();
      ptEle.makeLogarithmic();
    }

    for (const auto& v0type : v0types) {
      LOGF(info, "Adding histograms for %s", v0type);
      registry.add(Form("%s/TglTgl", v0type.data()), "tan(#lambda) vs. tan(#lambda)", HistType::kTH2F, {tglPos, tglEle});
      registry.add(Form("%s/PtPt", v0type.data()), "p_{T} vs. p_{T}", HistType::kTH2F, {{ptPos}, {ptEle}});

      registry.add(Form("%s/XX", v0type.data()), "x vs. x", HistType::kTH2F, {{xPos, xEle}});
      registry.add(Form("%s/ZZ", v0type.data()), "z vs. z", HistType::kTH2F, {{zPos, zEle}});
      registry.add(Form("%s/YY", v0type.data()), "y vs. y", HistType::kTH2F, {{yPos, yPosEle}});
      registry.add(Form("%s/XY", v0type.data()), "x vs. y", HistType::kTH2F, {{xPos, yPosEle}});
      registry.add(Form("%s/XZ", v0type.data()), "x vs. Z", HistType::kTH2F, {{xPos, zPosEle}});

      registry.add(Form("%s/ZPosTglPos", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zPos}, {tglPos}});
      registry.add(Form("%s/ZEleTglPos", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zEle}, {tglPos}});
      registry.add(Form("%s/ZPosTglEle", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zPos}, {tglEle}});
      registry.add(Form("%s/ZEleTglEle", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zEle}, {tglEle}});

      registry.add(Form("%s/PtPosTglPos", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptPos}, {tglPos}});
      registry.add(Form("%s/PtPosTglEle", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptPos}, {tglEle}});
      registry.add(Form("%s/PtEleTglPos", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptEle}, {tglPos}});
      registry.add(Form("%s/PtEleTglEle", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptEle}, {tglEle}});

      // AntiCor
      registry.add(Form("%s/anti/TglTgl", v0type.data()), "tan(#lambda) vs. tan(#lambda)", HistType::kTH2F, {tglPos, tglEle});
      registry.add(Form("%s/anti/PtPt", v0type.data()), "p_{T} vs. p_{T}", HistType::kTH2F, {{ptPos}, {ptEle}});

      registry.add(Form("%s/anti/XX", v0type.data()), "x vs. x", HistType::kTH2F, {{xPos, xEle}});
      registry.add(Form("%s/anti/ZZ", v0type.data()), "z vs. z", HistType::kTH2F, {{zPos, zEle}});
      registry.add(Form("%s/anti/YY", v0type.data()), "y vs. y", HistType::kTH2F, {{yPos, yPosEle}});
      registry.add(Form("%s/anti/XY", v0type.data()), "x vs. y", HistType::kTH2F, {{xPos, yPosEle}});
      registry.add(Form("%s/anti/XZ", v0type.data()), "x vs. Z", HistType::kTH2F, {{xPos, zPosEle}});

      registry.add(Form("%s/anti/ZPosTglPos", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zPos}, {tglPos}});
      registry.add(Form("%s/anti/ZEleTglPos", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zEle}, {tglPos}});
      registry.add(Form("%s/anti/ZPosTglEle", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zPos}, {tglEle}});
      registry.add(Form("%s/anti/ZEleTglEle", v0type.data()), "z vs. tan(#lambda)", HistType::kTH2F, {{zEle}, {tglEle}});

      registry.add(Form("%s/anti/PtPosTglPos", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptPos}, {tglPos}});
      registry.add(Form("%s/anti/PtPosTglEle", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptPos}, {tglEle}});
      registry.add(Form("%s/anti/PtEleTglPos", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptEle}, {tglPos}});
      registry.add(Form("%s/anti/PtEleTglEle", v0type.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {{ptEle}, {tglEle}});
    }
    registry.add("V0Counter", "V0 counter", HistType::kTH1F, {{6, 0.5, 6.5}});
  }

  void processV0(aod::Collisions const& /*collisions*/, aod::V0s const& v0s, FilteredTracks const& /*tracks*/)
  {
    for (const auto& v0 : v0s) {
      auto pos = v0.posTrack_as<Tracks>(); // positive daughter
      auto ele = v0.negTrack_as<Tracks>(); // negative daughter
      if (!checkV0leg(pos) || !checkV0leg(ele)) {
        continue;
      }
      registry.fill(HIST("V0Counter"), 1);

#define FILL(name)                                                        \
  {                                                                       \
    registry.fill(HIST(name "/TglTgl"), pos.tgl(), ele.tgl());            \
    registry.fill(HIST(name "/PtPt"), pos.pt(), ele.pt());                \
    registry.fill(HIST(name "/ZPosTglPos"), pos.z(), pos.tgl());          \
    registry.fill(HIST(name "/ZEleTglPos"), ele.z(), pos.tgl());          \
    registry.fill(HIST(name "/ZPosTglEle"), pos.z(), ele.tgl());          \
    registry.fill(HIST(name "/ZEleTglEle"), ele.z(), ele.tgl());          \
    registry.fill(HIST(name "/PtPosTglPos"), pos.pt(), pos.tgl());        \
    registry.fill(HIST(name "/PtPosTglEle"), pos.pt(), ele.tgl());        \
    registry.fill(HIST(name "/PtEleTglPos"), ele.pt(), pos.tgl());        \
    registry.fill(HIST(name "/PtEleTglEle"), ele.pt(), ele.tgl());        \
    registry.fill(HIST(name "/XX"), pos.x(), ele.x());                    \
    registry.fill(HIST(name "/ZZ"), pos.z(), ele.z());                    \
    registry.fill(HIST(name "/YY"), pos.y(), ele.y());                    \
    registry.fill(HIST(name "/XY"), pos.x(), pos.y() - ele.y());          \
    registry.fill(HIST(name "/XZ"), pos.x(), pos.z() - ele.z());          \
    if (isAntiCorTgl(pos, ele)) {                                         \
      registry.fill(HIST(name "/anti/TglTgl"), pos.tgl(), ele.tgl());     \
      registry.fill(HIST(name "/anti/PtPt"), pos.pt(), ele.pt());         \
      registry.fill(HIST(name "/anti/ZPosTglPos"), pos.z(), pos.tgl());   \
      registry.fill(HIST(name "/anti/ZEleTglPos"), ele.z(), pos.tgl());   \
      registry.fill(HIST(name "/anti/ZPosTglEle"), pos.z(), ele.tgl());   \
      registry.fill(HIST(name "/anti/ZEleTglEle"), ele.z(), ele.tgl());   \
      registry.fill(HIST(name "/anti/PtPosTglPos"), pos.pt(), pos.tgl()); \
      registry.fill(HIST(name "/anti/PtPosTglEle"), pos.pt(), ele.tgl()); \
      registry.fill(HIST(name "/anti/PtEleTglPos"), ele.pt(), pos.tgl()); \
      registry.fill(HIST(name "/anti/PtEleTglEle"), ele.pt(), ele.tgl()); \
      registry.fill(HIST(name "/anti/XX"), pos.x(), ele.x());             \
      registry.fill(HIST(name "/anti/ZZ"), pos.z(), ele.z());             \
      registry.fill(HIST(name "/anti/YY"), pos.y(), ele.y());             \
      registry.fill(HIST(name "/anti/XY"), pos.x(), pos.y() - ele.y());   \
      registry.fill(HIST(name "/anti/XZ"), pos.x(), pos.z() - ele.z());   \
    }                                                                     \
  }

      if (isTPConly_TPConly(ele, pos)) {
        FILL("TPConly_TPConly");
      }
      if (isITSonly_ITSonly(ele, pos)) {
        FILL("ITSonly_ITSonly");
      }
      if (isTPConly_ITSonly(ele, pos)) {
        FILL("TPConly_ITSonly");
      }
      if (isITSTPC_ITSTPC(ele, pos)) {
        FILL("ITSTPC_ITSTPC");
      }
      if (isITSTPC_TPConly(ele, pos)) {
        FILL("ITSTPC_TPConly");
      }
      if (isITSTPC_ITSonly(ele, pos)) {
        FILL("ITSTPC_ITSonly");
      }
    }
  }
  PROCESS_SWITCH(CheckMCV0, processV0, "process reconstructed info", true);

  void processDummy(aod::V0s const& v0s) {}
  PROCESS_SWITCH(CheckMCV0, processDummy, "process dummy", false);

  // Templates
  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (track.pt() < minpt && track.pt() > maxpt) {
      return false;
    }
    if (abs(track.eta()) > maxeta) {
      return false;
    }

    if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < minTPCNCls) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isITSTPCMatchedTrack(TTrack const& track)
  {
    return track.hasITS() && track.hasTPC();
  }

  template <typename TTrack>
  bool isTPConlyTrack(TTrack const& track)
  {
    return !track.hasITS() && track.hasTPC();
  }

  template <typename TTrack>
  bool isITSonlyTrack(TTrack const& track)
  {
    return track.hasITS() && !track.hasTPC();
  }

  template <typename TTrack>
  bool isITSTPC_ITSTPC(TTrack const& track0, TTrack const& track1)
  {
    return isITSTPCMatchedTrack(track0) && isITSTPCMatchedTrack(track1);
  }

  template <typename TTrack>
  bool isITSTPC_TPConly(TTrack const& track0, TTrack const& track1)
  {
    return (isITSTPCMatchedTrack(track0) && isTPConlyTrack(track1)) || (isITSTPCMatchedTrack(track1) && isTPConlyTrack(track0));
  }

  template <typename TTrack>
  bool isITSTPC_ITSonly(TTrack const& track0, TTrack const& track1)
  {
    return (isITSTPCMatchedTrack(track0) && isITSonlyTrack(track1)) || (isITSTPCMatchedTrack(track1) && isITSonlyTrack(track0));
  }

  template <typename TTrack>
  bool isTPConly_TPConly(TTrack const& track0, TTrack const& track1)
  {
    return isTPConlyTrack(track0) && isTPConlyTrack(track1);
  }

  template <typename TTrack>
  bool isITSonly_ITSonly(TTrack const& track0, TTrack const& track1)
  {
    return isITSonlyTrack(track0) && isITSonlyTrack(track1);
  }

  template <typename TTrack>
  bool isTPConly_ITSonly(TTrack const& track0, TTrack const& track1)
  {
    return (isTPConlyTrack(track0) && isITSonlyTrack(track1)) || (isTPConlyTrack(track1) && isITSonlyTrack(track0));
  }

  template <typename TTrack>
  bool isAntiCorTgl(TTrack const& track0, TTrack const& track1)
  {
    return (track0.tgl() < 0.0 && track1.tgl() > 0.0) || (track1.tgl() < 0.0 && track0.tgl() > 0.0);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckMCV0>(cfgc, TaskName{"check-mc-v0"})};
}
