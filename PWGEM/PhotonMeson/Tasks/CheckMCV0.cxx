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

/// \brief write relevant information for photon conversion analysis to an AO2D.root file.
/// dependencies: o2-analysis-lf-lambdakzeromcfinder
/// \author daiki.sekihata@cern.ch

// #include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
// #include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"
// #include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "DCAFitter/HelixHelper.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"

#include <TMath.h>
#include <TVector2.h>
#include "Math/Vector4D.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU>;
// using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;

struct CheckMCV0 {
  Configurable<float> minpt{"minpt", 0.01, "min pt for track in GeV/c"};
  Configurable<float> maxeta{"maxeta", 999.f, "eta acceptance"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    addhistograms();
  }

  static constexpr std::string_view v0types[6] = {"ITSTPC_ITSTPC", "TPConly_TPConly", "ITSonly_ITSonly", "ITSTPC_TPConly", "ITSTPC_ITSonly", "TPConly_ITSonly"};
  void addhistograms()
  {
    registry.add("V0Counter", "V0 counter", HistType::kTH1F, {{6, 0.5, 6.5}});
    for (int i = 0; i < 6; i++) {
      LOGF(info, "adding histogram for %s", v0types[i].data());
      registry.add(Form("%s/h2TglTgl", v0types[i].data()), "tgl vs. tgl;tan(#lambda) of e^{+};tan(#lambda) of e^{-}", HistType::kTH2F, {{100, -5, +5}, {100, -5, +5}});
    }
  }

  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (track.pt() < minpt || abs(track.eta()) > maxeta) {
      return false;
    }
    if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isITSTPCMatchedTrack(TTrack const& track)
  {
    return track.hasITS() & track.hasTPC();
  }

  template <typename TTrack>
  bool isTPConlyTrack(TTrack const& track)
  {
    return !track.hasITS() & track.hasTPC();
  }

  template <typename TTrack>
  bool isITSonlyTrack(TTrack const& track)
  {
    return track.hasITS() & !track.hasTPC();
  }

  // Preslice<aod::V0s> perCollision = aod::v0::collisionId;
  void processV0(aod::Collisions const& collisions, aod::V0s const& v0s, MyTracks const& tracks)
  {
    for (auto& v0 : v0s) {
      auto pos = v0.template posTrack_as<MyTracks>(); // positive daughter
      auto ele = v0.template negTrack_as<MyTracks>(); // negative daughter
      if (!checkV0leg(pos) || !checkV0leg(ele)) {
        continue;
      }
      registry.fill(HIST("V0Counter"), 1);

      if (isITSTPCMatchedTrack(pos) && isITSTPCMatchedTrack(ele)) {
        registry.fill(HIST("ITSTPC_ITSTPC/h2TglTgl"), pos.tgl(), ele.tgl());
      }
      if ((isITSTPCMatchedTrack(pos) && isTPConlyTrack(ele)) || (isITSTPCMatchedTrack(ele) && isTPConlyTrack(pos))) {
        registry.fill(HIST("ITSTPC_TPConly/h2TglTgl"), pos.tgl(), ele.tgl());
      }
      if ((isITSTPCMatchedTrack(pos) && isITSonlyTrack(ele)) || (isITSTPCMatchedTrack(ele) && isITSonlyTrack(pos))) {
        registry.fill(HIST("ITSTPC_ITSonly/h2TglTgl"), pos.tgl(), ele.tgl());
      }
      if (isTPConlyTrack(pos) && isTPConlyTrack(ele)) {
        registry.fill(HIST("TPConly_TPConly/h2TglTgl"), pos.tgl(), ele.tgl());
      }
      if (isITSonlyTrack(pos) && isITSonlyTrack(ele)) {
        registry.fill(HIST("ITSonly_ITSonly/h2TglTgl"), pos.tgl(), ele.tgl());
      }
      if ((isTPConlyTrack(pos) && isITSonlyTrack(ele)) || (isTPConlyTrack(ele) && isITSonlyTrack(pos))) {
        registry.fill(HIST("TPConly_ITSonly/h2TglTgl"), pos.tgl(), ele.tgl());
      }

    } // end of v0 loop
  }
  PROCESS_SWITCH(CheckMCV0, processV0, "process reconstructed info", true);

  void processDummy(aod::V0s const& v0s) {}
  PROCESS_SWITCH(CheckMCV0, processDummy, "process dummy", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckMCV0>(cfgc, TaskName{"check-mc-v0"})};
}
