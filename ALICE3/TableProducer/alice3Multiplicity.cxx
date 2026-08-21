// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Multiplicity task for ALICE3
/// \file alice3Multiplicity.cxx

#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksAlice3 = soa::Join<aod::Tracks, aod::TracksDCA, o2::aod::TracksAlice3, aod::TracksExtraA3>;

constexpr float EtaHalf = 0.5;
constexpr float Eta1 = 1.0;

struct Alice3Multiplicity {
  Produces<aod::PVMults> multPV;
  Produces<aod::MultsGlobal> multGlobal;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> minEta{"minEta", -2.5f, "Minimum eta in range for global track counting"};
  Configurable<float> maxEta{"maxEta", 2.5f, "Maximum eta in range for global track counting"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.0025f, "Max DCAxy for global track counting"};
  Configurable<float> maxDCAz{"maxDCAz", 0.0025f, "Max DCAz for global track counting"};
  Configurable<int> minSiliconHits{"minSiliconHits", 5, "Minimum number of hits in silicon detector for global track counting"};
  Configurable<bool> requireReconstructed{"requireReconstructed", false, "Require track to be reconstructed for global track counting"};
  Configurable<bool> doQA{"doQA", true, "Fill QA histograms"};

  ConfigurableAxis axisMult{"axisMult", {10000, 0, 10000}, "Reconstructed tracks"};

  Filter trackFilter = (aod::track::eta >= minEta) && (aod::track::eta <= maxEta) && (nabs(aod::track::dcaXY) <= maxDCAxy) && (nabs(aod::track::dcaZ) <= maxDCAz) && (aod::track_alice3::nSiliconHits >= minSiliconHits) && (!requireReconstructed || aod::track_alice3::isReconstructed);

  void init(InitContext&)
  {
    if (doQA) {
      histos.add("multiplicity/nTracksPV", "nTracksPV", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksPVeta1", "nTracksPVeta1", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksPVetaHalf", "nTracksPVetaHalf", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksGlobal", "nTracksGlobal", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksGlobalPV", "nTracksGlobalPV", kTH1D, {axisMult});
    }
  }

  void processGlobalTracks(const aod::Collision& /*collision*/, const soa::Filtered<TracksAlice3>& tracks)
  {
    int globalTracks = 0;
    int globalTracksPV = 0;

    for (const auto& track : tracks) {
      ++globalTracks;
      if (track.isPVContributor())
        ++globalTracksPV;
    }

    if (doQA) {
      histos.fill(HIST("multiplicity/nTracksGlobal"), globalTracks);
      histos.fill(HIST("multiplicity/nTracksGlobalPV"), globalTracksPV);
    }

    multGlobal(globalTracks, globalTracksPV, 0, 0);
  }

  void processPV(const aod::Collision& /*collision*/, const TracksAlice3& tracks)
  {
    int numTracksPV = 0;
    int numTracksPVeta1 = 0;
    int numTracksPVetaHalf = 0;

    for (const auto& track : tracks) {
      if (track.isPVContributor()) {
        ++numTracksPV;
        if (std::abs(track.eta()) < Eta1)
          ++numTracksPVeta1;
        if (std::abs(track.eta()) < EtaHalf)
          ++numTracksPVetaHalf;
      }
    }

    if (doQA) {
      histos.fill(HIST("multiplicity/nTracksPV"), numTracksPV);
      histos.fill(HIST("multiplicity/nTracksPVeta1"), numTracksPVeta1);
      histos.fill(HIST("multiplicity/nTracksPVetaHalf"), numTracksPVetaHalf);
    }

    multPV(numTracksPV, numTracksPVeta1, numTracksPVetaHalf);
  }

  void processDummy(const aod::Collision&)
  {
    // do nothing
  }

  PROCESS_SWITCH(Alice3Multiplicity, processGlobalTracks, "Process global track counter", false);
  PROCESS_SWITCH(Alice3Multiplicity, processPV, "Process primary vertex contributor tracks", false);
  PROCESS_SWITCH(Alice3Multiplicity, processDummy, "Dummy proccess function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3Multiplicity>(cfgc)};
}
