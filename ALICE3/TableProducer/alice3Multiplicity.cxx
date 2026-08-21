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

#include "ALICE3/DataModel/collisionAlice3.h"
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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksAlice3 = soa::Join<aod::Tracks, aod::TracksDCA, o2::aod::TracksAlice3, aod::TracksExtraA3>;

struct Alice3Multiplicity {
  Produces<aod::PVMults> multPV;
  Produces<aod::MultsGlobal> multGlobal;
  Produces<aod::MultsMCAlice3> multMC;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> minEta{"minEta", -2.5f, "Minimum eta in range for global track counting"};
  Configurable<float> maxEta{"maxEta", 2.5f, "Maximum eta in range for global track counting"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.0025f, "Max DCAxy for global track counting"};
  Configurable<float> maxDCAz{"maxDCAz", 0.0025f, "Max DCAz for global track counting"};
  Configurable<int> minSiliconHits{"minSiliconHits", 5, "Minimum number of hits in silicon detector for global track counting"};
  Configurable<bool> requireReconstructed{"requireReconstructed", false, "Require track to be reconstructed for global track counting"};
  Configurable<bool> doQA{"doQA", true, "Fill QA histograms"};

  ConfigurableAxis axisMult{"axisMult", {10000, 0, 10000}, "Reconstructed tracks"};

  Service<o2::framework::O2DatabasePDG> pdg;

  Filter trackFilter = (aod::track::eta >= minEta) && (aod::track::eta <= maxEta) && (nabs(aod::track::dcaXY) <= maxDCAxy) && (nabs(aod::track::dcaZ) <= maxDCAz) && (aod::track_alice3::nSiliconHits >= minSiliconHits) && (!requireReconstructed || aod::track_alice3::isReconstructed);

  void init(InitContext&)
  {
    if (doQA) {
      histos.add("multiplicity/nTracksPV", "nTracksPV", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksPVeta1", "nTracksPVeta1", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksPVetaHalf", "nTracksPVetaHalf", kTH1D, {axisMult});

      histos.add("multiplicity/nTracksGlobal", "nTracksGlobal", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksGlobalPV", "nTracksGlobalPV", kTH1D, {axisMult});

      histos.add("multiplicity/nTracksMC", "nTracksMC", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksMCEta25", "nTracksMCEta25", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksMCEta125", "nTracksMCEta125", kTH1D, {axisMult});
      histos.add("multiplicity/nTracksMCEta09", "nTracksMCEta09", kTH1D, {axisMult});
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
        if (std::abs(track.eta()) < 1.0)
          ++numTracksPVeta1;
        if (std::abs(track.eta()) < 0.5)
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

  void processMC(aod::McCollision const& /*mcCollision*/, aod::McParticles const& mcParticles)
  {
    int numMCParticles = 0;
    int numMCParticlesEta25 = 0;
    int numMCParticlesEta125 = 0;
    int numMCParticlesEta09 = 0;

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      auto charge = 0.;
      auto* p = pdg->GetParticle(mcParticle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 1e-3) {
        continue;
      }

      ++numMCParticles;
      if (std::abs(mcParticle.eta()) < 2.5)
        ++numMCParticlesEta25;
      if (std::abs(mcParticle.eta()) < 1.25)
        ++numMCParticlesEta125;
      if (std::abs(mcParticle.eta()) < 0.9)
        ++numMCParticlesEta09;
    }

    if (doQA) {
      histos.fill(HIST("multiplicity/nTracksMC"), numMCParticles);
      histos.fill(HIST("multiplicity/nTracksMCEta25"), numMCParticlesEta25);
      histos.fill(HIST("multiplicity/nTracksMCEta125"), numMCParticlesEta125);
      histos.fill(HIST("multiplicity/nTracksMCEta09"), numMCParticlesEta09);
    }

    multMC(numMCParticles, numMCParticlesEta25, numMCParticlesEta125, numMCParticlesEta09);
  }

  void processDummy(const aod::Collision&)
  {
    // do nothing
  }

  PROCESS_SWITCH(Alice3Multiplicity, processGlobalTracks, "Process global track counter", false);
  PROCESS_SWITCH(Alice3Multiplicity, processPV, "Process primary vertex contributor tracks", false);
  PROCESS_SWITCH(Alice3Multiplicity, processMC, "Process MC truth information", false);
  PROCESS_SWITCH(Alice3Multiplicity, processDummy, "Dummy proccess function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3Multiplicity>(cfgc)};
}
