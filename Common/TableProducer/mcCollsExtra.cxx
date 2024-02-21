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
//
// Quick and dirty task to correlate MC <-> data
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::FT0Mults>;

struct mcCollisionExtra {
  Produces<aod::McCollsExtra> mcCollsExtra;

  void processNoCentrality(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions)
  {
    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionIndex = collision.globalIndex();
      }
    }
    mcCollsExtra(collisions.size(), bestCollisionIndex, 0.0f);
  }
  void processWithCentrality(aod::McCollision const& mcCollision, soa::SmallGroups<FullCollisions> const& collisions)
  {
    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    float bestCollisionCentFT0C = 100.5f;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionIndex = collision.globalIndex();
      }
    }
    mcCollsExtra(collisions.size(), bestCollisionIndex, bestCollisionCentFT0C);
  }

  PROCESS_SWITCH(mcCollisionExtra, processNoCentrality, "process as if real data", false);
  PROCESS_SWITCH(mcCollisionExtra, processWithCentrality, "process as if real data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcCollisionExtra>(cfgc)};
}
