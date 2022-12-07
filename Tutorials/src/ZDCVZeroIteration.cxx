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
/// \brief These example tasks show how to access the ZDC and FV0A information which belongs to a collision.
///        The association is made through the BC column (and in Run 3 may not be unique!)
///        This example accesses the collisions and the related FV0A information.
///        Note that one has to subscribe to aod::FV0As const& to load
///        the relevant data even if you access the data itself through collisionMatched.fv0a().
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Here the "sparse" matcher is used which means, there can be collisions without FV0A information
// To find out, collisionMatched.has_fv0a() has to be called. Otherwise collisionMatched.fv0a() will fail.
// NOTE: subscribing to Collisions separately will lead to a circular dependency due to forwarding
struct IterateV0 {
  void process(aod::CollisionMatchedRun2Sparse const& collisionMatched, aod::FV0As const&)
  {
    LOGF(info, "Vertex = %f", collisionMatched.posZ());
    if (collisionMatched.has_fv0a()) {
      auto fv0a = collisionMatched.fv0a();
      for (unsigned int i = 0; i < fv0a.amplitude().size(); ++i) {
        LOGF(info, "V0A channel %d: %f", fv0a.channel()[i], fv0a.amplitude()[i]);
      }
    } else {
      LOGF(info, "No V0A info");
    }
  }
};

// This example is identical to IterateV0, but uses the exclusive match. This means that collisions where any
// of the tables asked for in Run2MatchedExclusive (see AnalysisDataModel.h) are missing are not there.
// Therefore, the syntax is more complicated because we cannot join against Collision
// (the tables have different number of entries)
// Only to be used if one is sure that all your events have the desired information
struct IterateV0Exclusive {
  void process(aod::Run2MatchedExclusive::iterator const& matcher, aod::Collisions const&, aod::FV0As const&)
  {
    LOGF(info, "Vertex = %f", matcher.collision().posZ());
    auto fv0a = matcher.fv0a();
    for (unsigned int i = 0; i < fv0a.amplitude().size(); ++i) {
      LOGF(info, "V0A channel %d: %f", fv0a.channel()[i], fv0a.amplitude()[i]);
    }
  }
};

// This example builds on IterateV0 and in addition accesses also the tracks grouped to the specific collision.
// The tracks are directly accessed through its pointer as usual
struct IterateV0Tracks {
  void process(aod::CollisionMatchedRun2Sparse const& collisionMatched, aod::FV0As const&, aod::Tracks const& tracks)
  {
    LOGF(info, "Vertex = %f. %d tracks", collisionMatched.posZ(), tracks.size());
    if (collisionMatched.has_fv0a()) {
      auto fv0a = collisionMatched.fv0a();
      for (unsigned int i = 0; i < fv0a.amplitude().size(); ++i) {
        LOGF(info, "V0A channel %d: %f", fv0a.channel()[i], fv0a.amplitude()[i]);
      }
    } else {
      LOGF(info, "No V0A info");
    }
  }
};

// This example accesses V0 and ZDC information
struct IterateV0ZDC {
  void process(aod::CollisionMatchedRun2Sparse const& collisionMatched, aod::FV0As const&, aod::Zdcs const&)
  {
    LOGF(info, "Vertex = %f", collisionMatched.posZ());
    if (collisionMatched.has_fv0a()) {
      auto fv0a = collisionMatched.fv0a();
      for (unsigned int i = 0; i < fv0a.amplitude().size(); ++i) {
        LOGF(info, "V0A channel %d: %f", fv0a.channel()[i], fv0a.amplitude()[i]);
      }
    } else {
      LOGF(info, "No V0A info");
    }
    if (collisionMatched.has_zdc()) {
      LOGF(info, "ZDC: E1 = %.3f; E2 = %.3f", collisionMatched.zdc().energyZEM1(), collisionMatched.zdc().energyZEM2());
    } else {
      LOGF(info, "No ZDC info");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<IterateV0>(cfgc),
    adaptAnalysisTask<IterateV0Exclusive>(cfgc),
    adaptAnalysisTask<IterateV0Tracks>(cfgc),
    adaptAnalysisTask<IterateV0ZDC>(cfgc),
  };
}
