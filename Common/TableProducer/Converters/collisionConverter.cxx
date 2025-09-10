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
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Swaps covariance matrix elements if the data is known to be bogus (collision_000 is bogus)
struct collisionConverter {
  Produces<aod::Collisions_001> Collisions_001;

  Configurable<bool> doNotSwap{"doNotSwap", false, "simple pass-through"};
  Configurable<bool> debug{"debug", false, "flag to save debug histo"};
  Configurable<int> nbins{"nbins", 1, "number of bins in debug histo"};
  Configurable<float> tolerance{"tolerance", 1e-3, "Tolerance for CYY check"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisCYYdebug{nbins, -1.0f, +1.0f, ""};
    histos.add("hCYY", "hCYY", kTH1F, {axisCYYdebug});
  }

  void process(aod::Collisions_000 const& collisionTable)
  {
    float negtolerance = -1.0f * tolerance;
    for (auto& collision : collisionTable) {
      float lYY = collision.covXZ();
      float lXZ = collision.covYY();
      if (doNotSwap) {
        lYY = collision.covYY();
        lXZ = collision.covXZ();
      }
      if (debug)
        histos.fill(HIST("hCYY"), lYY);
      if (lYY < negtolerance) {
        // This happened by accident!
        if (!doNotSwap && !debug) {
          LOGF(info, "Collision converter task found negative YY element!");
          LOGF(info, "CYY = %.10f, exceeds tolerance of %.10f", lYY, negtolerance);
          LOGF(info, "This is an indication that you're looping over data");
          LOGF(info, "produced with an O2 version of late December 2022.");
          LOGF(info, "Unfortunately, O2 versions of late December 2022");
          LOGF(info, "have a mistake in them for which a special mode");
          LOGF(info, "of this task exists. ");
          LOGF(info, "For this data, please operate the collision converter");
          LOGF(info, "with the configurable 'doNotSwap' set to true.");
          LOGF(info, "This program will now crash. Please adjust your settings!");
          LOGF(fatal, "FATAL: please set doNotSwap to true!");
        }
        if (!debug) {
          LOGF(info, "Collision converter task found negative YY element!");
          LOGF(info, "CYY = %.10f, exceeds tolerance of %.10f", lYY, negtolerance);
          LOGF(info, "You're running with 'doNotSwap' enabled, but the ");
          LOGF(info, "data your're analysing requires it to be disabled. ");
          LOGF(info, "This program will now crash. Please adjust your settings!");
          LOGF(fatal, "FATAL: please set doNotSwap to false!");
        }
      }
      // Repopulate new table
      Collisions_001(
        collision.bcId(),
        collision.posX(), collision.posY(), collision.posZ(),
        collision.covXX(),
        collision.covXY(),
        lYY, // <- this is the fixed part
        lXZ, // <- this is the fixed part
        collision.covYZ(),
        collision.covZZ(),
        collision.flags(), collision.chi2(), collision.numContrib(),
        collision.collisionTime(), collision.collisionTimeRes());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<collisionConverter>(cfgc),
  };
}
