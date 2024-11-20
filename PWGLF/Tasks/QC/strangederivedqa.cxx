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
// V0 analysis task
// ================
//
// This code does basic QA of strangeness derived data

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using std::array;

struct strangederivedqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis axisNCollisions{"axisNCollisions", {50000, -0.5f, 49999.5f}, "collisions"};
  ConfigurableAxis axisNV0s{"axisNV0s", {50000, -0.5f, 49999.5f}, "V0s"};

  Configurable<bool> verbose{"verbose", false, "do more printouts"};

  void init(InitContext const&)
  {
    auto h = histos.add<TH1>("hDFCounter", "hDFCounter", kTH1D, {{6, -0.5f, 5.5f}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "Ordered");
    h->GetXaxis()->SetBinLabel(3, "Unordered");

    auto h2 = histos.add<TH2>("hEventCounter", "hEventCounter", kTH2D, {{1, -0.5f, 0.5f}, {3, -0.5f, 2.5f}});
    auto h3 = histos.add<TH2>("hEventsPerDF", "hEventsPerDF", kTH2D, {axisNCollisions, {3, -0.5f, 2.5f}});
    auto h4 = histos.add<TH2>("hV0sPerDF", "hV0sPerDF", kTH2D, {axisNV0s, {3, -0.5f, 2.5f}});

    h2->GetYaxis()->SetBinLabel(1, "All");
    h2->GetYaxis()->SetBinLabel(2, "Ordered");
    h2->GetYaxis()->SetBinLabel(3, "Unordered");
    h3->GetYaxis()->SetBinLabel(1, "All");
    h3->GetYaxis()->SetBinLabel(2, "Ordered");
    h3->GetYaxis()->SetBinLabel(3, "Unordered");
    h4->GetYaxis()->SetBinLabel(1, "All");
    h4->GetYaxis()->SetBinLabel(2, "Ordered");
    h4->GetYaxis()->SetBinLabel(3, "Unordered");
  }

  // Real data processing
  void processOriginal(aod::Collisions const& collisions, aod::Origins const& origins, soa::Join<aod::V0Indices, aod::V0Cores> const& fullV0s)
  {
    histos.fill(HIST("hDFCounter"), 0.0f);
    histos.fill(HIST("hEventCounter"), 0.0f, 0.0f, collisions.size());
    histos.fill(HIST("hEventsPerDF"), collisions.size(), 0.0f);
    histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 0.0f);
    bool ordered = true;
    int previousIndex = -100;
    for (auto const& v0 : fullV0s) {
      if (v0.collisionId() < previousIndex) {
        ordered = false;
      }
      previousIndex = v0.collisionId();
    }
    if (ordered) {
      histos.fill(HIST("hEventCounter"), 0.0f, 1.0f, collisions.size());
      histos.fill(HIST("hEventsPerDF"), collisions.size(), 1.0f);
      histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 1.0f);

      if (verbose) {
        auto origin = origins.begin();
        LOGF(info, "Sorted DF ID: %lld collisions: %i V0s: %i", origin.dataframeID(), collisions.size(), fullV0s.size());
      }
    } else {
      histos.fill(HIST("hEventCounter"), 0.0f, 2.0f, collisions.size());
      histos.fill(HIST("hEventsPerDF"), collisions.size(), 2.0f);
      histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 2.0f);

      if (verbose) {
        auto origin = origins.begin();
        LOGF(info, "Unsorted DF ID: %lld collisions: %i V0s: %i", origin.dataframeID(), collisions.size(), fullV0s.size());
      }
    }
  }

  // Real data processing
  void processDerived(aod::StraCollisions const& collisions, aod::StraOrigins const& origins, soa::Join<aod::V0CollRefs, aod::V0Cores> const& fullV0s)
  {
    histos.fill(HIST("hDFCounter"), 0.0f);
    histos.fill(HIST("hEventCounter"), 0.0f, 0.0f, collisions.size());
    histos.fill(HIST("hEventsPerDF"), collisions.size(), 0.0f);
    histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 0.0f);
    bool ordered = true;
    int previousIndex = -100;
    for (auto const& v0 : fullV0s) {
      if (v0.straCollisionId() < previousIndex) {
        ordered = false;
      }
      previousIndex = v0.straCollisionId();
    }
    if (ordered) {
      histos.fill(HIST("hEventCounter"), 0.0f, 1.0f, collisions.size());
      histos.fill(HIST("hEventsPerDF"), collisions.size(), 1.0f);
      histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 1.0f);

      if (verbose) {
        auto origin = origins.begin();
        LOGF(info, "Sorted DF ID: %lld collisions: %i V0s: %i Origins size: %i", origin.dataframeID(), collisions.size(), fullV0s.size(), origins.size());
      }
    } else {
      histos.fill(HIST("hEventCounter"), 0.0f, 2.0f, collisions.size());
      histos.fill(HIST("hEventsPerDF"), collisions.size(), 2.0f);
      histos.fill(HIST("hV0sPerDF"), fullV0s.size(), 2.0f);

      if (verbose) {
        auto origin = origins.begin();
        LOGF(info, "Unsorted DF ID: %lld collisions: %i V0s: %i Origins size: %i", origin.dataframeID(), collisions.size(), fullV0s.size(), origins.size());
        uint64_t directoryName = origin.dataframeID();
        for (auto const& orig : origins) {
          LOGF(info, "Unsorted DF ID: %lld separate origin: %lld", directoryName, orig.dataframeID());
        }
      }
    }
  }

  PROCESS_SWITCH(strangederivedqa, processOriginal, "Process original data", false);
  PROCESS_SWITCH(strangederivedqa, processDerived, "Process derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangederivedqa>(cfgc)};
}
