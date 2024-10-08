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
using std::array;

using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs>;

struct strangederivedqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    ConfigurableAxis axisNCollisions{"axisNCollisions", {50000, -0.5f, 49999.5f}, "collisions"};

  void init(InitContext const&)
  {
    auto h = histos.add<TH1>("hDFCounter", "hDFCounter", kTH1D, {{6, -0.5f, 5.5f}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "Ordered");
    h->GetXaxis()->SetBinLabel(3, "Unordered");

    auto h = histos.add<TH1>("hEventCounter", "hEventCounter", kTH1D, {{1, -0.5f, 0.5f}});
    auto h = histos.add<TH1>("hEventsPerDF", "hEventsPerDF", kTH1D, {axisNCollisions});
  }

_____________________________________________
  // Real data processing - no MC subscription
  void process(aod::StraCollisions const& collisions, v0Candidates const& fullV0s)
  {
    histos.fill(HIST("hDFCounter"), 0.0f);
    histos.fill(HIST("hEventCounter"), 0.0, collisions.size());
    histos.fill(HIST("hEventsPerDF"), collisions.size());
    bool ordered = true;
    int previousIndex = -100;
    for (auto const& v0 : fullV0s) {
      if(v0.straCollisionId() < previousIndex){ 
        ordered = false;
      }
      previousIndex = v0.straCollisionId(); 
    }
    if(ordered){ 
      histos.fill(HIST("hDFCounter"), 1.0f);
    }else{ 
      histos.fill(HIST("hDFCounter"), 2.0f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangederivedqa>(cfgc)};
}
