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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsFiltered.h"
#include "Framework/runDataProcessing.h"
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

#define TWOPFILTERLOGCOLLISIONS info
#define TWOPFILTERLOGTRACKS info

struct TwoParticleCorrelationsFilter {
  Produces<aod::TwoPAcceptedCollisions> acceptedcollisions;
  Produces<aod::TwoPFilteredTracks> accepteddtracks;
  Produces<aod::TwoPAcceptedGenCollisions> acceptedgencollisions;
  Produces<aod::TwoPFilteredParticles> acceptedgentracks;

#include "skimmingconf.h"

  void init(InitContext const&)
  {
    using namespace twopfilter;
    using namespace twopskim;

    LOGF(info, "TwoParticleCorrelationsFilter::init()");
  }

  void processRun2(aod::TwoPSkimmedCollision const& collision, aod::TwoPSkimmedTracks const& tracks)
  {
    using namespace twopfilter;

    LOGF(TWOPFILTERLOGCOLLISIONS, "Received filtered collision with mask 0x%lx", collision.selflags());
  }
  PROCESS_SWITCH(TwoParticleCorrelationsFilter, processRun2, "Process Run 2 two particle correlations filtering", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelationsFilter>(cfgc)};
  return workflow;
}
