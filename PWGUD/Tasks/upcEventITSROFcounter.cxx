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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <TH1D.h>
#include <TH2D.h>
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


struct UpcEventITSROFcounter {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {

    const AxisSpec axisNtracks{30, -0.5, 29.5};

    histos.add("Events/hCountCollisions", ";;Number of collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountUPCcollisions", ";;Number of UPC (mult < 17) collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});


  } // end init


  void process(aod::BC const& bc, aod::Collisions const& collisions)
  {

    LOGF(info,"BC number %i, numContrib %i",bc.globalBC(),collisions.numContrib());

  }

  PROCESS_SWITCH(UpcEventITSROFcounter, process, "Runs over BCs and count collisions in ITSROF (594 BCs)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcEventITSROFcounter>(cfgc, TaskName{"upc-event-itsrof-counter"})};
}
