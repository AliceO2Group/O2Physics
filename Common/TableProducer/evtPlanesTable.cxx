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
/// \file   evtPlanesTable.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EvtPlanes.h"

// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;

struct evtPlanesTable {
  // Table.
  Produces<aod::EvtPlanes> evPlane;

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  // Helper variables.
  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
  } // void init(InitContext const&)

  void process(aod::Qvector const& qVec)
  {
    // Get the centrality bin, and fill the centrality distribution.
    int centBin = helperEP.GetCentBin(qVec.cent());
    if (centBin < 0 || centBin > 8) {
      return;
    }
    // Calculate the event plane for each detector, then save them in the
    // corresponding distribution. The order is the same as in detNames[].
    // TODO: Update the calculation of the event plane for the central barrel.

    float evtPlane[4];
    float evtPlaneBPos[4];
    float evtPlaneBNeg[4];

    evtPlane[0] = helperEP.GetEventPlane(qVec.qvecUncorRe(), qVec.qvecUncorIm());
    evtPlane[1] = helperEP.GetEventPlane(qVec.qvecRectrRe(), qVec.qvecRectrIm());
    evtPlane[2] = helperEP.GetEventPlane(qVec.qvecTwistRe(), qVec.qvecTwistIm());
    evtPlane[3] = helperEP.GetEventPlane(qVec.qvecFinalRe(), qVec.qvecFinalIm());

    evtPlaneBPos[0] = helperEP.GetEventPlane(qVec.qvecBPosUncorRe(), qVec.qvecBPosUncorIm());
    evtPlaneBPos[1] = helperEP.GetEventPlane(qVec.qvecBPosRectrRe(), qVec.qvecBPosRectrIm());
    evtPlaneBPos[2] = helperEP.GetEventPlane(qVec.qvecBPosTwistRe(), qVec.qvecBPosTwistIm());
    evtPlaneBPos[3] = helperEP.GetEventPlane(qVec.qvecBPosFinalRe(), qVec.qvecBPosFinalIm());

    evtPlaneBNeg[0] = helperEP.GetEventPlane(qVec.qvecBNegUncorRe(), qVec.qvecBNegUncorIm());
    evtPlaneBNeg[1] = helperEP.GetEventPlane(qVec.qvecBNegRectrRe(), qVec.qvecBNegRectrIm());
    evtPlaneBNeg[2] = helperEP.GetEventPlane(qVec.qvecBNegTwistRe(), qVec.qvecBNegTwistIm());
    evtPlaneBNeg[3] = helperEP.GetEventPlane(qVec.qvecBNegFinalRe(), qVec.qvecBNegFinalIm());

    // Fill the columns of the evtPlane table.
    evPlane(qVec.cent(),
            evtPlane[0], evtPlane[1], evtPlane[2], evtPlane[3],
            evtPlaneBPos[0], evtPlaneBPos[1], evtPlaneBPos[2], evtPlaneBPos[3],
            evtPlaneBNeg[0], evtPlaneBNeg[1], evtPlaneBNeg[2], evtPlaneBNeg[3],
            qVec.nTrkBPos(), qVec.nTrkBNeg());
  } // void process(aod::Qvector const& qVec)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<evtPlanesTable>(cfgc)};
}
