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
// This code calculates output histograms for centrality calibration
// as well as vertex-Z dependencies of raw variables (either for calibration
// of vtx-Z dependencies or for the calibration of those).
//
// This task is not strictly necessary in a typical analysis workflow,
// except for centrality calibration! The necessary task is the multiplicity
// tables.
//
// Comments, suggestions, questions? Please write to:
// - victor.gonzalez@cern.ch
// - david.dobrigkeit.chinellato@cern.ch
//

#include "Common/DataModel/Centrality.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct CentQA {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisCentrality{205, 0, 205, "centrality"};

    // Base histograms
    histos.add("hV0M", "V0M centrality", kTH1F, {axisCentrality});
    histos.add("hV0A", "V0A centrality", kTH1F, {axisCentrality});
    histos.add("hCL0", "CL0 centrality", kTH1F, {axisCentrality});
    histos.add("hCL1", "CL1 centrality", kTH1F, {axisCentrality});
    histos.add("hRefMult5", "RefMult .5 centrality", kTH1F, {axisCentrality});
    histos.add("hRefMult8", "RefMult .8 centrality", kTH1F, {axisCentrality});
  }

  void processV0M(soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator const& col)
  {
    histos.fill(HIST("hV0M"), col.centRun2V0M());
  }
  void processV0A(soa::Join<aod::Collisions, aod::CentRun2V0As>::iterator const& col)
  {
    histos.fill(HIST("hV0A"), col.centRun2V0A());
  }
  void processCL0(soa::Join<aod::Collisions, aod::CentRun2CL0s>::iterator const& col)
  {
    histos.fill(HIST("hCL0"), col.centRun2CL0());
  }
  void processCL1(soa::Join<aod::Collisions, aod::CentRun2CL1s>::iterator const& col)
  {
    histos.fill(HIST("hCL1"), col.centRun2CL1());
  }
  void processRefMult5(soa::Join<aod::Collisions, aod::CentRun2RefMult5s>::iterator const& col)
  {
    histos.fill(HIST("hRefMult5"), col.centRun2RefMult5());
  }
  void processRefMult8(soa::Join<aod::Collisions, aod::CentRun2RefMult8s>::iterator const& col)
  {
    histos.fill(HIST("hRefMult8"), col.centRun2RefMult8());
  }

  PROCESS_SWITCH(CentQA, processV0M, "QA V0M centrality", true);
  PROCESS_SWITCH(CentQA, processV0A, "QA V0A centrality", false);
  PROCESS_SWITCH(CentQA, processCL0, "QA CL0 centrality", false);
  PROCESS_SWITCH(CentQA, processCL1, "QA CL1 centrality", false);
  PROCESS_SWITCH(CentQA, processRefMult5, "QA RefMult5 centrality", false);
  PROCESS_SWITCH(CentQA, processRefMult8, "QA RefMult8 centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentQA>(cfgc)};
}
