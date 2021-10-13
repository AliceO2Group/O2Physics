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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/T0infoTable.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "TH1F.h"
#include "TH2F.h"
using namespace o2;
using namespace o2::framework;

struct FT0InfoQaTask {
  OutputObj<TH1F> hT0A{TH1F("hT0A", "T0A; ns", 200, -1, 1)};
  OutputObj<TH1F> hT0C{TH1F("hT0C", "T0C; ns", 200, -1, 1)};
  OutputObj<TH1F> hT0AC{TH1F("hT0AC", "(T0C+T0A)/2; ns", 200, -1, 1)};
  OutputObj<TH1F> hT0res{TH1F("hT0res", "T0 resolution; ns", 200, -0.5, 0.5)};
  OutputObj<TH1F> hT0Vertex{TH1F("hT0vertex", "T0 vertex; cm", 200, -30, 30.)};
  OutputObj<TH2F> hVertex_T0_PV{TH2F(" hVertex_T0_PV", "T0 vertex vs Primary Vertex; T0 vertex, cm; PV, cm", 200, -30, 30., 200, -30, 30.)};
  OutputObj<TH2F> hResMult{TH2F("hResMult", "T0 resolution vs event multiplicity", 600, 0., 2000, 100, -500, 500.)};
  OutputObj<TH2F> hT0V0mult{TH2F("hT0V0mult", "T0 vs V0 multiplicity ;V0Mult;T0Mmult", 200, 0., 50000., 200, 0., 20000.)};

  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::T0info>::iterator const& col /*, aod::T0info const &ft0s*/)
  {
    int64_t foundFT0 = col.foundFT0();
    if (foundFT0 != -1) {
      LOGF(info, "multV0A=%5.0f;  multT0A=%5.0f; PV=%f; T0vertex=%f;T0A=%f; T0C=%f; T0AC=%f ", col.multV0A(), col.multT0A(), col.posZ(), col.t0A(), col.t0C(), col.t0AC(), col.t0vertex());
      hT0A->Fill(col.t0A());
      hT0C->Fill(col.t0C());
      hT0AC->Fill(col.t0AC());
      hT0Vertex->Fill(col.t0vertex());
      hVertex_T0_PV->Fill(col.t0vertex(), col.posZ());
      hT0res->Fill((col.t0A() - col.t0C()) / 2.);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FT0InfoQaTask>(cfgc, TaskName{"ft0-qa"})};
}
