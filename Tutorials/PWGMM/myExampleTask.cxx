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
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisHelpers.h>

using namespace o2;
using namespace o2::framework;

namespace test_col
{
    DECLARE_SOA_COLUMN(Eta, eta, float);
    DECLARE_SOA_COLUMN(Phi, phi, float);
    DECLARE_SOA_COLUMN(Pt, pt, float);
}
DECLARE_SOA_TABLE(Test_table, "AOD", "Blabla",
                  test_col::Eta, test_col::Phi, test_col::Pt);


struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Produces<aod::Test_table> TestTable;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispT{nBinsPt, 0, 10, "p_{T}"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("pTHistogram", "pTHistogram", kTH1F, {axispT});
  }

  void process(aod::TracksIU const& tracks)
  {
      LOG(info) << "Hello";
      for (auto& track : tracks) {
          histos.fill(HIST("etaHistogram"), track.eta());
          histos.fill(HIST("pTHistogram"), track.pt());

          TestTable(track.eta(),track.phi(),track.pt());
      }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}
