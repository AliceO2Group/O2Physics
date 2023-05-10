#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
using namespace o2;
using namespace o2::framework;
struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
  }
  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      histos.fill(HIST("etaHistogram"), track.eta());
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)
  };
}