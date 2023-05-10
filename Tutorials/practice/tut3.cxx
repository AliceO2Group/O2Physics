#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

// In this program we will change the histogram binning via configurable.

struct configexample {
  // configurable for number of bins.
    Configurable<int> nBinsEta{"nBinsEta", 100, "n bins in eta histogram"};
    Configurable<int> binspt{"binspt", 100, "n bins in pt histogram"};

  OutputObj<TH1F> etaHistogram{"etaHistogram"}; // this is like prototype I think!!!
  OutputObj<TH1F> ptHistogram{"ptHistogram"};

  void init(InitContext const&) // This is an example of initialization of histograms differently
  {
    // now we will complete the definintion of OutputObj
    etaHistogram.setObject(new TH1F("etaHistogram", " PseudoRapidity ", nBinsEta, -1, 1));
    ptHistogram.setObject(new TH1F("ptHistogram", " Transverse Momentum", binspt, 0, 10));
  }

  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
      ptHistogram->Fill(track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<configexample>(cfgc)};
}
