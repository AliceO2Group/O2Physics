#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct histo {
  // histogram created with OutputObj<TH1F>
  // OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "psedorapidity", 200, -1., +1)};
  // OutputObj<TH1F> ptHistogram{TH1F("ptHistogram", "transverse momentum", 100, 0, 20)}; 

  OutputObj<TH1F> etaHistogram{"etaHistogram"};  // this is like prototype I think!!!
  OutputObj<TH1F> ptHistogram{"ptHistogram"};

  void init(InitContext const&)  // This is an example of initialization of histograms differently
  {
    // now we will complete the definintion of OutputObj
    etaHistogram.setObject(new TH1F("etaHistogram", " PseudoRapidity ", 200, -1, 1));
    ptHistogram.setObject(new TH1F("ptHistogram", " Transverse Momentum", 100, 0, 10));
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
    adaptAnalysisTask<histo>(cfgc)
  };
}
