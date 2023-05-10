

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

//STEP 0
//This is an empty analysis skeleton: the starting point! 
struct histo {
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};
  OutputObj<TH1F> ptHistogram{TH1F("ptHistogram", "transverse momentum", 100, 0, 20)}; 


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
