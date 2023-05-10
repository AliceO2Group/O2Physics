#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

// In this program we will use Histogram registry to store the output.
struct histregistry {
  Configurable<int> etabin{"etabin", 500, "N bins in eta histogram"};
  Configurable<int> ptbin{"ptbin", 200, " No. of bins in pt histogram"};

  HistogramRegistry registry{
    "registry",
    {{"etaHistogram", " pseudorapidity", {HistType::kTH1F, {{etabin, -1, 1}}}},
     {"ptHistogram", "transverse momentum", {HistType::kTH1F, {{ptbin, 0, 10}}}}}};
     
  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      registry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<histregistry>(cfgc)};
}
