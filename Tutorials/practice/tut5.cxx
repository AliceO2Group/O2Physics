#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

// In this program we will subscribe to more data tables.
struct moretables {
  Configurable<int> etabin{"etabin", 500, "N bins in eta histogram"};
  Configurable<int> ptbin{"ptbin", 200, " No. of bins in pt histogram"};

  HistogramRegistry registry{
    "registry",
    {{"etaHistogram", " pseudorapidity", {HistType::kTH1F, {{etabin, -1, 1}}}},
     {"ptHistogram", "transverse momentum", {HistType::kTH1F, {{ptbin, 0, 10}}}}}};
     
  void process(aod::Collision const& collision, soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks)
  {
    // Fill the counter here
    // Check the getter in this link:  https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))-> Fill(collision.posZ());
    for (auto& track : tracks) {
        if(track.tpcNClsCrossedRows() <70) continue;  // this is to skip stuff not tracked well by TPC
      registry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<moretables>(cfgc)};
}
