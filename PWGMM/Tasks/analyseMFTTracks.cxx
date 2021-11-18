#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;

using Particles = aod::McParticles;



//First approach to analysing an AO2D.root file
//written thanks to https://aliceo2group.github.io/analysis-framework/docs/tutorials/analysistask.html


HistogramRegistry registry
{
   "registry",
   {
     {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},            //
     {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},         //
   }                                                                                                                   //
 };

struct analyseMFTTracks : AnalysisTask
{
  //init()?
  void process(o2::aod::MFTTracks const& tracks)
    {
        for (auto& track : tracks)
        {
          registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
        }

    }
    //end of process

    void processGen(Particles const& particles)
      {
          for (auto& particle : particles)
          {
            registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
          }

      }
      //end of processGen
}
//end of MyTask

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyTask>(cfgc),
  };
}
