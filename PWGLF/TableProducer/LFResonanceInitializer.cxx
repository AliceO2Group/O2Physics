#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct mcqa {
  // necessary for particle charges
  Service<O2DatabasePDG> pdgc;
  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                             // mass of 323, 333, 9010221, 10221, 9030221, 10331, 113, 213, 3224, 3124, 3324, 10323, 123314, 123324
                              {"hMass323", "hMass323", {HistType::kTH1F, {{400, 0.7f, 1.1f}}}},
                              {"hMass333", "hMass333", {HistType::kTH1F, {{400, 0.7f, 1.1f}}}},
                              {"hMass9010221", "hMass9010221", {HistType::kTH1F, {{400, 0.8f, 1.2f}}}},
                              {"hMass10221", "hMass10221", {HistType::kTH1F, {{400, 1.1f, 1.5f}}}},
                              {"hMass9030221", "hMass9030221", {HistType::kTH1F, {{400, 1.3f, 1.7f}}}},
                              {"hMass10331", "hMass10331", {HistType::kTH1F, {{400, 1.5f, 1.9f}}}},
                              {"hMass113", "hMass113", {HistType::kTH1F, {{400, 0.5f, 0.9f}}}},
                              {"hMass213", "hMass213", {HistType::kTH1F, {{400, 0.5f, 0.9f}}}},
                              {"hMass3224", "hMass3224", {HistType::kTH1F, {{400, 1.2f, 1.6f}}}},
                              {"hMass3124", "hMass3124", {HistType::kTH1F, {{400, 1.3f, 1.7f}}}},
                              {"hMass3324", "hMass3324", {HistType::kTH1F, {{400, 1.3f, 1.7f}}}},
                              {"hMass10323", "hMass10323", {HistType::kTH1F, {{400, 1.0f, 1.4f}}}},
                              {"hMass123314", "hMass123314", {HistType::kTH1F, {{400, 1.7f, 2.1f}}}},
                              {"hMass123324", "hMass123324", {HistType::kTH1F, {{400, 1.7f, 2.1f}}}},
                              // pT of 323, 333, 9010221, 10221, 9030221, 10331, 113, 213, 3224, 3124, 3324, 10323, 123314, 123324
                              {"hpT323", "hpT323", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT333", "hpT333", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT9010221", "hpT9010221", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT10221", "hpT10221", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT9030221", "hpT9030221", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT10331", "hpT10331", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT113", "hpT113", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT213", "hpT213", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT3224", "hpT3224", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT3124", "hpT3124", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT3324", "hpT3324", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT10323", "hpT10323", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT123314", "hpT123314", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
                              {"hpT123324", "hpT123324", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}}
                              }};

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter mcFilter = (nabs(aod::mcparticle::pdgCode) == 313)        // K*
                                                    || (nabs(aod::mcparticle::pdgCode) == 323)     // K*pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 333)     // phi
                                                    || (nabs(aod::mcparticle::pdgCode) == 9010221) // f_0(980)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10221)   // f_0(1370)
                                                    || (nabs(aod::mcparticle::pdgCode) == 9030221) // f_0(1500)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10331)   // f_0(1710)
                                                    || (nabs(aod::mcparticle::pdgCode) == 113)     // rho(770)
                                                    || (nabs(aod::mcparticle::pdgCode) == 213)     // rho(770)pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 3224)    // Sigma(1385)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 3124)    // Lambda(1520)
                                                    || (nabs(aod::mcparticle::pdgCode) == 3324)    // Xi(1530)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 10323)   // K1(1270)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 123314)  // Xi(1820)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 123324); // Xi(1820)-0
  // Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == 313)        // K*
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 323)     // K*pm
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 333)     // phi
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 9010221) // f_0(980)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 10221)   // f_0(1370)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 9030221) // f_0(1500)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 10331)   // f_0(1710)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 113)     // rho(770)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 213)     // rho(770)pm
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 3224)    // Sigma(1385)+
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 3124)    // Lambda(1520)
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 3324)    // Xi(1530)0
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 10323)   // K1(1270)+
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 123314)  // Xi(1820)0
  //                                                   || (nabs(aod::mcparticle::pdgCode) == 123324); // Xi(1820)-0

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>::iterator const& collision,
               aod::McCollisions const& mcCols, soa::Filtered<aod::McParticles> const& mcParticles)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());

    // Loop over all MC particles
    for (auto& mcPart : mcParticles) {
      std::vector<int> daughterPDGs;
      if (mcPart.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }
      auto pdg = abs(mcPart.pdgCode());
      auto pdgInfo = pdgc->GetParticle(pdg);
      if (pdgInfo == nullptr) {
        continue;
      }
      // fill in case of 323, 333, 9010221, 10221, 9030221, 10331, 113, 213, 3224, 3124, 3324, 10323, 123314, 123324
      if (pdg == 323) {
        registry.fill(HIST("hMass323"), pdgInfo->Mass());
        registry.fill(HIST("hpT323"), mcPart.pt());
      }
      if (pdg == 333) {
        registry.fill(HIST("hMass333"), pdgInfo->Mass());
        registry.fill(HIST("hpT333"), mcPart.pt());
      }
      if (pdg == 9010221) {
        registry.fill(HIST("hMass9010221"), pdgInfo->Mass());
        registry.fill(HIST("hpT9010221"), mcPart.pt());
      }
      if (pdg == 10221) {
        registry.fill(HIST("hMass10221"), pdgInfo->Mass());
        registry.fill(HIST("hpT10221"), mcPart.pt());
      }
      if (pdg == 9030221) {
        registry.fill(HIST("hMass9030221"), pdgInfo->Mass());
        registry.fill(HIST("hpT9030221"), mcPart.pt());
      }
      if (pdg == 10331) {
        registry.fill(HIST("hMass10331"), pdgInfo->Mass());
        registry.fill(HIST("hpT10331"), mcPart.pt());
      }
      if (pdg == 113) {
        registry.fill(HIST("hMass113"), pdgInfo->Mass());
        registry.fill(HIST("hpT113"), mcPart.pt());
      }
      if (pdg == 213) {
        registry.fill(HIST("hMass213"), pdgInfo->Mass());
        registry.fill(HIST("hpT213"), mcPart.pt());
      }
      if (pdg == 3224) {
        registry.fill(HIST("hMass3224"), pdgInfo->Mass());
        registry.fill(HIST("hpT3224"), mcPart.pt());
      }
      if (pdg == 3124) {
        registry.fill(HIST("hMass3124"), pdgInfo->Mass());
        registry.fill(HIST("hpT3124"), mcPart.pt());
      }
      if (pdg == 3324) {
        registry.fill(HIST("hMass3324"), pdgInfo->Mass());
        registry.fill(HIST("hpT3324"), mcPart.pt());
      }
      if (pdg == 10323) {
        registry.fill(HIST("hMass10323"), pdgInfo->Mass());
        registry.fill(HIST("hpT10323"), mcPart.pt());
      }
      if (pdg == 123314) {
        registry.fill(HIST("hMass123314"), pdgInfo->Mass());
        registry.fill(HIST("hpT123314"), mcPart.pt());
      }
      if (pdg == 123324) {
        registry.fill(HIST("hMass123324"), pdgInfo->Mass());
        registry.fill(HIST("hpT123324"), mcPart.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mcqa>(cfgc, TaskName{"lf-mcqa"})};
}
