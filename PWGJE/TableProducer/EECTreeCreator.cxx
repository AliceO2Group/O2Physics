#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>
#include "Framework/O2DatabasePDGPlugin.h"

#include <vector>
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace tree
{
DECLARE_SOA_COLUMN(ZVTX, zvtx, float);
DECLARE_SOA_COLUMN(WEIGHT, weight, float);
DECLARE_SOA_COLUMN(PTHAT, ptHat, float);
DECLARE_SOA_COLUMN(MULTIPLICITY, multiplicity, float);

DECLARE_SOA_COLUMN(DETPt, detPt, std::vector<float>);
DECLARE_SOA_COLUMN(DETEta, detEta, std::vector<float>);
DECLARE_SOA_COLUMN(DETPhi, detPhi, std::vector<float>);
DECLARE_SOA_COLUMN(DETCharge, detCharge, std::vector<int>);

DECLARE_SOA_COLUMN(TruthPt, truthPt, std::vector<float>);
DECLARE_SOA_COLUMN(TruthEta, truthEta, std::vector<float>);
DECLARE_SOA_COLUMN(TruthPhi, truthPhi, std::vector<float>);
DECLARE_SOA_COLUMN(TruthE, truthE, std::vector<float>);
DECLARE_SOA_COLUMN(PDG, pdg, std::vector<int>);
}

DECLARE_SOA_TABLE(TREE, "AOD", "TREE",
    tree::ZVTX,
    tree::WEIGHT,
    tree::PTHAT,
    tree::MULTIPLICITY,
    tree::DETPt,
    tree::DETEta,
    tree::DETPhi,
    tree::DETCharge,
    tree::TruthPt,
    tree::TruthEta,
    tree::TruthPhi,
    tree::TruthE,
    tree::PDG);
}

struct EECTreeCreatorTask
{
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "vertex Z cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", ""};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", ""};

  Produces<aod::TREE> tree;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  Preslice<aod::JMcParticles> particlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  bool isChargedParticle(int code)
  {
    const float chargeUnit = 3.;
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= chargeUnit;
  }

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections);
  }

  void processMC(aod::JetCollisionsMCD::iterator const& collision, aod::JetTracks const& tracks, aod::JMcParticles const& mcParticles,aod::JetMcCollisions const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) return;
    if (std::abs(collision.posZ()) > vertexZCut) return;

    float w = collision.has_mcCollision() ? collision.mcCollision().weight() : 1.f;
    float pthard = collision.has_mcCollision() ? collision.mcCollision().ptHard() : 1.f;

    std::vector<float> detPt, detEta, detPhi;
    std::vector<int> detCharge;

    std::vector<float> truthPt, truthEta, truthPhi, truthE;
    std::vector<int> pdg;

    for (auto const& track : tracks)
    {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) continue;

      if( fabs(track.eta()) > 0.9) continue;
      if( track.pt() < 0.15 ) continue;

      detPt.push_back(track.pt());
      detEta.push_back(track.eta());
      detPhi.push_back(track.phi());
      detCharge.push_back(track.sign());

    }

    int mcId = collision.has_mcCollision() ? collision.mcCollisionId() : -1;

    if (mcId >= 0)
    {
      auto particles = mcParticles.sliceBy(particlesPerMcCollision, mcId);

      for (auto const& p : particles)
      {
        if( !p.isPhysicalPrimary()) continue;
        if( fabs(p.eta()) > 0.9) continue;
        if( p.pt() < 0.15) continue;
        if( !isChargedParticle(p.pdgCode())) continue;
        
        truthPt.push_back(p.pt());
        truthEta.push_back(p.eta());
        truthPhi.push_back(p.phi());
        truthE.push_back(p.e());
        pdg.push_back(p.pdgCode());
      }
    }

    tree( collision.posZ(), w, pthard, collision.multFT0C(), detPt, detEta, detPhi, detCharge, truthPt, truthEta, truthPhi, truthE, pdg);}

  PROCESS_SWITCH(EECTreeCreatorTask, processMC, "MC processing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{ adaptAnalysisTask<EECTreeCreatorTask>(cfgc, TaskName{"eec-tree-creator-mc"})};
}