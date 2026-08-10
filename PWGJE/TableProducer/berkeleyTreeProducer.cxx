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

/// \file berkeleyTreeProducer.cxx
/// \brief Task to save MC events into the BerkeleyTree format
/// \author Vassu Doomra <vdoomra@berkeley.edu>
/// \author Tucker Hwang <mhwang@cern.ch>

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
DECLARE_SOA_COLUMN(VtxZ, vtxZ, float);
DECLARE_SOA_COLUMN(Weight, weight, float);
DECLARE_SOA_COLUMN(PtHat, ptHat, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_BITMAP_COLUMN(Rct, rct, 32);

DECLARE_SOA_COLUMN(DetPt, detPt, std::vector<float>);
DECLARE_SOA_COLUMN(DetEta, detEta, std::vector<float>);
DECLARE_SOA_COLUMN(DetPhi, detPhi, std::vector<float>);
DECLARE_SOA_COLUMN(DetCharge, detCharge, std::vector<int>);

DECLARE_SOA_COLUMN(GenPt, genPt, std::vector<float>);
DECLARE_SOA_COLUMN(GenEta, genEta, std::vector<float>);
DECLARE_SOA_COLUMN(GenPhi, genPhi, std::vector<float>);
DECLARE_SOA_COLUMN(GenE, genE, std::vector<float>);
DECLARE_SOA_COLUMN(GenCharge, genCharge, std::vector<int>);
DECLARE_SOA_COLUMN(PdgId, pdgId, std::vector<int>);
}

DECLARE_SOA_TABLE(TREE, "AOD", "TREE",
    tree::VtxZ,
    tree::Weight,
    tree::PtHat,
    tree::Multiplicity,
    tree::Rct,
    tree::DetPt,
    tree::DetEta,
    tree::DetPhi,
    tree::DetCharge,
    tree::GenPt,
    tree::GenEta,
    tree::GenPhi,
    tree::GenE,
    tree::GenCharge,
    tree::PdgId);
}

struct BerkeleyTreeProducer
{
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "maximum Z vertex"};
  Configurable<float> etaMaxDet{"etaMaxDet", 0.9f, "maximum eta for det-level tracks"};
  Configurable<float> etaMaxGen{"etaMaxGen", 0.9f, "maximum eta for gen-level particles"};
  Configurable<float> ptMinDet{"ptMinDet", 0.15f, "minimum pt (GeV) for det-level tracks"};
  Configurable<float> ptMinGen{"ptMinGen", 0.15f, "minimum pt (GeV) for gen-level particles"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", ""};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", ""};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "skip MB gap events"};

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

  int getCharge(int code) {
    auto p = pdg->GetParticle(code);
    if (!p) {
      LOG(fatal) << "Cannot find particle with PDG code " << code;
      return 0;
    }
    auto charge = p->Charge() / 3.0;
    return std::lround(charge);
  }

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections);
  }

  void processMC(aod::JetCollisionsMCD::iterator const& collision, aod::JetTracks const& tracks, aod::JMcParticles const& mcParticles,aod::JetMcCollisions const&)
  {
    // do not do any RCT selections, will be done on analysis level
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, false, "", false, false)) return;
    if (std::abs(collision.posZ()) > vertexZCut) return;

    float w = collision.has_mcCollision() ? collision.mcCollision().weight() : 1.f;
    float pthard = collision.has_mcCollision() ? collision.mcCollision().ptHard() : 1.f;

    std::vector<float> detPt, detEta, detPhi;
    std::vector<int> detCharge;

    std::vector<float> genPt, genEta, genPhi, genE;
    std::vector<int> genCharge, pdgId;

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) continue;

      if (std::fabs(track.eta()) > etaMaxDet) continue;
      if (track.pt() < ptMinDet) continue;

      detPt.push_back(track.pt());
      detEta.push_back(track.eta());
      detPhi.push_back(track.phi());
      detCharge.push_back(track.sign());
    }

    int mcId = collision.has_mcCollision() ? collision.mcCollisionId() : -1;

    if (mcId >= 0) {
      auto particles = mcParticles.sliceBy(particlesPerMcCollision, mcId);

      for (auto const& p : particles) {
        if(!p.isPhysicalPrimary()) continue;
        if(std::fabs(p.eta()) > etaMaxGen) continue;
        if(p.pt() < ptMinGen) continue;
        if(!isChargedParticle(p.pdgCode())) continue;

        genPt.push_back(p.pt());
        genEta.push_back(p.eta());
        genPhi.push_back(p.phi());
        genE.push_back(p.e());
        genCharge.push_back(getCharge(p.pdgCode()));
        pdgId.push_back(p.pdgCode());
      }
    }

    tree( collision.posZ(), w, pthard, collision.multFT0C(), collision.rct_raw(), detPt, detEta, detPhi, detCharge, genPt, genEta, genPhi, genE, genCharge, pdgId);
  }

  PROCESS_SWITCH(BerkeleyTreeProducer, processMC, "MC processing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{ adaptAnalysisTask<BerkeleyTreeProducer>(cfgc)};
}