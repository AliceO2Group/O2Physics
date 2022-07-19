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

#include <cmath>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//AxisSpec PtAxis = {1001, -0.005, 10.005};

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct EffPtMFT {
  Service<TDatabasePDG> pdg;

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  ConfigurableAxis PtAxis{"PtAxis", {1001, -0.0005, 1.0005}, "pt axis for histograms"};

  HistogramRegistry registry{
    "registry",
    {
      {"TracksPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}, //

    } //
  };

  void init(InitContext&)
  {

    if (doprocessGen) {
      registry.add({"TracksToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}});
      registry.add({"TracksToPartPtEtaFakeMcColl", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}});
      registry.add({"TracksToPartPtEtaPrim", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}}); //
      registry.add({"TracksPtEtaGen", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}});
      registry.add({"TracksPtEtaGen_t", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, {18, -4.6, -1.}}}});
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks)
  {
    if (!useEvSel || (useEvSel && collision.sel8())) {

      if ((collision.posZ() < 5) && (collision.posZ() > -5)) {
        for (auto& track : tracks) {

          registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
        }
      }
    }
  }

  using Particles = aod::McParticles;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;

  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, soa::Filtered<Particles> const& particles)
  {
    bool atLeastOne = false;

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (auto& collision : collisions) {

      if (!useEvSel || (useEvSel && collision.sel8())) {
        atLeastOne = true;
      }
    }

    if ((mcCollision.posZ() < 5) && (mcCollision.posZ() > -5)) {
      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        auto charge = 0;
        if (p != nullptr) {
          charge = (int)p->Charge();
        }
        if (std::abs(charge) < 3.) {
          continue;
        }

        if (atLeastOne) {

          registry.fill(HIST("TracksPtEtaGen"), particle.pt(), particle.eta());
        }

        registry.fill(HIST("TracksPtEtaGen_t"), particle.pt(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(EffPtMFT, processGen, "Process generator-level info", false);

  void processGenPt(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    //In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead

    if (collision.has_mcCollision()) {
      if ((collision.mcCollision().posZ() < 5) && (collision.mcCollision().posZ() > -5)) {
        if (!useEvSel || (useEvSel && collision.sel8())) {

          for (auto& track : tracks) {
            if (!track.has_mcParticle()) {
              continue;
            }
            auto particle = track.mcParticle(); //this mcParticle doesn't necessarly come from the right mcCollision
            registry.fill(HIST("TracksToPartPtEta"), particle.pt(), particle.eta());
            if ((particle.mcCollisionId() != collision.mcCollision().globalIndex()) || (!particle.isPhysicalPrimary())) {
              if (particle.mcCollisionId() != collision.mcCollision().globalIndex()) {
                registry.fill(HIST("TracksToPartPtEtaFakeMcColl"), particle.pt(), particle.eta());
              }
              continue;
            }
            registry.fill(HIST("TracksToPartPtEtaPrim"), particle.pt(), particle.eta());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(EffPtMFT, processGenPt, "Process particle-level info of pt", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EffPtMFT>(cfgc)};
}
