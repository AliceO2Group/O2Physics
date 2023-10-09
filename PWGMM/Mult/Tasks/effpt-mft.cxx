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
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"
#include "Index.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec EtaAxis = {18, -4.6, -1.};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec ZAxis = {301, -30.1, 30.1};

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct EffPtMFT {
  SliceCache cache;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  ConfigurableAxis PtAxis{"PtAxis", {1001, -0.0005, 1.0005}, "pt axis for histograms"};
  Configurable<float> zMax{"zMax", 5., "value for Zvtx cut"};

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

    if (doprocessTrackEfficiencyIndexed) {
      registry.add({"TracksPtEtaGenI", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"TracksPtEtaZvtxGenI", " ; p_{T} (GeV/c); #eta; #it{z}_{vtx} (cm)", {HistType::kTH3F, {PtAxis, EtaAxis, ZAxis}}});
      registry.add({"TracksPtEtaPrimI", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"TracksPtEtaZvtxPrimI", " ; p_{T} (GeV/c); #eta; #it{z}_{vtx} (cm)", {HistType::kTH3F, {PtAxis, EtaAxis, ZAxis}}});
      registry.add({"TracksPtEtaDuplicates", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"TracksPhiEtaGenDuplicates", " ; #phi; #eta", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"TracksPhiEtaDuplicates", " ; #phi; #eta", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks)
  {
    if (!useEvSel || (useEvSel && collision.sel8())) {

      if ((collision.posZ() < zMax) && (collision.posZ() > -zMax)) {
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

    if ((mcCollision.posZ() < zMax) && (mcCollision.posZ() > -zMax)) {
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
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead

    if (collision.has_mcCollision()) {
      if ((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)) {
        if (!useEvSel || (useEvSel && collision.sel8())) {

          for (auto& track : tracks) {
            if (!track.has_mcParticle()) {
              continue;
            }
            auto particle = track.mcParticle(); // this mcParticle doesn't necessarily come from the right mcCollision
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

  // using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using ParticlesI = soa::Join<aod::McParticles, aod::ParticlesToMftTracks>;
  // expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<ParticlesI> primariesI = ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  void processTrackEfficiencyIndexed(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
                                     aod::McCollisions const&,
                                     ParticlesI const&,
                                     MFTTracksLabeled const& tracks)
  {
    for (auto& collision : collisions) {
      if (useEvSel && !collision.sel8()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision();
      auto particlesI = primariesI->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      particlesI.bindExternalIndices(&tracks);

      for (auto& particle : particlesI) {
        auto charge = 0.;
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          charge = p->Charge();
        }
        if (std::abs(charge) < 3.) {
          continue;
        }
        registry.fill(HIST("TracksPtEtaZvtxGenI"), particle.pt(), particle.eta(), mcCollision.posZ()); // ptEtaGenPrimary
        if ((mcCollision.posZ() > -12) && (mcCollision.posZ() < 12)) {
          registry.fill(HIST("TracksPtEtaGenI"), particle.pt(), particle.eta()); // ptEtaGenPrimary
        }

        if (particle.has_mfttracks()) {
          auto counted = false;
          auto counter = 0;
          auto relatedTracks = particle.mfttracks_as<MFTTracksLabeled>();
          for (auto& track : relatedTracks) {

            ++counter;

            if (!counted) // this particle was not already counted
            {
              registry.fill(HIST("TracksPtEtaZvtxPrimI"), particle.pt(), particle.eta(), mcCollision.posZ());
              if ((mcCollision.posZ() > -12) && (mcCollision.posZ() < 12)) {
                registry.fill(HIST("TracksPtEtaPrimI"), particle.pt(), particle.eta());
              }
              counted = true;
            }
            if (counter > 1) {
              registry.fill(HIST("TracksPtEtaDuplicates"), particle.pt(), particle.eta());
              registry.fill(HIST("TracksPhiEtaDuplicates"), track.phi(), track.eta());
            } // Adding TracksPtEtaDuplicates and TracksPtEtaPrimI gives you the total
          }

          if (relatedTracks.size() > 1) {
            registry.fill(HIST("TracksPhiEtaGenDuplicates"), particle.phi(), particle.eta());
          }
        } // the particle has a track
      }   // loop on particlesI
    }     // loop on collisions
  }       // end of processTrackEfficiencyIndexed

  PROCESS_SWITCH(EffPtMFT, processTrackEfficiencyIndexed, "Calculate tracking efficiency vs pt (indexed)", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EffPtMFT>(cfgc)};
}
