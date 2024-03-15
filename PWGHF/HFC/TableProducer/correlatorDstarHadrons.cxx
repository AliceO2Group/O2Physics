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

/// \file correlatorDstarHadron.cxx
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

// O2
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// O2Physics
#include "Common/DataModel/Multiplicity.h"

// PWGHF
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// flaging a collision if D* meson is found.
struct HfCollisionSelector {
  Produces<aod::DmesonSelection> collisionWDstar;

  Configurable<bool> selectionFlagDstar{"selectionFlagDstar", true, "selection flag for Dstar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  using DstarCandidates = soa::Join<aod::HfCandDstar, aod::HfSelDstarToD0Pi>;
  using FilteredCandidates = soa::Filtered<DstarCandidates>;

  SliceCache cache;
  Preslice<DstarCandidates> perColDstarCand = aod::hf_cand::collisionId;

  // candidates who passed the slection criteria defined in "CandidateSelectionTables.h"
  Filter candidateFilter = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar;

  void processCollisionSelWDstar(aod::Collisions const& collisions,
                                 FilteredCandidates const& candidates)
  {

    for (const auto& collision : collisions) {
      bool isDstarFound = false;
      auto candidatesPerCol = candidates.sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      if (!(candidatesPerCol.size() > 0)) {
        collisionWDstar(isDstarFound); // compatible with collision table (filled collision by collision)
        continue;
      }
      for (const auto& candidate : candidatesPerCol) {
        auto yDstar = candidate.y(constants::physics::MassDStar);
        auto ptDstar = candidate.pt();
        if (yCandMax >= 0 && std::abs(yDstar) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0 && ptDstar < ptCandMin) {
          continue;
        }
        isDstarFound = true;
        break;
      } // candidate loop
      LOG(info) << "processCollisionSelWDstar: isDstarFound = " << isDstarFound;
      collisionWDstar(isDstarFound); // compatible with collision table (filled collision by collision)
    }                                // collision loop
  }
  PROCESS_SWITCH(HfCollisionSelector, processCollisionSelWDstar, "process only data for dstar hadron correlation", true);
};

struct HfCorrelatorDstarHadron {
  Produces<aod::DstarHadronPair> rowsDstarHadronPair;
  // Dstar candidate related configurable
  Configurable<bool> selectOnlyCollisionWDstar{"selectOnlyCollisionWDstar", true, " select on collisions which have atleast a Dstar candidate"};
  Configurable<bool> selectionFlagDstar{"selectionFlagDstar", true, "selection flag for Dstar"};
  Configurable<float> pTMinDstar{"pTMinDstar", 1.5, "min pT of dstar candidate"};
  Configurable<float> pTMaxDstar{"pTMaxDstar", 50, "max pT of dstar Candidate"};
  // Configurable<float> etaAbsMaxDstar{"etaAbsMaxDstar",1.0,"max Abs(eta) cut on Dstar candidate"};
  Configurable<float> yMaxDstar{"yMaxDstar", 0.8, "max. cand. rapidity"};
  // track related configurable
  Configurable<float> etaAbsMaxAssoTrack{"etaAbsMaxAssoTrack", 0.8, "max Abs(eta) cut on Associated Track"};
  Configurable<float> dcaxyMinAssoTrack{"dcaxyMinAssoTrack", 0.0, "min DCAxy of Associated Track"};
  Configurable<float> dcaxyMaxAssoTrack{"dcaxyMaxAssoTrack", 10.0, "max DCAxy of Associated Track"};
  Configurable<float> dcazMinAssoTrack{"dcazMinAssoTrack", 0.0, "min DCAz of Associated Track"};
  Configurable<float> dcazMaxAssoTrack{"dcazMaxAssoTrack", 10.0, "max DCAz of Associated Track"};
  Configurable<float> pTMinAssoTrack{"pTMinAssoTrack", 0.5, "min Pt of Associated Track"};
  Configurable<float> pTMaxAssoTrack{"pTMaxAssoTrack", 50.0, "max pT of Associated Track"};

  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};

  // ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>> binningScheme{{binsZVtx, binsMultiplicity},true};
  // ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>> binningScheme;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  BinningType binningScheme{{binsZVtx, binsMultiplicity}, true};

  // Collision Table
  using CollisionsWDstar = soa::Join<aod::Collisions, aod::Mults, aod::DmesonSelection>;
  using FilteredCollisions = soa::Filtered<CollisionsWDstar>;

  // candidate table
  using DstarCandidates = soa::Join<aod::HfCandDstar, aod::HfSelDstarToD0Pi>; // Added extra cloumns in HfCandDstar (Prong0Id, Prong1Id), so no need to add table HfCandD0Fromdstar
  using FilteredCandidates = soa::Filtered<DstarCandidates>;

  using FilteredTracks = soa::Filtered<aod::TracksWDca>;

  // collision table filter
  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == selectOnlyCollisionWDstar;
  // candidate filter
  Filter candidateFilter = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar;
  // track table filter
  Filter trackFilter = nabs(aod::track::eta) <= etaAbsMaxAssoTrack && aod::track::pt >= pTMinAssoTrack && aod::track::pt <= pTMaxAssoTrack &&
                       aod::track::dcaXY >= dcaxyMinAssoTrack && aod::track::dcaXY <= dcaxyMaxAssoTrack &&
                       aod::track::dcaZ >= dcazMinAssoTrack && aod::track::dcaZ <= dcazMaxAssoTrack;
  SliceCache cache;
  // Preslice<DstarCandidates> perColCandidates = aod::hf_cand::collisionId;
  Preslice<FilteredCandidates> perColCandidates = aod::hf_cand::collisionId;
  // Preslice<aod::TracksWDca> perColTracks = aod::track::collisionId;
  Preslice<FilteredTracks> perColTracks = aod::track::collisionId;

  HistogramRegistry registry{
    "registry",
    {{"hTriggerColCandPairCounts", "Counts of Trigger Collision, Trigger Candidates and Pair Counts", {HistType::kTH1F, {{3, 0.0, 3.0}}}}}};

  void init(InitContext&)
  {
    binningScheme = {{binsZVtx, binsMultiplicity}, true};
  }

  void processDataSameEvent(FilteredCollisions const& collisions, // only collisions who have altleast one D*
                            FilteredTracks const& tracks,
                            FilteredCandidates const& candidates,
                            aod::BCsWithTimestamps const&)
  {

    for (const auto& collision : collisions) {
      registry.fill(HIST("hTriggerColCandPairCounts"), 0); // counting trigger collision

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto timestamp = bc.timestamp();
      auto candidatesPerCol = candidates.sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto tracksPerCol = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      if (candidatesPerCol.size() && tracksPerCol.size() == 0) {
        continue;
      } // endif

      registry.fill(HIST("hTriggerColCandPairCounts"), 1); // counting number of trigger particle

      // Pair creation
      for (const auto& [triggerParticle, assocParticle] : soa::combinations(soa::CombinationsFullIndexPolicy(candidatesPerCol, tracksPerCol))) {
        auto gItriggerParticle = triggerParticle.globalIndex();
        auto gIassocParticle = assocParticle.globalIndex();

        // Track rejection based on daughter index
        if ((triggerParticle.prong0Id() == gIassocParticle) || (triggerParticle.prong1Id() == gIassocParticle) || (triggerParticle.prongPiId() == gIassocParticle)) {
          continue; // rejected pair if associated particle is same as any of daughter particle
        }           // endif

        // Trigger Particle Rejection
        if (triggerParticle.pt() > pTMaxDstar || triggerParticle.pt() < pTMinDstar) {
          continue;
        } // endif
        auto yDstar = triggerParticle.y(constants::physics::MassDStar);
        if (std::abs(yDstar) > yMaxDstar) {
          continue;
        } // endif

        registry.fill(HIST("hTriggerColCandPairCounts"), 2); // counting number of pairs

        auto binNumber = binningScheme.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));

        // Fill table
        if (triggerParticle.signSoftPi() > 0) { // Fill Dstar candidate
          rowsDstarHadronPair(collision.globalIndex(),
                              gItriggerParticle,
                              triggerParticle.phi(),
                              triggerParticle.eta(),
                              triggerParticle.pt(),
                              triggerParticle.invMassDstar(),
                              gIassocParticle,
                              assocParticle.phi(),
                              assocParticle.eta(),
                              assocParticle.pt(),
                              timestamp,
                              binNumber);
        } else { // Fill AntiDstar candidate
          rowsDstarHadronPair(collision.globalIndex(),
                              gItriggerParticle,
                              triggerParticle.phi(),
                              triggerParticle.eta(),
                              triggerParticle.pt(),
                              triggerParticle.invMassAntiDstar(),
                              gIassocParticle,
                              assocParticle.phi(),
                              assocParticle.eta(),
                              assocParticle.pt(),
                              timestamp,
                              binNumber);
        } // endif

      } // D-H pair loop

    } // collision loop

  } // processDataSameEvent
  PROCESS_SWITCH(HfCorrelatorDstarHadron, processDataSameEvent, "process only same event data", true);

  void processDataWithMixedEvent(FilteredCollisions const& collisions, // only collisions who have altleast one D*
                                 FilteredTracks const& tracks,
                                 FilteredCandidates const& candidates,
                                 aod::BCsWithTimestamps const&)
  {

    auto dstarHadronTuple = std::make_tuple(candidates, tracks);
    Pair<FilteredCollisions, FilteredCandidates, FilteredTracks, BinningType> pairData{binningScheme, 5, -1, collisions, dstarHadronTuple, &cache};

    for (const auto& [c1, candidatesPerCol, c2, tracksPerCol] : pairData) {

      auto bc = c2.bc_as<aod::BCsWithTimestamps>();
      auto timestamp = bc.timestamp();

      for (const auto& [triggerParticle, assocParticle] : soa::combinations(soa::CombinationsFullIndexPolicy(candidatesPerCol, tracksPerCol))) {

        auto gItriggerParticle = triggerParticle.globalIndex();
        auto gIassocParticle = assocParticle.globalIndex();

        auto yDstar = triggerParticle.y(constants::physics::MassDStar);
        if (std::abs(yDstar) > yMaxDstar) {
          continue;
        } // endif

        int binNumber = binningScheme.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
        // Fill table
        if (triggerParticle.signSoftPi() > 0) { // Fill Dstar candidate
          rowsDstarHadronPair(c2.globalIndex(), // taking c2, why not c1?
                              gItriggerParticle,
                              triggerParticle.phi(),
                              triggerParticle.eta(),
                              triggerParticle.pt(),
                              triggerParticle.invMassDstar(),
                              gIassocParticle,
                              assocParticle.phi(),
                              assocParticle.eta(),
                              assocParticle.pt(),
                              timestamp,
                              binNumber);
        } else { // Fill AntiDstar candidate
          rowsDstarHadronPair(c2.globalIndex(),
                              gItriggerParticle,
                              triggerParticle.phi(),
                              triggerParticle.eta(),
                              triggerParticle.pt(),
                              triggerParticle.invMassAntiDstar(),
                              gIassocParticle,
                              assocParticle.phi(),
                              assocParticle.eta(),
                              assocParticle.pt(),
                              timestamp,
                              binNumber);
        } // endif

      } // D-H loop

    } // Event Mixing loop

  } // processDataWithMixedEvent
  PROCESS_SWITCH(HfCorrelatorDstarHadron, processDataWithMixedEvent, "process only mixed events data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCollisionSelector>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDstarHadron>(cfgc)};
}