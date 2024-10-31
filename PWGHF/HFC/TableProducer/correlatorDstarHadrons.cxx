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

/// \file correlatorDstarHadrons.cxx
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
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// flaging a collision if D* meson is found.
struct HfCorrelatorDstarHadronsCollisionSelector {
  Produces<aod::DmesonSelection> collisionWDstar;

  Configurable<bool> selectionFlagDstar{"selectionFlagDstar", true, "selection flag for Dstar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  using DstarCandidates = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
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
      // LOG(info) << "processCollisionSelWDstar: isDstarFound = " << isDstarFound;
      collisionWDstar(isDstarFound); // compatible with collision table (filled collision by collision)
    }                                // collision loop
  }
  PROCESS_SWITCH(HfCorrelatorDstarHadronsCollisionSelector, processCollisionSelWDstar, "process only data for dstar hadron correlation", true);
};

struct HfCorrelatorDstarHadrons {
  Produces<aod::DstarHadronPair> rowsDstarHadronPair;
  Produces<aod::Dstar> rowsDstar;
  Produces<aod::Hadron> rowsAssoTrack;

  // Enable separate tables for Dstar and Track for offline Event mixing
  Configurable<bool> enableSeparateTables{"enableSeparateTables", false, "Enable separate tables for Dstar and Track for offline Event mixing"};

  // Dstar candidate related configurable
  Configurable<bool> selectOnlyCollisionWDstar{"selectOnlyCollisionWDstar", true, " select on collisions which have atleast a Dstar candidate"};
  Configurable<bool> selectionFlagDstar{"selectionFlagDstar", true, "selection flag for Dstar"};
  Configurable<float> ptDstarMin{"ptDstarMin", 1.5, "min pT of dstar candidate"};
  Configurable<float> ptDstarMax{"ptDstarMax", 50, "max pT of dstar Candidate"};
  // Configurable<float> etaAbsDstarMax{"etaAbsDstarMax",1.0,"max Abs(eta) cut on Dstar candidate"};
  Configurable<float> yAbsDstarMax{"yAbsDstarMax", 0.8, "max. cand. rapidity"};
  // track related configurable
  Configurable<float> etaAbsAssoTrackMax{"etaAbsAssoTrackMax", 0.8, "max Abs(eta) cut on Associated Track"};
  Configurable<float> dcaxyAssoTrackMin{"dcaxyAssoTrackMin", 0.0, "min DCAxy of Associated Track"};
  Configurable<float> dcaxyAssoTrackMax{"dcaxyAssoTrackMax", 10.0, "max DCAxy of Associated Track"};
  Configurable<float> dcazAssoTrackMin{"dcazAssoTrackMin", 0.0, "min DCAz of Associated Track"};
  Configurable<float> dcazAssoTrackMax{"dcazAssoTrackMax", 10.0, "max DCAz of Associated Track"};
  Configurable<float> ptAssoTrackMin{"ptAssoTrackMin", 0.5, "min Pt of Associated Track"};
  Configurable<float> ptAssoTrackMax{"ptAssoTrackMax", 50.0, "max pT of Associated Track"};

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
  using DstarCandidates = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>; // Added extra cloumns in HfCandDstar (Prong0Id, Prong1Id), so no need to add table HfCandD0Fromdstar
  using FilteredCandidates = soa::Filtered<DstarCandidates>;

  using FilteredTracks = soa::Filtered<aod::TracksWDca>;

  // collision table filter
  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == selectOnlyCollisionWDstar;
  // candidate filter
  Filter candidateFilter = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar;
  // track table filter
  Filter trackFilter = nabs(aod::track::eta) <= etaAbsAssoTrackMax && aod::track::pt >= ptAssoTrackMin && aod::track::pt <= ptAssoTrackMax &&
                       aod::track::dcaXY >= dcaxyAssoTrackMin && aod::track::dcaXY <= dcaxyAssoTrackMax &&
                       aod::track::dcaZ >= dcazAssoTrackMin && aod::track::dcaZ <= dcazAssoTrackMax;
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
    std::array<bool, 2> processes = {doprocessDataSameEvent, doprocessDataWithMixedEvent};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

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
        if (triggerParticle.pt() > ptDstarMax || triggerParticle.pt() < ptDstarMin) {
          continue;
        } // endif
        auto yDstar = triggerParticle.y(constants::physics::MassDStar);
        if (std::abs(yDstar) > yAbsDstarMax) {
          continue;
        } // endif

        registry.fill(HIST("hTriggerColCandPairCounts"), 2); // counting number of pairs

        auto binNumber = binningScheme.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));

        // Inv Mass of Dstar and D0 Candidate
        float invMassDstarParticle = -999.;
        float invMassD0Particle = -999.;

        if (triggerParticle.signSoftPi() > 0) {
          invMassDstarParticle = triggerParticle.invMassDstar();
          invMassD0Particle = triggerParticle.invMassD0();
        } else {
          invMassDstarParticle = triggerParticle.invMassAntiDstar();
          invMassD0Particle = triggerParticle.invMassD0Bar();
        }

        // Fill Tables
        rowsDstarHadronPair(collision.globalIndex(),
                            gItriggerParticle,
                            triggerParticle.phi(),
                            triggerParticle.eta(),
                            triggerParticle.pt(),
                            invMassDstarParticle,
                            invMassD0Particle,
                            gIassocParticle,
                            assocParticle.phi(),
                            assocParticle.eta(),
                            assocParticle.pt(),
                            timestamp,
                            binNumber);

        if (enableSeparateTables) {
          rowsDstar(collision.globalIndex(),
                    gItriggerParticle,
                    triggerParticle.phi(),
                    triggerParticle.eta(),
                    triggerParticle.pt(),
                    invMassDstarParticle,
                    invMassD0Particle,
                    timestamp,
                    binNumber);

          rowsAssoTrack(assocParticle.phi(),
                        assocParticle.eta(),
                        assocParticle.pt(),
                        binNumber,
                        collision.globalIndex(),
                        timestamp);
        }

      } // D-H pair loop

    } // collision loop

  } // processDataSameEvent
  PROCESS_SWITCH(HfCorrelatorDstarHadrons, processDataSameEvent, "process only same event data", true);

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
        if (std::abs(yDstar) > yAbsDstarMax) {
          continue;
        } // endif

        int binNumber = binningScheme.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));

        // Inv Mass of Dstar and D0 Candidate
        float invMassDstarParticle = -999.;
        float invMassD0Particle = -999.;

        if (triggerParticle.signSoftPi() > 0) {
          invMassDstarParticle = triggerParticle.invMassDstar();
          invMassD0Particle = triggerParticle.invMassD0();
        } else {
          invMassDstarParticle = triggerParticle.invMassAntiDstar();
          invMassD0Particle = triggerParticle.invMassD0Bar();
        }

        // Fill Table
        rowsDstarHadronPair(c2.globalIndex(), // taking c2, why not c1?
                            gItriggerParticle,
                            triggerParticle.phi(),
                            triggerParticle.eta(),
                            triggerParticle.pt(),
                            invMassDstarParticle,
                            invMassD0Particle,
                            gIassocParticle,
                            assocParticle.phi(),
                            assocParticle.eta(),
                            assocParticle.pt(),
                            timestamp,
                            binNumber);

      } // D-H loop

    } // Event Mixing loop

  } // processDataWithMixedEvent
  PROCESS_SWITCH(HfCorrelatorDstarHadrons, processDataWithMixedEvent, "process only mixed events data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDstarHadronsCollisionSelector>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDstarHadrons>(cfgc)};
}
