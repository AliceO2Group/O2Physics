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

/// \file jetmatchingmc.cxx
/// \brief matching detector level and generator level jets
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_
#define PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_

#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#include <vector>

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename CandidatesBase, typename CandidatesTag, typename ClustersBase>
struct JetMatchingMc {

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  o2::framework::Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  o2::framework::Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsTag>();

  o2::framework::Preslice<JetsBase> baseJetsPerCollision = jetsBaseIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;
  o2::framework::Preslice<JetsTag> tagJetsPerCollision = jetsTagIsMc ? o2::aod::jet::mcCollisionId : o2::aod::jet::collisionId;

  o2::framework::PresliceUnsorted<o2::aod::JetCollisionsMCD> CollisionsPerMcCollision = o2::aod::jmccollisionlb::mcCollisionId;

  void init(o2::framework::InitContext const&)
  {
  }

  void processJets(o2::aod::JetMcCollisions const& mcCollisions, o2::aod::JetCollisionsMCD const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   o2::aod::JetTracksMCD const& tracks,
                   ClustersBase const& clusters,
                   o2::aod::JetParticles const& particles,
                   CandidatesBase const& candidatesBase,
                   CandidatesTag const& candidatesTag)
  {
    // initialise objects used to store the matching index arrays (array in case a mcCollision is split) before filling the matching tables
    std::vector<std::vector<int>> jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF;
    std::vector<std::vector<int>> jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF;
    //  waiting for framework fix to make sliced collection of same type as original collection:
    jetsBasetoTagMatchingGeo.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingPt.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingHF.assign(jetsBase.size(), {});
    jetsTagtoBaseMatchingGeo.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingPt.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingHF.assign(jetsTag.size(), {});

    for (const auto& mcCollision : mcCollisions) {

      const auto collisionsPerMcColl = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());

      for (const auto& collision : collisionsPerMcColl) {

        const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, jetsBaseIsMc ? mcCollision.globalIndex() : collision.globalIndex());
        const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMc ? mcCollision.globalIndex() : collision.globalIndex());

        jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidatesBase, tracks, clusters, candidatesTag, particles, particles, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
      }
    }
    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingMc, processJets, "Perform jet matching", true);
};

#endif // PWGJE_TABLEPRODUCER_MATCHING_JETMATCHINGMC_H_
