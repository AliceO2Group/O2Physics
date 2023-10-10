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

/// \file jetmatching.cxx
/// \brief Unified implementation of jet matching based on different criteria
/// expanding on previously separate implementations of geometric matching
/// (by Raymond Ehlers) and heavy-flavour matching
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename BaseJetCollection, typename TagJetCollection,
          typename BaseToTagMatchingTable, typename TagToBaseMatchingTable, typename McParticles, typename HfCandidates>
struct JetMatching {

  Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  Configurable<bool> doMatchingHf{"doMatchingHf", true, "Enable HF matching"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.4f, "Max matching distance"};
  Configurable<float> minPtFraction{"minPtFraction", 0.f, "Minimum pt fraction for pt matching"};

  Produces<BaseToTagMatchingTable> jetsBaseToTag;
  Produces<TagToBaseMatchingTable> jetsTagToBase;

  // preslicing jet collections, only for MC-based collection
  static constexpr bool jetsBaseIsMC = o2::soa::relatedByIndex<aod::McCollisions, BaseJetCollection>();
  static constexpr bool jetsTagIsMC = o2::soa::relatedByIndex<aod::McCollisions, TagJetCollection>();
  Preslice<BaseJetCollection> baseJetsPerCollision = jetsBaseIsMC ? aod::jet::mcCollisionId : aod::jet::collisionId;
  Preslice<TagJetCollection> tagJetsPerCollision = jetsTagIsMC ? aod::jet::mcCollisionId : aod::jet::collisionId;

  using Collisions = soa::Join<aod::McCollisionLabels, aod::Collisions>;
  PresliceUnsorted<Collisions> CollisionsCollectionPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  using Tracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

  // initialise objects used to store the matching index arrays (array in case a mcCollision is split) before filling the matching tables
  std::vector<std::vector<int>> geojetidBaseToTag, ptjetidBaseToTag, hfjetidBaseToTag;
  std::vector<std::vector<int>> geojetidTagToBase, ptjetidTagToBase, hfjetidTagToBase;

  static constexpr int8_t getHfFlag()
  {
    if (std::is_same<BaseToTagMatchingTable, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;

    if (std::is_same<BaseToTagMatchingTable, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;

    if (std::is_same<BaseToTagMatchingTable, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_bplus::DecayType::BplusToD0Pi;

    return -1;
  }

  void init(InitContext const&)
  {
  }

  // function that does the geometric matching of jets from jetsBasePerColl and jets from jetsTagPerColl
  template <typename T, typename U>
  void MatchGeo(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int>& baseToTagGeo, std::vector<int>& tagToBaseGeo)
  {
    std::vector<double> jetsBasePhi;
    std::vector<double> jetsBaseEta;
    for (const auto& jet : jetsBasePerColl) {
      jetsBasePhi.emplace_back(jet.phi());
      jetsBaseEta.emplace_back(jet.eta());
    }
    std::vector<double> jetsTagPhi;
    std::vector<double> jetsTagEta;
    for (const auto& jet : jetsTagPerColl) {
      jetsTagPhi.emplace_back(jet.phi());
      jetsTagEta.emplace_back(jet.eta());
    }
    std::tie(baseToTagGeo, tagToBaseGeo) = JetUtilities::MatchJetsGeometrically(jetsBasePhi, jetsBaseEta, jetsTagPhi, jetsTagEta, maxMatchingDistance);
    LOGF(debug, "geometric matching: %d - %d jets", baseToTagGeo.size(), tagToBaseGeo.size());
    for (std::size_t i = 0; i < baseToTagGeo.size(); ++i) {
      LOGF(debug, "bjet %i -> %i", i, baseToTagGeo[i]);
    }
    for (std::size_t i = 0; i < tagToBaseGeo.size(); ++i) {
      LOGF(debug, "tjet %i -> %i", i, tagToBaseGeo[i]);
    }
  }

  // function that does the HF matching of jets from jetsBasePerColl and jets from jetsTagPerColl
  template <typename T, typename U>
  void MatchHF(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int>& baseToTagHF, std::vector<int>& tagToBaseHF)
  {
    int index_bjet = 0;
    int index_tjet = 0;
    for (const auto& bjet : jetsBasePerColl) {
      LOGF(debug, "jet index: %d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           bjet.index(), bjet.collisionId(), bjet.pt(), bjet.phi(), bjet.tracks().size(), bjet.hfcandidates().size());
      const auto hfcand = bjet.template hfcandidates_first_as<HfCandidates>();
      if (hfcand.flagMcMatchRec() & getHfFlag()) {
        const auto hfCandMC = hfcand.template prong0_as<Tracks>().template mcParticle_as<McParticles>();
        const auto hfCandMcId = hfCandMC.template mothers_first_as<McParticles>().globalIndex();
        for (const auto& tjet : jetsTagPerColl) {
          const auto cand = tjet.template hfcandidates_first_as<McParticles>();
          if (cand.globalIndex() == hfCandMcId) {
            LOGF(debug, "Found HF match: %d (pt %g) <-> %d (pt %g)",
                 bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
            baseToTagHF[index_bjet] = index_tjet;
            tagToBaseHF[index_tjet] = index_bjet;
          }
          index_tjet++;
        }
      }
      index_bjet++;
    }
  }

  template <typename T, typename U>
  float getPtSum_FirstArgIsMC(T const& btracks, U const& ttracks)
  {
    float ptSum = 0;
    for (const auto& btrack : btracks) {
      for (const auto& ttrack : ttracks) {
        if (ttrack.has_mcParticle() && btrack.globalIndex() == ttrack.template mcParticle_as<McParticles>().globalIndex()) {
          ptSum += ttrack.pt();
          break;
        }
      }
    }
    return ptSum;
  }
  template <typename T, typename U>
  float getPtSum_SecondArgIsMC(T const& btracks, U const& ttracks)
  {
    float ptSum = 0;
    for (const auto& btrack : btracks) {
      for (const auto& ttrack : ttracks) {
        if (btrack.has_mcParticle() && ttrack.globalIndex() == btrack.template mcParticle_as<McParticles>().globalIndex()) {
          ptSum += ttrack.pt();
          break;
        }
      }
    }
    return ptSum;
  }

  // function that does the pT matching of jets from jetsBasePerColl and jets from jetsTagPerColl
  template <typename T, typename U>
  void MatchPt(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int>& baseToTagPt, std::vector<int>& tagToBasePt)
  {
    float ptSum;
    int index_bjet = 0;
    int index_tjet = 0;

    for (const auto& bjet : jetsBasePerColl) {
      for (const auto& tjet : jetsTagPerColl) {
        auto btracksMcParts = bjet.template tracks_as<McParticles>();
        auto btracksTracks = bjet.template tracks_as<Tracks>();
        auto ttracksMcParts = tjet.template tracks_as<McParticles>();
        auto ttracksTracks = tjet.template tracks_as<Tracks>();

        jetsBaseIsMC ? ptSum = getPtSum_FirstArgIsMC(btracksMcParts, ttracksTracks) : ptSum = getPtSum_SecondArgIsMC(btracksTracks, ttracksMcParts);

        if (ptSum > tjet.pt() * minPtFraction) {
          LOGF(debug, "Found pt match: %d (pt %g) <-> %d (pt %g)",
               bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
          baseToTagPt[index_bjet] = index_tjet;
        } else {
          LOGF(debug, "DID NOT find pt match: %d (pt %g) <-> %d (pt %g)",
               bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
        }
        index_tjet++;
      }
      index_bjet++;
    }

    index_bjet = 0;
    index_tjet = 0;
    for (const auto& tjet : jetsTagPerColl) {
      for (const auto& bjet : jetsBasePerColl) {
        auto btracksMcParts = bjet.template tracks_as<McParticles>();
        auto btracksTracks = bjet.template tracks_as<Tracks>();
        auto ttracksMcParts = tjet.template tracks_as<McParticles>();
        auto ttracksTracks = tjet.template tracks_as<Tracks>();

        jetsBaseIsMC ? ptSum = getPtSum_SecondArgIsMC(ttracksTracks, btracksMcParts) : ptSum = getPtSum_FirstArgIsMC(ttracksMcParts, btracksTracks);

        if (ptSum > bjet.pt() * minPtFraction) {
          LOGF(debug, "Found pt match: %d (pt %g) <-> %d (pt %g)",
               tjet.globalIndex(), tjet.pt(), bjet.globalIndex(), bjet.pt());
          tagToBasePt[index_tjet] = index_bjet;
        } else {
          LOGF(debug, "DID NOT find pt match: %d (pt %g) <-> %d (pt %g)",
               bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
        }
        index_bjet++;
      }
      index_tjet++;
    }
  }

  // // old pt matching function taken as is from jochen initial matching workflow; it assumed bjet was mcd and tjet was mcp
  // template <typename T, typename U>
  // void MatchPt(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int> & baseToTagPt, std::vector<int> & tagToBasePt)
  // {
  //   int index_bjet = 0;
  //   int index_tjet = 0;
  //   for (const auto& bjet : jetsBasePerColl) {
  //     for (const auto& tjet : jetsTagPerColl) {
  //       float ptSum = 0;
  //       for (const auto& btrack : bjet.template tracks_as<Tracks>()) { // assumes bjet is mcd and tjet is mcp?
  //         for (const auto& ttrack : tjet.template tracks_as<McParticles>()) {
  //           if (btrack.has_mcParticle() &&
  //               ttrack.globalIndex() == btrack.template mcParticle_as<McParticles>().globalIndex()) {
  //             ptSum += ttrack.pt();
  //             break;
  //           }
  //         }
  //       }
  //       if (ptSum > tjet.pt() * minPtFraction) {
  //         LOGF(debug, "Found pt match: %d (pt %g) <-> %d (pt %g)",
  //              bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
  //         baseToTagPt[index_bjet] = index_tjet;
  //       }
  //       index_tjet++;
  //     }
  //     index_bjet++;
  //   }

  //   index_bjet = 0;
  //   index_tjet = 0;
  //   for (const auto& tjet : jetsTagPerColl) {
  //     for (const auto& bjet : jetsBasePerColl) {
  //       float ptSum = 0;
  //       for (const auto& ttrack : tjet.template tracks_as<McParticles>()) {
  //         for (const auto& btrack : bjet.template tracks_as<Tracks>()) { // assumes bjet is mcd and tjet is mcp?s
  //           if (btrack.has_mcParticle() &&
  //               ttrack.globalIndex() == btrack.template mcParticle_as<McParticles>().globalIndex()) {
  //             ptSum += btrack.pt();
  //             break;
  //           }
  //         }
  //       }
  //       if (ptSum > bjet.pt() * minPtFraction) {
  //         LOGF(debug, "Found pt match: %d (pt %g) <-> %d (pt %g)",
  //              tjet.globalIndex(), tjet.pt(), bjet.globalIndex(), bjet.pt());
  //         tagToBasePt[index_tjet] = index_bjet;
  //       }
  //       index_bjet++;
  //     }
  //     index_tjet++;
  //   }
  //   for (std::size_t i = 0; i < baseToTagPt.size(); ++i) {
  //     LOGF(debug, "pt matched base table: bjet %i -> %i", i, baseToTagPt[i]);
  //   }
  //   for (std::size_t i = 0; i < tagToBasePt.size(); ++i) {
  //     LOGF(debug, "pt matched tag table: tjet %i -> %i", i, tagToBasePt[i]);
  //   }
  // }

  // function that calls all the Match functions
  template <typename T, typename U>
  void doAllMatching(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int>& baseToTagGeo, std::vector<int>& baseToTagHF, std::vector<int>& baseToTagPt, std::vector<int>& tagToBaseGeo, std::vector<int>& tagToBaseHF, std::vector<int>& tagToBasePt)
  {
    // geometric matching
    if (doMatchingGeo) {
      MatchGeo(jetsBasePerColl, jetsTagPerColl, baseToTagGeo, tagToBaseGeo);
    }
    // HF matching
    if constexpr (getHfFlag() > 0) {
      if (doMatchingHf) {
        MatchHF(jetsBasePerColl, jetsTagPerColl, baseToTagHF, tagToBaseHF);
      }
    }
    // pt matching
    if (doMatchingPt) {
      MatchPt(jetsBasePerColl, jetsTagPerColl, baseToTagPt, tagToBasePt);
    }
  }

  // function that fills the jetidTagToBase vectors, where the vector that is the i-th entry (corresponding to the tag jet with global index i) sees added to itself the global index of the matched base jet
  template <typename T, typename U>
  void fillJetIdArraysTagToBase(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int> const& tagToBaseGeo, std::vector<int> const& tagToBasePt, std::vector<int> const& tagToBaseHF)
  {
    int geojetidTemp;
    int ptjetidTemp;
    int hfjetidTemp;
    std::vector<int> geojetidTempVector;
    std::vector<int> ptjetidTempVector;
    std::vector<int> hfjetidTempVector;
    for (const auto& jet : jetsTagPerColl) {
      geojetidTemp = tagToBaseGeo[jet.index()];
      if (geojetidTemp > -1 && geojetidTemp < jetsBasePerColl.size()) {
        geojetidTemp = jetsBasePerColl.iteratorAt(geojetidTemp).globalIndex();
        geojetidTagToBase[jet.globalIndex()].push_back(geojetidTemp);
      } else {
        geojetidTemp = -1;
      }
      // Pt matching
      ptjetidTemp = tagToBasePt[jet.index()];
      if (ptjetidTemp > -1 && ptjetidTemp < jetsBasePerColl.size()) {
        ptjetidTemp = jetsBasePerColl.iteratorAt(ptjetidTemp).globalIndex();
        ptjetidTagToBase[jet.globalIndex()].push_back(ptjetidTemp);
      } else {
        ptjetidTemp = -1;
      }
      // HF matching
      hfjetidTemp = tagToBaseHF[jet.index()];
      if (hfjetidTemp > -1 && hfjetidTemp < jetsBasePerColl.size()) {
        hfjetidTemp = jetsBasePerColl.iteratorAt(hfjetidTemp).globalIndex();
        hfjetidTagToBase[jet.globalIndex()].push_back(hfjetidTemp);
      }
      LOGF(debug, "registering matches for tag jet %d (%d): geo -> %d (%d), pT -> %d (%d), HF -> %d",
           jet.index(), jet.globalIndex(), geojetidTemp, tagToBaseGeo[jet.index()], ptjetidTemp, tagToBasePt[jet.index()], tagToBaseHF[jet.index()]);
    }
  }

  // function that fills the jetidBaseToTag vectors, where the vector that is the i-th entry (corresponding to the base jet with global index i) sees added to itself the global index of the matched tag jet
  template <typename T, typename U>
  void fillJetIdArraysBaseToTag(T const& jetsBasePerColl, U const& jetsTagPerColl, std::vector<int> const& baseToTagGeo, std::vector<int> const& baseToTagPt, std::vector<int> const& baseToTagHF)
  {
    int geojetidTemp;
    int ptjetidTemp;
    int hfjetidTemp;
    std::vector<int> geojetidTempVector;
    std::vector<int> ptjetidTempVector;
    std::vector<int> hfjetidTempVector;
    for (const auto& jet : jetsBasePerColl) {
      geojetidTemp = baseToTagGeo[jet.index()];
      if (geojetidTemp > -1 && geojetidTemp < jetsTagPerColl.size()) {
        geojetidTemp = jetsTagPerColl.iteratorAt(geojetidTemp).globalIndex();
        geojetidBaseToTag[jet.globalIndex()].push_back(geojetidTemp);
      } else {
        geojetidTemp = -1;
      }
      // Pt matching
      ptjetidTemp = baseToTagPt[jet.index()];
      if (ptjetidTemp > -1 && ptjetidTemp < jetsTagPerColl.size()) {
        ptjetidTemp = jetsTagPerColl.iteratorAt(ptjetidTemp).globalIndex();
        ptjetidBaseToTag[jet.globalIndex()].push_back(ptjetidTemp);
      }
      // HF matching
      hfjetidTemp = baseToTagHF[jet.index()];
      if (hfjetidTemp > -1 && hfjetidTemp < jetsTagPerColl.size()) {
        hfjetidTemp = jetsTagPerColl.iteratorAt(hfjetidTemp).globalIndex();
        hfjetidBaseToTag[jet.globalIndex()].push_back(hfjetidTemp);
      }
      LOGF(debug, "registering matches for base jet %d (%d): geo -> %d (%d), pT -> %d (%d), HF -> %d",
           jet.index(), jet.globalIndex(), geojetidTemp, baseToTagGeo[jet.index()], ptjetidTemp, baseToTagPt[jet.index()], baseToTagHF[jet.index()]);
    }
  }

  void processDummy(aod::McCollisions const& mcCollisions)
  {
  }
  PROCESS_SWITCH(JetMatching, processDummy, "Dummy process", true);

  // for now:
  // BaseJetCollection and TagJetCollection must have MC info
  void processJets(aod::McCollisions const& mcCollisions, Collisions const& collisions,
                   BaseJetCollection const& jetsBase, TagJetCollection const& jetsTag, McParticles const& particlesMC,
                   Tracks const& tracks, HfCandidates const& hfcandidates)
  {
    // waiting for framework fix to make sliced collection of same type as original collection:
    geojetidBaseToTag.assign(jetsBase.size(), {});
    ptjetidBaseToTag.assign(jetsBase.size(), {});
    hfjetidBaseToTag.assign(jetsBase.size(), {});
    geojetidTagToBase.assign(jetsTag.size(), {});
    ptjetidTagToBase.assign(jetsTag.size(), {});
    hfjetidTagToBase.assign(jetsTag.size(), {});

    for (const auto& mcCollision : mcCollisions) {

      const auto collisionsPerMcColl = collisions.sliceBy(CollisionsCollectionPerMcCollision, mcCollision.globalIndex());

      int i_collision = 0;
      for (const auto& collision : collisionsPerMcColl) {
        ++i_collision;

        const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, jetsBaseIsMC ? mcCollision.globalIndex() : collision.globalIndex());
        const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMC ? mcCollision.globalIndex() : collision.globalIndex());

        // mini jet matching tables for the collection of jets from the current mcCollision (if mcp jet) or collision (if mcd or data jet)
        std::vector<int> baseToTagGeo(jetsBasePerColl.size(), -1), baseToTagHF(jetsBasePerColl.size(), -1), baseToTagPt(jetsBasePerColl.size(), -1);
        std::vector<int> tagToBaseGeo(jetsTagPerColl.size(), -1), tagToBaseHF(jetsTagPerColl.size(), -1), tagToBasePt(jetsTagPerColl.size(), -1);

        LOGF(debug, "performing geometric matching for mcCollision %d and collision %d (%d / %d jets)",
             mcCollision.globalIndex(), collision.globalIndex(), jetsBasePerColl.size(), jetsTagPerColl.size());

        doAllMatching(jetsBasePerColl, jetsTagPerColl, baseToTagGeo, baseToTagHF, baseToTagPt, tagToBaseGeo, tagToBaseHF, tagToBasePt);

        // time to fill the tables
        if (i_collision == 1) { // do this once regardless; if both tag and base are MC, or if one is MC and the other is not and there is no split vertex: it's the only iteration
          fillJetIdArraysBaseToTag(jetsBasePerColl, jetsTagPerColl, baseToTagGeo, baseToTagPt, baseToTagHF);
          fillJetIdArraysTagToBase(jetsBasePerColl, jetsTagPerColl, tagToBaseGeo, tagToBasePt, tagToBaseHF);
        }
        if (i_collision > 1) { // collision is split
          if (!jetsBaseIsMC || !jetsTagIsMC) { // if both are MC, allowing those two lines below to happen would make duplicates in the matching table; this is why thie i_collision>1 case is not treated the same way as the i_collision==1 case
            fillJetIdArraysBaseToTag(jetsBasePerColl, jetsTagPerColl, baseToTagGeo, baseToTagPt, baseToTagHF);
            fillJetIdArraysTagToBase(jetsBasePerColl, jetsTagPerColl, tagToBaseGeo, tagToBasePt, tagToBaseHF);
          }
        }
      }
    }
    for (std::size_t i = 0; i < jetsTag.size(); i++) {
      jetsTagToBase(geojetidTagToBase[i], ptjetidTagToBase[i], hfjetidTagToBase[i]); // is (and needs to) be filled in order
    }
    for (std::size_t i = 0; i < jetsBase.size(); ++i) {
      jetsBaseToTag(geojetidBaseToTag[i], ptjetidBaseToTag[i], hfjetidBaseToTag[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatching, processJets, "Perform jet matching", false);
};

using ChargedJetMatching = JetMatching<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                       soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>,
                                       aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                       aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                                       aod::McParticles,
                                       aod::FwdTracks>;
// using ChargedJetMatching = JetMatching<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>,     // test of the mcd <-> mcp switch in the base/tag categories
//                                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
//                                        aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
//                                        aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
//                                        aod::McParticles,
//                                        aod::DummyTracks>;

using D0ChargedJetMatching = JetMatching<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>,
                                         aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                                         aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                                         soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>,
                                         soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
using LcChargedJetMatching = JetMatching<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,
                                         aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets,
                                         aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets,
                                         soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>,
                                         soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
using BplusChargedJetMatching = JetMatching<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,
                                            soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,
                                            aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets,
                                            aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets,
                                            soa::Join<aod::McParticles, aod::HfCandBplusMcGen>,
                                            soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatching>(cfgc, TaskName{"jet-matching-d0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatching>(cfgc, TaskName{"jet-matching-lc-ch"}));
  tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatching>(cfgc, TaskName{"jet-matching-bplus-ch"}));

  return WorkflowSpec{tasks};
}
