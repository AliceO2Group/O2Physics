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

// jet analysis tasks (subscribing to jet finder task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetMatchingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

template <typename JetsBase, typename JetsTag, typename SplittingsBasetoTagMatchingTable, typename SplittingsTagtoBaseMatchingTable, typename PairsBasetoTagMatchingTable, typename PairsTagtoBaseMatchingTable, typename SplittingsBase, typename SplittingsTag, typename PairsBase, typename PairsTag, typename CandidatesBase, typename CandidatesTag, typename TracksBase, typename TracksTag, typename ClustersBase>
struct JetSubstructureMatching {

  Produces<SplittingsBasetoTagMatchingTable> splittingsBasetoTagMatchingTable;
  Produces<SplittingsTagtoBaseMatchingTable> splittingsTagtoBaseMatchingTable;

  Produces<PairsBasetoTagMatchingTable> pairsBasetoTagMatchingTable;
  Produces<PairsTagtoBaseMatchingTable> pairsTagtoBaseMatchingTable;

  Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};
  Configurable<bool> requireGeoMatchedJets{"requireGeoMatchedJets", false, "require jets are geo matched as well"};
  Configurable<bool> requirePtMatchedJets{"requirePtMatchedJets", false, "require jets are pT matched as well"};
  Configurable<bool> requireHFMatchedJets{"requireHFMatchedJets", false, "require jets are HF matched as well"};

  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<aod::JetMcCollisions, JetsTag>();

  void init(InitContext const&)
  {
  }

  PresliceOptional<aod::ChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetInclusive = aod::chargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::ChargedMCParticleLevelSPs> TagSplittingsPerTagJetInclusive = aod::chargedmcparticlelevelsplitting::jetId;
  PresliceOptional<aod::D0ChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetD0 = aod::d0chargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::D0ChargedMCParticleLevelSPs> TagSplittingsPerTagJetD0 = aod::d0chargedmcparticlelevelsplitting::jetId;
  PresliceOptional<aod::DplusChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDplus = aod::dpluschargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::DplusChargedMCParticleLevelSPs> TagSplittingsPerTagJetDplus = aod::dpluschargedmcparticlelevelsplitting::jetId;
  PresliceOptional<aod::LcChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetLc = aod::lcchargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::LcChargedMCParticleLevelSPs> TagSplittingsPerTagJetLc = aod::lcchargedmcparticlelevelsplitting::jetId;
  PresliceOptional<aod::BplusChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetBplus = aod::bpluschargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::BplusChargedMCParticleLevelSPs> TagSplittingsPerTagJetBplus = aod::bpluschargedmcparticlelevelsplitting::jetId;
  PresliceOptional<aod::DielectronChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDielectron = aod::dielectronchargedmcdetectorlevelsplitting::jetId;
  PresliceOptional<aod::DielectronChargedMCParticleLevelSPs> TagSplittingsPerTagJetDielectron = aod::dielectronchargedmcparticlelevelsplitting::jetId;

  PresliceOptional<aod::ChargedMCDetectorLevelPRs> BasePairsPerBaseJetInclusive = aod::chargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::ChargedMCParticleLevelPRs> TagPairsPerTagJetInclusive = aod::chargedmcparticlelevelpair::jetId;
  PresliceOptional<aod::D0ChargedMCDetectorLevelPRs> BasePairsPerBaseJetD0 = aod::d0chargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::D0ChargedMCParticleLevelPRs> TagPairsPerTagJetD0 = aod::d0chargedmcparticlelevelpair::jetId;
  PresliceOptional<aod::DplusChargedMCDetectorLevelPRs> BasePairsPerBaseJetDplus = aod::dpluschargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::DplusChargedMCParticleLevelPRs> TagPairsPerTagJetDplus = aod::dpluschargedmcparticlelevelpair::jetId;
  PresliceOptional<aod::LcChargedMCDetectorLevelPRs> BasePairsPerBaseJetLc = aod::lcchargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::LcChargedMCParticleLevelPRs> TagPairsPerTagJetLc = aod::lcchargedmcparticlelevelpair::jetId;
  PresliceOptional<aod::BplusChargedMCDetectorLevelPRs> BasePairsPerBaseJetBplus = aod::bpluschargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::BplusChargedMCParticleLevelPRs> TagPairsPerTagJetBplus = aod::bpluschargedmcparticlelevelpair::jetId;
  PresliceOptional<aod::DielectronChargedMCDetectorLevelPRs> BasePairsPerBaseJetDielectron = aod::dielectronchargedmcdetectorlevelpair::jetId;
  PresliceOptional<aod::DielectronChargedMCParticleLevelPRs> TagPairsPerTagJetDielectron = aod::dielectronchargedmcparticlelevelpair::jetId;

  // workaround till binding nodes can be passed as template arguments
  template <typename CandidateTable, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q>
  auto slicedPerJetForMatching(T const& table, U const& jet, V const& perIncluisveJet, M const& perD0Jet, N const& perDplusJet, O const& perLcJet, P const& perBplusJet, Q const& perDielectronJet)
  {
    if constexpr (jethfutilities::isHFTable<CandidateTable>() || jethfutilities::isHFMcTable<CandidateTable>()) {
      return jethfutilities::slicedPerHFJet<CandidateTable>(table, jet, perD0Jet, perDplusJet, perLcJet, perBplusJet);
    } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>() || jetdqutilities::isDielectronMcTable<CandidateTable>()) {
      return jetdqutilities::slicedPerDielectronJet<CandidateTable>(table, jet, perDielectronJet);
    } else {
      return table.sliceBy(perIncluisveJet, jet.globalIndex());
    }
  }

  template <typename T>
  auto defaultMatchedJets(T const& jetTag)
  {
    if constexpr (jetcandidateutilities::isCandidateTable<CandidatesBase>()) {
      return jetTag.template matchedJetCand_as<JetsBase>();
    } else {
      return jetTag.template matchedJetGeo_as<JetsBase>();
    }
  }

  void processData(JetsTag const& jetsTag,
                   JetsBase const&,
                   SplittingsBase const& jetsBaseSplittings,
                   SplittingsTag const& jetsTagSplittings,
                   PairsBase const& jetsBasePairs,
                   PairsTag const& jetsTagPairs,
                   CandidatesBase const& candidatesBase,
                   CandidatesTag const& candidatesTag,
                   ClustersBase const& clustersBase,
                   TracksBase const& tracksBase, TracksTag const& tracksTag)
  {
    std::vector<std::vector<int>> jetsBasetoTagSplittingsMatchingGeo, jetsBasetoTagSplittingsMatchingPt, jetsBasetoTagSplittingsMatchingHF;
    jetsBasetoTagSplittingsMatchingGeo.assign(jetsBaseSplittings.size(), {});
    jetsBasetoTagSplittingsMatchingPt.assign(jetsBaseSplittings.size(), {});
    jetsBasetoTagSplittingsMatchingHF.assign(jetsBaseSplittings.size(), {});

    std::vector<std::vector<int>> jetsTagtoBaseSplittingsMatchingGeo, jetsTagtoBaseSplittingsMatchingPt, jetsTagtoBaseSplittingsMatchingHF;
    jetsTagtoBaseSplittingsMatchingGeo.assign(jetsTagSplittings.size(), {});
    jetsTagtoBaseSplittingsMatchingPt.assign(jetsTagSplittings.size(), {});
    jetsTagtoBaseSplittingsMatchingHF.assign(jetsTagSplittings.size(), {});

    std::vector<int> jetTagSplittingsMap;
    std::vector<int> jetBaseSplittingsMap;
    jetTagSplittingsMap.resize(jetsTagSplittings.size(), -1);
    jetBaseSplittingsMap.resize(jetsBaseSplittings.size(), -1);

    std::vector<std::vector<int>> jetsBasetoTagPairsMatching;
    std::vector<std::vector<int>> jetsTagtoBasePairsMatching;
    jetsBasetoTagPairsMatching.assign(jetsBasePairs.size(), {});
    jetsTagtoBasePairsMatching.assign(jetsTagPairs.size(), {});

    std::vector<int> jetTagPairsMap;
    std::vector<int> jetBasePairsMap;
    jetTagPairsMap.resize(jetsTagPairs.size(), -1);
    jetBasePairsMap.resize(jetsBasePairs.size(), -1);

    for (auto jetTag : jetsTag) {
      bool hasMatchedJet = false;
      if constexpr (jetcandidateutilities::isCandidateTable<CandidatesBase>()) {
        hasMatchedJet = jetTag.has_matchedJetCand();
      } else {
        hasMatchedJet = jetTag.has_matchedJetGeo();
      }
      if (hasMatchedJet) {
        // auto const& jetTagSplittings = jetsTagSplittings.sliceBy(TagSplittingsPerTagJet, jetTag.globalIndex());
        auto const& jetTagSplittings = slicedPerJetForMatching<CandidatesTag>(jetsTagSplittings, jetTag, TagSplittingsPerTagJetInclusive, TagSplittingsPerTagJetD0, TagSplittingsPerTagJetDplus, TagSplittingsPerTagJetLc, TagSplittingsPerTagJetBplus, TagSplittingsPerTagJetDielectron);
        int tagSplittingIndex = 0;
        for (auto const& jetTagSplitting : jetTagSplittings) {
          jetTagSplittingsMap[jetTagSplitting.globalIndex()] = tagSplittingIndex;
          tagSplittingIndex++;
        }
        // auto const& jetTagPairs = jetsTagPairs.sliceBy(TagPairsPerTagJet, jetTag.globalIndex());
        auto const& jetTagPairs = slicedPerJetForMatching<CandidatesTag>(jetsTagPairs, jetTag, TagPairsPerTagJetInclusive, TagPairsPerTagJetD0, TagPairsPerTagJetDplus, TagPairsPerTagJetLc, TagPairsPerTagJetBplus, TagPairsPerTagJetDielectron);
        int tagPairIndex = 0;
        for (auto const& jetTagPair : jetTagPairs) {
          jetTagPairsMap[jetTagPair.globalIndex()] = tagPairIndex;
          tagPairIndex++;
        }
        for (auto& jetBase : defaultMatchedJets(jetTag)) {
          if (requireGeoMatchedJets) {
            bool jetsMatchedWithGeo = false;
            for (auto& jetBaseForMatchGeo : jetTag.template matchedJetGeo_as<JetsBase>()) {
              if (jetBaseForMatchGeo.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithGeo = true;
              }
            }
            if (!jetsMatchedWithGeo) {
              continue;
            }
          }
          if (requirePtMatchedJets) {
            bool jetsMatchedWithPt = false;
            for (auto& jetBaseForMatchPt : jetTag.template matchedJetPt_as<JetsBase>()) {
              if (jetBaseForMatchPt.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithPt = true;
              }
            }
            if (!jetsMatchedWithPt) {
              continue;
            }
          }
          if (requireHFMatchedJets) {
            bool jetsMatchedWithHF = false;
            for (auto& jetBaseForMatchHF : jetTag.template matchedJetCand_as<JetsBase>()) {
              if (jetBaseForMatchHF.globalIndex() == jetBase.globalIndex()) {
                jetsMatchedWithHF = true;
              }
            }
            if (!jetsMatchedWithHF) {
              continue;
            }
          }
          // auto const& jetBaseSplittings = jetsBaseSplittings.sliceBy(BaseSplittingsPerBaseJet, jetBase.globalIndex());
          auto const& jetBaseSplittings = slicedPerJetForMatching<CandidatesBase>(jetsBaseSplittings, jetBase, BaseSplittingsPerBaseJetInclusive, BaseSplittingsPerBaseJetD0, BaseSplittingsPerBaseJetDplus, BaseSplittingsPerBaseJetLc, BaseSplittingsPerBaseJetBplus, BaseSplittingsPerBaseJetDielectron);
          int baseSplittingIndex = 0;
          for (auto const& jetBaseSplitting : jetBaseSplittings) {
            jetBaseSplittingsMap[jetBaseSplitting.globalIndex()] = baseSplittingIndex;
            baseSplittingIndex++;
          }
          jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetBaseSplittings, jetTagSplittings, jetsBasetoTagSplittingsMatchingGeo, jetsBasetoTagSplittingsMatchingPt, jetsBasetoTagSplittingsMatchingHF, jetsTagtoBaseSplittingsMatchingGeo, jetsTagtoBaseSplittingsMatchingPt, jetsTagtoBaseSplittingsMatchingHF, candidatesBase, tracksBase, clustersBase, candidatesTag, tracksTag, tracksTag, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
          // auto const& jetBasePairs = jetsBasePairs.sliceBy(BasePairsPerBaseJet, jetBase.globalIndex());
          auto const& jetBasePairs = slicedPerJetForMatching<CandidatesBase>(jetsBasePairs, jetBase, BasePairsPerBaseJetInclusive, BasePairsPerBaseJetD0, BasePairsPerBaseJetDplus, BasePairsPerBaseJetLc, BasePairsPerBaseJetBplus, BasePairsPerBaseJetDielectron);
          int basePairIndex = 0;
          for (auto const& jetBasePair : jetBasePairs) {
            jetBasePairsMap[jetBasePair.globalIndex()] = basePairIndex;
            basePairIndex++;
          }
          jetmatchingutilities::doPairMatching<jetsBaseIsMc, jetsTagIsMc>(jetBasePairs, jetTagPairs, jetsBasetoTagPairsMatching, jetsTagtoBasePairsMatching, candidatesBase, tracksBase, candidatesTag, tracksTag);
        }
      }
    }
    for (auto jetsTagSplitting : jetsTagSplittings) {
      std::vector<int> tagToBaseMatchingGeoIndex;
      std::vector<int> tagToBaseMatchingPtIndex;
      std::vector<int> tagToBaseMatchingHFIndex;
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingGeo[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingGeoIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingPt[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingPtIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      for (auto jetBaseSplittingIndex : jetsTagtoBaseSplittingsMatchingHF[jetsTagSplitting.globalIndex()]) {
        tagToBaseMatchingHFIndex.push_back(jetBaseSplittingsMap[jetBaseSplittingIndex]);
      }
      splittingsTagtoBaseMatchingTable(tagToBaseMatchingGeoIndex, tagToBaseMatchingPtIndex, tagToBaseMatchingHFIndex);
    }
    for (auto jetsBaseSplitting : jetsBaseSplittings) {
      std::vector<int> baseToTagMatchingGeoIndex;
      std::vector<int> baseToTagMatchingPtIndex;
      std::vector<int> baseToTagMatchingHFIndex;
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingGeo[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingGeoIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingPt[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingPtIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      for (auto jetTagSplittingIndex : jetsBasetoTagSplittingsMatchingHF[jetsBaseSplitting.globalIndex()]) {
        baseToTagMatchingHFIndex.push_back(jetTagSplittingsMap[jetTagSplittingIndex]);
      }
      splittingsBasetoTagMatchingTable(baseToTagMatchingGeoIndex, baseToTagMatchingPtIndex, baseToTagMatchingHFIndex);
    }
    for (auto jetsTagPair : jetsTagPairs) {
      std::vector<int> tagToBaseMatchingIndex;
      for (auto jetBasePairIndex : jetsTagtoBasePairsMatching[jetsTagPair.globalIndex()]) {
        tagToBaseMatchingIndex.push_back(jetBasePairsMap[jetBasePairIndex]);
      }
      pairsTagtoBaseMatchingTable(tagToBaseMatchingIndex);
    }
    for (auto jetsBasePair : jetsBasePairs) {
      std::vector<int> baseToTagMatchingIndex;
      for (auto jetTagPairIndex : jetsBasetoTagPairsMatching[jetsBasePair.globalIndex()]) {
        baseToTagMatchingIndex.push_back(jetTagPairsMap[jetTagPairIndex]);
      }
      pairsBasetoTagMatchingTable(baseToTagMatchingIndex);
    }
  }
  PROCESS_SWITCH(JetSubstructureMatching, processData, "charged jet substructure", true);
};
