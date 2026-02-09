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

#ifndef PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHING_H_
#define PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHING_H_

#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h> // IWYU pragma: export

#include <vector>

template <typename JetsBase, typename JetsTag, typename SplittingsBasetoTagMatchingTable, typename SplittingsTagtoBaseMatchingTable, typename PairsBasetoTagMatchingTable, typename PairsTagtoBaseMatchingTable, typename SplittingsBase, typename SplittingsTag, typename PairsBase, typename PairsTag, typename CandidatesBase, typename CandidatesTag, typename TracksBase, typename TracksTag, typename ClustersBase>
struct JetSubstructureMatching {

  o2::framework::Produces<SplittingsBasetoTagMatchingTable> splittingsBasetoTagMatchingTable;
  o2::framework::Produces<SplittingsTagtoBaseMatchingTable> splittingsTagtoBaseMatchingTable;

  o2::framework::Produces<PairsBasetoTagMatchingTable> pairsBasetoTagMatchingTable;
  o2::framework::Produces<PairsTagtoBaseMatchingTable> pairsTagtoBaseMatchingTable;

  o2::framework::Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  o2::framework::Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  o2::framework::Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  o2::framework::Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  o2::framework::Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};
  o2::framework::Configurable<bool> requireGeoMatchedJets{"requireGeoMatchedJets", false, "require jets are geo matched as well"};
  o2::framework::Configurable<bool> requirePtMatchedJets{"requirePtMatchedJets", false, "require jets are pT matched as well"};
  o2::framework::Configurable<bool> requireHFMatchedJets{"requireHFMatchedJets", false, "require jets are HF matched as well"};

  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<o2::aod::JetMcCollisions, JetsTag>();

  void init(o2::framework::InitContext const&)
  {
  }

  o2::framework::PresliceOptional<o2::aod::ChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetInclusive = o2::aod::chargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::ChargedMCParticleLevelSPs> TagSplittingsPerTagJetInclusive = o2::aod::chargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetD0 = o2::aod::d0chargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedMCParticleLevelSPs> TagSplittingsPerTagJetD0 = o2::aod::d0chargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDplus = o2::aod::dpluschargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedMCParticleLevelSPs> TagSplittingsPerTagJetDplus = o2::aod::dpluschargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDs = o2::aod::dschargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedMCParticleLevelSPs> TagSplittingsPerTagJetDs = o2::aod::dschargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDstar = o2::aod::dstarchargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedMCParticleLevelSPs> TagSplittingsPerTagJetDstar = o2::aod::dstarchargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetLc = o2::aod::lcchargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedMCParticleLevelSPs> TagSplittingsPerTagJetLc = o2::aod::lcchargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetB0 = o2::aod::b0chargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedMCParticleLevelSPs> TagSplittingsPerTagJetB0 = o2::aod::b0chargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetBplus = o2::aod::bpluschargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedMCParticleLevelSPs> TagSplittingsPerTagJetBplus = o2::aod::bpluschargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetXicToXiPiPi = o2::aod::xictoxipipichargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedMCParticleLevelSPs> TagSplittingsPerTagJetXicToXiPiPi = o2::aod::xictoxipipichargedmcparticlelevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedMCDetectorLevelSPs> BaseSplittingsPerBaseJetDielectron = o2::aod::dielectronchargedmcdetectorlevelsplitting::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedMCParticleLevelSPs> TagSplittingsPerTagJetDielectron = o2::aod::dielectronchargedmcparticlelevelsplitting::jetId;

  o2::framework::PresliceOptional<o2::aod::ChargedMCDetectorLevelPRs> BasePairsPerBaseJetInclusive = o2::aod::chargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::ChargedMCParticleLevelPRs> TagPairsPerTagJetInclusive = o2::aod::chargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedMCDetectorLevelPRs> BasePairsPerBaseJetD0 = o2::aod::d0chargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::D0ChargedMCParticleLevelPRs> TagPairsPerTagJetD0 = o2::aod::d0chargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedMCDetectorLevelPRs> BasePairsPerBaseJetDplus = o2::aod::dpluschargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DplusChargedMCParticleLevelPRs> TagPairsPerTagJetDplus = o2::aod::dpluschargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedMCDetectorLevelPRs> BasePairsPerBaseJetDs = o2::aod::dschargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DsChargedMCParticleLevelPRs> TagPairsPerTagJetDs = o2::aod::dschargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedMCDetectorLevelPRs> BasePairsPerBaseJetDstar = o2::aod::dstarchargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DstarChargedMCParticleLevelPRs> TagPairsPerTagJetDstar = o2::aod::dstarchargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedMCDetectorLevelPRs> BasePairsPerBaseJetLc = o2::aod::lcchargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::LcChargedMCParticleLevelPRs> TagPairsPerTagJetLc = o2::aod::lcchargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedMCDetectorLevelPRs> BasePairsPerBaseJetB0 = o2::aod::b0chargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::B0ChargedMCParticleLevelPRs> TagPairsPerTagJetB0 = o2::aod::b0chargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedMCDetectorLevelPRs> BasePairsPerBaseJetBplus = o2::aod::bpluschargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::BplusChargedMCParticleLevelPRs> TagPairsPerTagJetBplus = o2::aod::bpluschargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedMCDetectorLevelPRs> BasePairsPerBaseJetXicToXiPiPi = o2::aod::xictoxipipichargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::XicToXiPiPiChargedMCParticleLevelPRs> TagPairsPerTagJetXicToXiPiPi = o2::aod::xictoxipipichargedmcparticlelevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedMCDetectorLevelPRs> BasePairsPerBaseJetDielectron = o2::aod::dielectronchargedmcdetectorlevelpair::jetId;
  o2::framework::PresliceOptional<o2::aod::DielectronChargedMCParticleLevelPRs> TagPairsPerTagJetDielectron = o2::aod::dielectronchargedmcparticlelevelpair::jetId;

  // workaround till binding nodes can be passed as template arguments
  template <typename CandidateTable, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R, typename S, typename A, typename B>
  auto slicedPerJetForMatching(T const& table, U const& jet, V const& perInclusiveJet, M const& perD0Jet, N const& perDplusJet, O const& perDsJet, P const& perDstarJet, Q const& perLcJet, R const& perB0Jet, S const& perBplusJet, A const& perXicToXiPiPiJet, B const& perDielectronJet)
  {
    if constexpr (jethfutilities::isHFTable<CandidateTable>() || jethfutilities::isHFMcTable<CandidateTable>()) {
      return jethfutilities::slicedPerHFJet<CandidateTable>(table, jet, perD0Jet, perDplusJet, perDsJet, perDstarJet, perLcJet, perB0Jet, perBplusJet, perXicToXiPiPiJet);
    } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>() || jetdqutilities::isDielectronMcTable<CandidateTable>()) {
      return jetdqutilities::slicedPerDielectronJet<CandidateTable>(table, jet, perDielectronJet);
    } else {
      return table.sliceBy(perInclusiveJet, jet.globalIndex());
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
        auto const& jetTagSplittings = slicedPerJetForMatching<CandidatesTag>(jetsTagSplittings, jetTag, TagSplittingsPerTagJetInclusive, TagSplittingsPerTagJetD0, TagSplittingsPerTagJetDplus, TagSplittingsPerTagJetDs, TagSplittingsPerTagJetDstar, TagSplittingsPerTagJetLc, TagSplittingsPerTagJetB0, TagSplittingsPerTagJetBplus, TagSplittingsPerTagJetXicToXiPiPi, TagSplittingsPerTagJetDielectron);
        int tagSplittingIndex = 0;
        for (auto const& jetTagSplitting : jetTagSplittings) {
          jetTagSplittingsMap[jetTagSplitting.globalIndex()] = tagSplittingIndex;
          tagSplittingIndex++;
        }
        // auto const& jetTagPairs = jetsTagPairs.sliceBy(TagPairsPerTagJet, jetTag.globalIndex());
        auto const& jetTagPairs = slicedPerJetForMatching<CandidatesTag>(jetsTagPairs, jetTag, TagPairsPerTagJetInclusive, TagPairsPerTagJetD0, TagPairsPerTagJetDplus, TagPairsPerTagJetDs, TagPairsPerTagJetDstar, TagPairsPerTagJetLc, TagPairsPerTagJetB0, TagPairsPerTagJetBplus, TagPairsPerTagJetXicToXiPiPi, TagPairsPerTagJetDielectron);
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
          auto const& jetBaseSplittings = slicedPerJetForMatching<CandidatesBase>(jetsBaseSplittings, jetBase, BaseSplittingsPerBaseJetInclusive, BaseSplittingsPerBaseJetD0, BaseSplittingsPerBaseJetDplus, BaseSplittingsPerBaseJetDs, BaseSplittingsPerBaseJetDstar, BaseSplittingsPerBaseJetLc, BaseSplittingsPerBaseJetB0, BaseSplittingsPerBaseJetBplus, BaseSplittingsPerBaseJetXicToXiPiPi, BaseSplittingsPerBaseJetDielectron);
          int baseSplittingIndex = 0;
          for (auto const& jetBaseSplitting : jetBaseSplittings) {
            jetBaseSplittingsMap[jetBaseSplitting.globalIndex()] = baseSplittingIndex;
            baseSplittingIndex++;
          }
          jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetBaseSplittings, jetTagSplittings, jetsBasetoTagSplittingsMatchingGeo, jetsBasetoTagSplittingsMatchingPt, jetsBasetoTagSplittingsMatchingHF, jetsTagtoBaseSplittingsMatchingGeo, jetsTagtoBaseSplittingsMatchingPt, jetsTagtoBaseSplittingsMatchingHF, candidatesBase, tracksBase, clustersBase, candidatesTag, tracksTag, tracksTag, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
          // auto const& jetBasePairs = jetsBasePairs.sliceBy(BasePairsPerBaseJet, jetBase.globalIndex());
          auto const& jetBasePairs = slicedPerJetForMatching<CandidatesBase>(jetsBasePairs, jetBase, BasePairsPerBaseJetInclusive, BasePairsPerBaseJetD0, BasePairsPerBaseJetDplus, BasePairsPerBaseJetDs, BasePairsPerBaseJetDstar, BasePairsPerBaseJetLc, BasePairsPerBaseJetB0, BasePairsPerBaseJetBplus, BasePairsPerBaseJetXicToXiPiPi, BasePairsPerBaseJetDielectron);
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

#endif // PWGJE_TABLEPRODUCER_MATCHING_SUBSTRUCTURE_JETSUBSTRUCTUREMATCHING_H_
