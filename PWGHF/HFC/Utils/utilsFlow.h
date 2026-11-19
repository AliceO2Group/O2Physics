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

/// \file utilsFlow.h
/// \brief Utilities for HFC analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>

#ifndef PWGHF_HFC_UTILS_UTILSFLOW_H_
#define PWGHF_HFC_UTILS_UTILSFLOW_H_

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <csignal>
#include <cstdlib>
#include <utility>
#include <vector>

namespace o2::hfc_utils_flow
{

enum TriggerChannels {
  DplusToPiKPi = 0,
  DsToKKPi,
  DsToPiKK,
  D0ToPiK,
  D0ToKPi,
  LcToPKPi,
  LcToPiKP,
  Hadron
};

enum AssTracksTypes {
  Pi = 0,
  Ka,
  Pr,
  El,
  Mu,
  NAssTracksTypes
};

class HfcFlowHelper
{

 public:
  HfcFlowHelper() = default;

  /// Initialize the HFC flow helper
  /// \param ccdbService is the CCDB service
  /// \param histRegistry is the histogram registry
  /// \param mlClasses are the machine learning classes
  /// \param maxRapidity is the maximum rapidity
  void init(std::vector<unsigned int> mlClasses,
            float maxRapidity)
  {
    this->classMl = std::move(mlClasses);
    this->yCandMax = maxRapidity;
    LOG(info) << "HFC flow helper initialized";
  }

  /// Get charm hadron candidate mass
  /// \param candidate is the charm hadron candidate
  template <TriggerChannels CandType, typename TCand>
  double getCandMass(const TCand& candidate)
  {
    if constexpr (CandType == TriggerChannels::DsToKKPi) {
      return HfHelper::invMassDsToKKPi(candidate);
    }
    if constexpr (CandType == TriggerChannels::DsToPiKK) {
      return HfHelper::invMassDsToPiKK(candidate);
    }
    if constexpr (CandType == TriggerChannels::DplusToPiKPi) {
      return HfHelper::invMassDplusToPiKPi(candidate);
    }
    if constexpr (CandType == TriggerChannels::D0ToPiK) {
      return HfHelper::invMassD0ToPiK(candidate);
    }
    if constexpr (CandType == TriggerChannels::D0ToKPi) {
      return HfHelper::invMassD0barToKPi(candidate);
    }
    if constexpr (CandType == TriggerChannels::LcToPKPi) {
      return HfHelper::invMassLcToPKPi(candidate);
    }
    if constexpr (CandType == TriggerChannels::LcToPiKP) {
      return HfHelper::invMassLcToPiKP(candidate);
    }
    return -1.;
  }

  /// Get charm hadron bdt scores
  /// \param candidate is the charm hadron candidate
  template <TriggerChannels CandType, typename TCand>
  std::array<float, 2> getCandMlScores(const TCand& candidate)
  {
    std::array<float, 2> outputMl{-999.f, -999.f};
    if constexpr (CandType == TriggerChannels::DsToKKPi) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::DsToPiKK) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::DplusToPiKPi) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::D0ToPiK) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::D0ToKPi) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0bar()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::LcToPKPi) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl.at(iclass)];
      }
    }
    if constexpr (CandType == TriggerChannels::LcToPiKP) {
      for (unsigned int iclass = 0; iclass < classMl.size(); iclass++) {
        outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl.at(iclass)];
      }
    }
    return outputMl;
  }

  /// Cut on rapidity of the candidate
  /// \param candidate is the charm hadron candidate
  /// \return true if the candidate passes the rapidity cut, false otherwise
  template <TriggerChannels CandType, typename TCand>
  bool cutCandRapidity(const TCand& candidate)
  {
    double y = 0.0;
    if constexpr (CandType == TriggerChannels::DsToKKPi || CandType == TriggerChannels::DsToPiKK) {
      y = HfHelper::yDs(candidate);
    } else if constexpr (CandType == TriggerChannels::DplusToPiKPi) {
      y = HfHelper::yDplus(candidate);
    } else if constexpr (CandType == TriggerChannels::D0ToPiK || CandType == TriggerChannels::D0ToKPi) {
      y = HfHelper::yD0(candidate);
    } else if constexpr (CandType == TriggerChannels::LcToPKPi || CandType == TriggerChannels::LcToPiKP) {
      y = HfHelper::yLc(candidate);
    } else {
      return true;
    }
    return std::fabs(y) < yCandMax;
  }

  /// Runtime dispatch to the compile-time checkCollision
  /// \param collision is the collision
  /// \param centEst is the centrality estimator
  /// \param cent is the centrality value
  /// \param mult is the multiplicity value
  template <typename BCsInfo, typename GenColl, typename RecoCollSlice>
  bool selectGenCollision(o2::hf_evsel::HfEventSelectionMc& hfEvSelMc,
                          GenColl const& genColl,
                          RecoCollSlice const& recoCollsSlice,
                          o2::hf_centrality::CentralityEstimator centEst,
                          float& cent,
                          float& mult)
  {
    int nSplitColl = recoCollsSlice.size();
    o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
    switch (centEst) {
      case o2::hf_centrality::CentralityEstimator::FT0A:
        collRejMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::FT0A>(genColl, recoCollsSlice, cent);
        mult = o2::hf_centrality::getCentralityGenColl(recoCollsSlice, o2::hf_centrality::CentralityEstimator::FT0A);
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::FT0A>(genColl, collRejMask, nSplitColl);
        break;
      case o2::hf_centrality::CentralityEstimator::FT0C:
        collRejMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::FT0C>(genColl, recoCollsSlice, cent);
        mult = o2::hf_centrality::getCentralityGenColl(recoCollsSlice, o2::hf_centrality::CentralityEstimator::FT0C);
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::FT0C>(genColl, collRejMask, nSplitColl);
        break;
      case o2::hf_centrality::CentralityEstimator::FT0M:
        collRejMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::FT0M>(genColl, recoCollsSlice, cent);
        mult = o2::hf_centrality::getCentralityGenColl(recoCollsSlice, o2::hf_centrality::CentralityEstimator::FT0M);
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::FT0M>(genColl, collRejMask, nSplitColl);
        break;
      case o2::hf_centrality::CentralityEstimator::FV0A:
        collRejMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::FV0A>(genColl, recoCollsSlice, cent);
        mult = o2::hf_centrality::getCentralityGenColl(recoCollsSlice, o2::hf_centrality::CentralityEstimator::FV0A);
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::FV0A>(genColl, collRejMask, nSplitColl);
        break;
      case o2::hf_centrality::CentralityEstimator::None:
        collRejMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::None>(genColl, recoCollsSlice, cent);
        mult = -1.f;
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::None>(genColl, collRejMask, nSplitColl);
        break;
      default:
        LOG(fatal) << "Centrality estimator not recognized for collision selection";
        std::abort();
    }

    return (collRejMask == 0);
  }

  /// Apply event selection on reconstructed collisions
  /// \param collision is the collision
  /// \param centEst is the centrality estimator
  /// \param cent is the centrality value
  /// \param mult is the multiplicity value
  template <typename Coll>
  bool selectRecoCollision(o2::hf_evsel::HfEventSelection& hfEvSel,
                           o2::framework::Service<o2::ccdb::BasicCCDBManager>& ccdb,
                           o2::framework::HistogramRegistry& registry,
                           Coll const& collision,
                           o2::hf_centrality::CentralityEstimator centEst,
                           float& cent,
                           float& mult)
  {
    o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
    switch (centEst) {
      case o2::hf_centrality::CentralityEstimator::FT0A:
        collRejMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::FT0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        mult = collision.multFT0A();
        break;
      case o2::hf_centrality::CentralityEstimator::FT0C:
        collRejMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        mult = collision.multFT0C();
        break;
      case o2::hf_centrality::CentralityEstimator::FT0M:
        collRejMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        mult = collision.multFT0M();
        break;
      case o2::hf_centrality::CentralityEstimator::FV0A:
        collRejMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::FV0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        mult = collision.multFV0A();
        break;
      case o2::hf_centrality::CentralityEstimator::None:
        collRejMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
        mult = -1.f;
        break;
      default:
        LOG(fatal) << "Centrality estimator not recognized for collision selection";
        std::abort();
    }

    hfEvSel.fillHistograms(collision, collRejMask, cent);

    return (collRejMask == 0);
  }

 private:
  std::vector<unsigned int> classMl;
  float yCandMax{-999.f};
};

} // namespace o2::hfc_utils_flow
#endif // PWGHF_HFC_UTILS_UTILSFLOW_H_
