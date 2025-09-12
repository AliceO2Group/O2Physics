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

/// \file derivedDataCreatorCorrelationsReduced.cxx
/// \brief CharmHadrons-Hadrons correlator tree creator for data and MC-reco analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, CERN
/// \author Wu Chuntai <chuntai.wu@cern.ch>, CCNU, INFN Padova, and Padova University

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum DecayChannel {
  DplusToPiKPi = 0,
  DsToKKPi,
  DsToPiKK,
  D0ToPiK,
  D0ToKPi
};

/// Code to select collisions with at least one Ds meson
struct HfDerivedDataCreatorCorrelationsReduced {
  Produces<aod::HfcRedFlowColls> rowCollisions;
  Produces<aod::HfcRedCharmTrigs> rowCharmCandidates;
  Produces<aod::HfcRedCharmMls> rowCharmCandidatesMl;
  Produces<aod::HfcRedTrkAssocs> rowAssocTrackReduced;
  Produces<aod::HfcRedTrkSels> rowAssocTrackSelInfo;

  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<bool> forceCharmInCollision{"forceCharmInCollision", false, "Flag to force charm in collision"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 24., "max. cand. pT"};
  Configurable<float> etaTrackMax{"etaTrackMax", 1., "max. track eta"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.15, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 5., "max. track pT"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. track DCA Z"};

  HfHelper hfHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  SliceCache cache;

  double massCharm{0.};

  using CollsWithCentMult = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
  using CandDsDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectTrackData = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  Preslice<CandDsDataWMl> candsDsPerCollWMl = aod::hf_cand::collisionId;
  Preslice<CandDplusDataWMl> candsDplusPerCollWMl = aod::hf_cand::collisionId;
  Preslice<CandD0DataWMl> candsD0PerCollWMl = aod::hf_cand::collisionId;
  Preslice<TracksData> trackIndicesPerColl = aod::track::collisionId;

  Partition<CandDsDataWMl> selectedDsToKKPiWMl = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsDataWMl> selectedDsToPiKKWMl = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (doprocessDplusWithMl) {
      massCharm = o2::constants::physics::MassDPlus;
    } else if (doprocessDsWithMl) {
      massCharm = o2::constants::physics::MassDS;
    } else if (doprocessD0WithMl) {
      massCharm = o2::constants::physics::MassD0;
    } else {
      LOG(fatal) << "No decay channel selected to process";
    }

    hfEvSel.addHistograms(registry); // collision monitoring
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Check event selections for collision and fill the collision table
  /// \param collision is the collision
  template <typename Coll>
  bool checkAndFillCollision(Coll const& collision)
  {
    float cent{-1.f};
    float mult{-1.f};
    o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
    if (centEstimator == CentralityEstimator::FT0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0A();
    } else if (centEstimator == CentralityEstimator::FT0C) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0C();
    } else if (centEstimator == CentralityEstimator::FT0M) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0M();
    } else if (centEstimator == CentralityEstimator::FV0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FV0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFV0A();
    } else {
      LOG(fatal) << "Centrality estimator not recognized for collision selection";
      std::abort();
    }
    hfEvSel.fillHistograms(collision, collRejMask, cent);
    if (collRejMask != 0) {
      return false;
    }
    rowCollisions(mult, collision.numContrib(), cent, collision.posZ());
    return true;
  }

  /// Get charm hadron candidate mass
  /// \param candidate is the charm hadron candidate
  template <DecayChannel channel, typename TCand>
  double getCandMass(const TCand& candidate)
  {
    if constexpr (channel == DecayChannel::DsToKKPi) {
      return hfHelper.invMassDsToKKPi(candidate);
    }
    if constexpr (channel == DecayChannel::DsToPiKK) {
      return hfHelper.invMassDsToPiKK(candidate);
    }
    if constexpr (channel == DecayChannel::DplusToPiKPi) {
      return hfHelper.invMassDplusToPiKPi(candidate);
    }
    if constexpr (channel == DecayChannel::D0ToPiK) {
      return hfHelper.invMassD0ToPiK(candidate);
    }
    if constexpr (channel == DecayChannel::D0ToKPi) {
      return hfHelper.invMassD0barToKPi(candidate);
    }
    return -1.;
  }

  /// Get charm hadron bdt scores
  /// \param candidate is the charm hadron candidate
  template <DecayChannel channel, typename TCand>
  std::vector<float> getCandMlScores(const TCand& candidate)
  {
    std::vector<float> outputMl{-999., -999.};
    if constexpr (channel == DecayChannel::DsToKKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (channel == DecayChannel::DsToPiKK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
      }
    }
    if constexpr (channel == DecayChannel::DplusToPiKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (channel == DecayChannel::D0ToPiK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
      }
    }
    if constexpr (channel == DecayChannel::D0ToKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
      }
    }
    return outputMl;
  }

  /// Fill charm hadron tables
  /// \param candidates are the selected charm hadron candidates
  template <DecayChannel channel, typename TCand>
  void fillCharmHadronTables(TCand const& candidates)
  {
    int indexRedColl = rowCollisions.lastIndex();
    for (const auto& candidate : candidates) {
      if (std::abs(candidate.y(massCharm)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      double massCand = getCandMass<channel>(candidate);
      if constexpr (channel == DecayChannel::D0ToKPi || channel == DecayChannel::D0ToPiK) {
        rowCharmCandidates(indexRedColl, candidate.phi(), candidate.eta(), candidate.pt(), massCand, candidate.prong0Id(), candidate.prong1Id(), -1);
      } else {
        rowCharmCandidates(indexRedColl, candidate.phi(), candidate.eta(), candidate.pt(), massCand, candidate.prong0Id(), candidate.prong1Id(), candidate.prong2Id());
      }
      std::vector<float> outputMl = getCandMlScores<channel>(candidate);
      rowCharmCandidatesMl(indexRedColl, outputMl[0], outputMl[1]);
    }
  }

  /// Fill tracks tables
  /// \param tracks are the selected tracks
  template <typename TTrack>
  void fillTracksTables(TTrack const& tracks)
  {
    int indexRedColl = rowCollisions.lastIndex();
    for (const auto& track : tracks) {
      if (!track.isGlobalTrackWoDCA()) {
        continue;
      }
      rowAssocTrackReduced(indexRedColl, track.globalIndex(), track.phi(), track.eta(), track.pt());
      rowAssocTrackSelInfo(indexRedColl, track.tpcNClsCrossedRows(), track.itsClusterMap(), track.itsNCls(), track.dcaXY(), track.dcaZ());
    }
  }

  // Dplus with ML selections
  void processDplusWithMl(CollsWithCentMult const& colls,
                          CandDplusDataWMl const& candsDplus,
                          TracksData const& tracks)
  {
    for (const auto& coll : colls) {
      auto thisCollId = coll.globalIndex();
      auto candsCThisColl = candsDplus.sliceBy(candsDplusPerCollWMl, thisCollId);
      if (forceCharmInCollision && candsCThisColl.size() < 1) {
        continue;
      }
      if (!checkAndFillCollision(coll)) {
        continue;
      }
      auto trackIdsThisColl = tracks.sliceBy(trackIndicesPerColl, thisCollId);
      fillCharmHadronTables<DecayChannel::DplusToPiKPi>(candsCThisColl);
      fillTracksTables(trackIdsThisColl);
    }
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDplusWithMl, "Process Dplus candidates with ML info", false);

  // Ds with ML selections
  void processDsWithMl(CollsWithCentMult const& colls,
                       TracksData const& tracks,
                       CandDsDataWMl const&)
  {
    for (const auto& coll : colls) {
      auto thisCollId = coll.globalIndex();
      auto candsDsToKKPiWMl = selectedDsToKKPiWMl->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache);
      auto candsDsToPiKKWMl = selectedDsToPiKKWMl->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache);
      if (forceCharmInCollision && candsDsToKKPiWMl.size() < 1 && candsDsToPiKKWMl.size() < 1) {
        continue;
      }
      if (!checkAndFillCollision(coll)) {
        continue;
      }
      auto trackIdsThisColl = tracks.sliceBy(trackIndicesPerColl, thisCollId);
      fillCharmHadronTables<DecayChannel::DsToPiKK>(candsDsToPiKKWMl);
      fillCharmHadronTables<DecayChannel::DsToKKPi>(candsDsToKKPiWMl);
      fillTracksTables(trackIdsThisColl);
    }
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDsWithMl, "Process Ds candidates with ML info", false);

  // D0 with ML selections
  void processD0WithMl(CollsWithCentMult const& colls,
                       TracksData const& tracks,
                       CandD0DataWMl const&)
  {
    for (const auto& coll : colls) {
      auto thisCollId = coll.globalIndex();
      auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache);
      auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache);
      if (forceCharmInCollision && candsD0ToPiKWMl.size() < 1 && candsD0ToKPiWMl.size() < 1) {
        continue;
      }
      if (!checkAndFillCollision(coll)) {
        continue;
      }
      auto trackIdsThisColl = tracks.sliceBy(trackIndicesPerColl, thisCollId);
      fillCharmHadronTables<DecayChannel::D0ToPiK>(candsD0ToPiKWMl);
      fillCharmHadronTables<DecayChannel::D0ToKPi>(candsD0ToKPiWMl);
      fillTracksTables(trackIdsThisColl);
    }
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processD0WithMl, "Process D0 candidates with ML info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorCorrelationsReduced>(cfgc)};
}
