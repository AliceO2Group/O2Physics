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

/// \file derivedDataCreatorXicToXiPiPi.cxx
/// \brief Producer of derived tables of Ξc± → (Ξ∓ → (Λ → p π∓) π∓) π± π± candidates, collisions and MC particles
/// \note Based on derivedDataCreatorBplusToD0Pi.cxx
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

// #include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;
using namespace o2::analysis::hf_derived;
using namespace o2::aod::hf_cand_xic_to_xi_pi_pi;
// using namespace o2::hf_decay::hf_cand_beauty;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorXicToXiPiPi {
  HfProducesDerivedData<
    o2::aod::HfXicToXiPiPiBases,
    o2::aod::HfXicToXiPiPiCollBases,
    o2::aod::HfXicToXiPiPiCollIds,
    o2::aod::HfXicToXiPiPiMcCollBases,
    o2::aod::HfXicToXiPiPiMcCollIds,
    o2::aod::HfXicToXiPiPiMcRCollIds,
    o2::aod::HfXicToXiPiPiPBases,
    o2::aod::HfXicToXiPiPiPIds>
    rowsCommon;
  // Candidates
  Produces<o2::aod::HfXicToXiPiPiPars> rowCandidatePar;
  Produces<o2::aod::HfXicToXiPiPiParXis> rowCandidateParXi;
  Produces<o2::aod::HfXicToXiPiPiParEs> rowCandidateParE;
  Produces<o2::aod::HfXicToXiPiPiSels> rowCandidateSel;
  Produces<o2::aod::HfXicToXiPiPiMls> rowCandidateMl;
  Produces<o2::aod::HfXicToXiPiPiMlXis> rowCandidateMlXi;
  Produces<o2::aod::HfXicToXiPiPiIds> rowCandidateId;
  Produces<o2::aod::HfXicToXiPiPiMcs> rowCandidateMc;

  // Switches for filling tables
  HfConfigurableDerivedData confDerData;
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParXi{"fillCandidateParXi", true, "Fill Xi candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateMlXi{"fillCandidateMlXi", true, "Fill Xi candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;
  SliceCache cache;
  static constexpr double mass{o2::constants::physics::MassXiCPlus};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandXicMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using THfCandDaughtersMl = soa::Join<aod::Cascades>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToXiPiPi & static_cast<int8_t>(BIT(aod::SelectionStep::RecoMl - 1))) != 0;
  Filter filterMcGenMatching = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchGen) == static_cast<int8_t>(DecayType::XicToXiPiPi);

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == static_cast<int8_t>(DecayType::XicToXiPiPi);
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != static_cast<int8_t>(DecayType::XicToXiPiPi);
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == static_cast<int8_t>(DecayType::XicToXiPiPi);
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != static_cast<int8_t>(DecayType::XicToXiPiPi);

  void init(InitContext const&)
  {
    std::array<bool, 9> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll, doprocessMcGenOnly};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T, typename U, typename V>
  void fillTablesCandidate(const T& candidate, const U& prongCascade, const V& prongBachelor0, int candFlag, double invMass,
                           double ct, double y, int8_t flagMc, int8_t origin, float mlScore, const std::vector<float>& mlScoresCascade)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.chi2PCA(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1());
    }
    if (fillCandidateParXi) {
      std::array<std::array<std::array<float, 3>, 2>, 2> sigmas{}; // PID nSigma [Expected][Hypothesis][TPC/TOF/TPC+TOF]
      if (candFlag == 0) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], prongCascade, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], prongCascade, 0, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], prongCascade, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], prongCascade, 1, Ka)
      } else if (candFlag == 1) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], prongCascade, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], prongCascade, 1, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], prongCascade, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], prongCascade, 0, Ka)
      }
      rowCandidateParXi(
        prongCascade.cpa(),
        prongCascade.decayLength(),
        prongCascade.impactParameter0(),
        prongCascade.impactParameter1(),
        prongCascade.impactParameterProduct(),
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][0],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][1],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][2],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][0],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][1],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][2],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][0],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][1],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][2],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][0],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][1],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][2]);
    }
    if (fillCandidateParE) {
      rowCandidateParE(
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.errorDecayLength(),
        candidate.errorDecayLengthXY(),
        candidate.rSecondaryVertex(),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.errorImpactParameter1(),
        ct);
    }
    if (fillCandidateSel) {
      rowCandidateSel(
        BIT(candFlag));
    }
    if (fillCandidateMl) {
      rowCandidateMl(
        mlScore);
    }
    if (fillCandidateMlXi) {
      rowCandidateMlXi(
        mlScoresCascade);
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        prongCascade.prong0Id(),
        prongCascade.prong1Id(),
        candidate.prong1Id());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
    }
  }

  template <bool isMl, bool isMc, bool onlyBkg, bool onlySig, typename CollType, typename CandType, typename CandCascadeType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
                         CandCascadeType const&,
                         TracksWPid const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    if constexpr (isMc) {
      if (confDerData.fillMcRCollId) {
        rowsCommon.matchedCollisions.clear();
      }
    }
    auto sizeTableColl = collisions.size();
    rowsCommon.reserveTablesColl(sizeTableColl);
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candidatesThisColl = candidates->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache); // FIXME
      auto sizeTableCand = candidatesThisColl.size();
      LOGF(debug, "Rec. collision %d has %d candidates", thisCollId, sizeTableCand);
      // Skip collisions without HF candidates (and without HF particles in matched MC collisions if saving indices of reconstructed collisions matched to MC collisions)
      bool mcCollisionHasMcParticles{false};
      if constexpr (isMc) {
        mcCollisionHasMcParticles = confDerData.fillMcRCollId && collision.has_mcCollision() && rowsCommon.hasMcParticles[collision.mcCollisionId()];
        LOGF(debug, "Rec. collision %d has MC collision %d with MC particles? %s", thisCollId, collision.mcCollisionId(), mcCollisionHasMcParticles ? "yes" : "no");
      }
      if (sizeTableCand == 0 && (!confDerData.fillMcRCollId || !mcCollisionHasMcParticles)) {
        LOGF(debug, "Skipping rec. collision %d", thisCollId);
        continue;
      }
      LOGF(debug, "Filling rec. collision %d at derived index %d", thisCollId, rowsCommon.rowCollBase.lastIndex() + 1);
      rowsCommon.fillTablesCollision<isMc>(collision);

      // Fill candidate properties
      rowsCommon.reserveTablesCandidates(sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParXi, fillCandidateParXi, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateMlXi, fillCandidateMlXi, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (isMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (isMl) {
          if (!TESTBIT(candidate.isSelXicToXiPiPi(), aod::SelectionStep::RecoMl)) {
            continue;
          }
        }
        if constexpr (isMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originRec();
          if constexpr (onlyBkg) {
            if (std::abs(flagMcRec) == DecayType::XicToXiPiPi) {
              continue;
            }
            if (downSampleBkgFactor < 1.) {
              float pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
              if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
                continue;
              }
            }
          }
          if constexpr (onlySig) {
            if (std::abs(flagMcRec) != DecayType::XicToXiPiPi) {
              continue;
            }
          }
        }
        // FIXME
        auto prongCascade = candidate.template cascade_as<CandCascadeType>();
        auto prongBachelor0 = candidate.template pi0_as<TracksWPid>();
        auto prongBachelor1 = candidate.template pi1_as<TracksWPid>();
        auto prongBachelorWhat = candidate.template bachelor_as<TracksWPid>();
        auto prongPosTrack = candidate.template posTrack_as<TracksWPid>();
        auto prongNegTrack = candidate.template negTrack_as<TracksWPid>();

        double ct = hfHelper.ctXic(candidate);
        double y = hfHelper.yXic(candidate);
        float massXicToXiPiPi = candidate.invMassXicPlus();
        float mlScoreXicToXiPiPi{-1.f};
        std::vector<float> mlScoresXi;
        bool isXi = prongBachelor0.sign() < 0;
        if (isXi) { // is Xi
          std::copy(prongCascade.mlProbXi().begin(), prongCascade.mlProbXi().end(), std::back_inserter(mlScoresXi));
        } else { // is Xibar
          std::copy(prongCascade.mlProbXibar().begin(), prongCascade.mlProbXibar().end(), std::back_inserter(mlScoresXi));
        }
        if constexpr (isMl) {
          mlScoreXicToXiPiPi = candidate.mlProbXicToXiPiPi();
        }
        // flag = 0 for Xi pi-, flag = 1 for Xibar pi+
        fillTablesCandidate(candidate, prongCascade, prongBachelor0, isXi ? 0 : 1, massXicToXiPiPi, ct, y, flagMcRec, origin, mlScoreXicToXiPiPi, mlScoresXi);
      }
    }
  }

  void processData(CollisionsWCentMult const& collisions,
                   SelectedCandidates const&,
                   THfCandDaughtersMl const& candidatesDaughters,
                   TracksWPid const& tracks,
                   aod::BCs const& bcs)
  {
    processCandidates<false, false, false, false>(collisions, candidatesAll, candidatesDaughters, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processData, "Process data", true);

  void processMcSig(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughtersMl const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, true>(collisions, candidatesMcSig, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcSig, "Process MC only for signals", false);

  void processMcBkg(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughtersMl const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, true, false>(collisions, candidatesMcBkg, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcBkg, "Process MC only for background", false);

  void processMcAll(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughtersMl const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, false>(collisions, candidatesMcAll, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     THfCandDaughtersMl const& candidatesDaughters,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, candidatesDaughters, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processDataMl, "Process data with ML", false);

  void processMcMlSig(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, true>(collisions, candidatesMcMlSig, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlSig, "Process MC with ML only for signals", false);

  void processMcMlBkg(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, true, false>(collisions, candidatesMcMlBkg, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlBkg, "Process MC with ML only for background", false);

  void processMcMlAll(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, false>(collisions, candidatesMcMlAll, candidatesDaughters, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlAll, "Process MC with ML", false);

  void processMcGenOnly(TypeMcCollisions const& mcCollisions,
                        MatchedGenCandidatesMc const& mcParticles)
  {
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcGenOnly, "Process MC gen. only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorXicToXiPiPi>(cfgc)};
}
