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

#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"
#include "PWGLF/DataModel/mcCentrality.h"

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
  Produces<o2::aod::HfXicToXiPiPiParEs> rowCandidateParE;
  Produces<o2::aod::HfXicToXiPiPiSels> rowCandidateSel;
  Produces<o2::aod::HfXicToXiPiPiMls> rowCandidateMl;
  Produces<o2::aod::HfXicToXiPiPiIds> rowCandidateId;
  Produces<o2::aod::HfXicToXiPiPiMcs> rowCandidateMc;

  // Switches for filling tables
  HfConfigurableDerivedData confDerData;
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  SliceCache cache;
  static constexpr double Mass{o2::constants::physics::MassXiCPlus};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandXicMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using THfCandDaughtersMl = aod::Cascades;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToXiPiPi & static_cast<int>(BIT(o2::aod::hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoMl - 1))) != 0;
  Filter filterMcGenMatching = aod::hf_cand_mc_flag::flagMcMatchGen != 0;

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
  Partition<SelectedCandidatesMc> candidatesMcSig = aod::hf_cand_mc_flag::flagMcMatchRec != 0;
  Partition<SelectedCandidatesMc> candidatesMcBkg = aod::hf_cand_mc_flag::flagMcMatchRec == 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = aod::hf_cand_mc_flag::flagMcMatchRec != 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = aod::hf_cand_mc_flag::flagMcMatchRec == 0;

  void init(InitContext const&)
  {
    std::array<bool, 9> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll, doprocessMcGenOnly};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T>
  void fillTablesCandidate(const T& candidate, int candFlag, double invMass,
                           double ct, double y, int8_t flagMc, int8_t origin, const std::vector<float>& mlScores)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.sign(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.invMassXi(),
        candidate.invMassLambda(),
        candidate.invMassXiPi0(),
        candidate.invMassXiPi1(),
        candidate.chi2PCA(),
        ct,
        candidate.decayLength(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXY(),
        candidate.decayLengthXYNormalised(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.cpaXi(),
        candidate.cpaXYXi(),
        candidate.cpaLambda(),
        candidate.cpaXYLambda(),
        candidate.impactParameter0(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameter1(),
        candidate.impactParameterNormalised1(),
        candidate.impactParameter2(),
        candidate.impactParameterNormalised2(),
        candidate.maxNormalisedDeltaIP());
    }
    if (fillCandidateParE) {
      rowCandidateParE(
        candidate.cpaLambdaToXi(),
        candidate.cpaXYLambdaToXi(),
        candidate.pProng1(),
        candidate.pProng2(),
        candidate.pBachelorPi(),
        candidate.pPiFromLambda(),
        candidate.pPrFromLambda(),
        candidate.dcaXiDaughters(),
        candidate.dcaV0Daughters(),
        candidate.dcaPosToPV(),
        candidate.dcaNegToPV(),
        candidate.dcaBachelorToPV(),
        candidate.dcaXYCascToPV(),
        candidate.dcaZCascToPV(),
        candidate.nSigTpcPiFromXicPlus0(),
        candidate.nSigTpcPiFromXicPlus1(),
        candidate.nSigTpcBachelorPi(),
        candidate.nSigTpcPiFromLambda(),
        candidate.nSigTpcPrFromLambda(),
        candidate.nSigTofPiFromXicPlus0(),
        candidate.nSigTofPiFromXicPlus1(),
        candidate.nSigTofBachelorPi(),
        candidate.nSigTofPiFromLambda(),
        candidate.nSigTofPrFromLambda());
    }
    if (fillCandidateSel) {
      rowCandidateSel(
        BIT(candFlag));
    }
    if (fillCandidateMl) {
      rowCandidateMl(
        mlScores);
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        candidate.pi0Id(),
        candidate.pi1Id(),
        candidate.bachelorId(),
        candidate.posTrackId(),
        candidate.negTrackId());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
    }
  }

  template <bool IsMl, bool IsMc, bool OnlyBkg, bool OnlySig, typename CollType, typename CandType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
                         TracksWPid const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    if constexpr (IsMc) {
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
      if constexpr (IsMc) {
        mcCollisionHasMcParticles = confDerData.fillMcRCollId && collision.has_mcCollision() && rowsCommon.hasMcParticles[collision.mcCollisionId()];
        LOGF(debug, "Rec. collision %d has MC collision %d with MC particles? %s", thisCollId, collision.mcCollisionId(), mcCollisionHasMcParticles ? "yes" : "no");
      }
      if (sizeTableCand == 0 && (!confDerData.fillMcRCollId || !mcCollisionHasMcParticles)) {
        LOGF(debug, "Skipping rec. collision %d", thisCollId);
        continue;
      }
      LOGF(debug, "Filling rec. collision %d at derived index %d", thisCollId, rowsCommon.rowCollBase.lastIndex() + 1);
      rowsCommon.fillTablesCollision<IsMc>(collision);

      // Fill candidate properties
      rowsCommon.reserveTablesCandidates(sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (IsMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (IsMl) {
          if (!TESTBIT(candidate.isSelXicToXiPiPi(), o2::aod::hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoMl)) {
            continue;
          }
        }
        if constexpr (IsMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          if constexpr (OnlyBkg) {
            if (TESTBIT(std::abs(flagMcRec), DecayType::XicToXiPiPi)) {
              continue;
            }
            if (downSampleBkgFactor < 1.) {
              float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
              if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
                continue;
              }
            }
          }
          if constexpr (OnlySig) {
            if (!TESTBIT(std::abs(flagMcRec), DecayType::XicToXiPiPi)) {
              continue;
            }
          }
        }
        float const massXicToXiPiPi = candidate.invMassXicPlus();
        double const ct = HfHelper::ctXic(candidate);
        double const y = HfHelper::yXic(candidate);
        std::vector<float> mlScoresXicToXiPiPi;
        if constexpr (IsMl) {
          std::copy(candidate.mlProbXicToXiPiPi().begin(), candidate.mlProbXicToXiPiPi().end(), std::back_inserter(mlScoresXicToXiPiPi));
        }
        // FIXME: Remove candFlag?
        fillTablesCandidate(candidate, 1, massXicToXiPiPi, ct, y, flagMcRec, origin, mlScoresXicToXiPiPi);
      }
    }
  }

  void processData(CollisionsWCentMult const& collisions,
                   SelectedCandidates const&,
                   TracksWPid const& tracks,
                   aod::BCs const& bcs)
  {
    processCandidates<false, false, false, false>(collisions, candidatesAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processData, "Process data", true);

  void processMcSig(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcSig, "Process MC only for signals", false);

  void processMcBkg(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcBkg, "Process MC only for background", false);

  void processMcAll(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processDataMl, "Process data with ML", false);

  void processMcMlSig(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, true>(collisions, candidatesMcMlSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlSig, "Process MC with ML only for signals", false);

  void processMcMlBkg(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, true, false>(collisions, candidatesMcMlBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlBkg, "Process MC with ML only for background", false);

  void processMcMlAll(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, false>(collisions, candidatesMcMlAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcMlAll, "Process MC with ML", false);

  void processMcGenOnly(TypeMcCollisions const& mcCollisions,
                        MatchedGenCandidatesMc const& mcParticles)
  {
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorXicToXiPiPi, processMcGenOnly, "Process MC gen. only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorXicToXiPiPi>(cfgc)};
}
