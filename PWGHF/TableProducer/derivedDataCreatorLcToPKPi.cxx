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

/// \file derivedDataCreatorLcToPKPi.cxx
/// \brief Producer of derived tables of Lc candidates, collisions and MC particles
/// \note Based on treeCreatorLcToPKPi.cxx and derivedDataCreatorD0ToKPi.cxx
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include <algorithm>
#include <map>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGLF/DataModel/mcCentrality.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_derived;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorLcToPKPi {
  HfProducesDerivedData<
    o2::aod::HfLcBases,
    o2::aod::HfLcCollBases,
    o2::aod::HfLcCollIds,
    o2::aod::HfLcMcCollBases,
    o2::aod::HfLcMcCollIds,
    o2::aod::HfLcMcRCollIds,
    o2::aod::HfLcPBases,
    o2::aod::HfLcPIds>
    rowsCommon;
  // Candidates
  Produces<o2::aod::HfLcPars> rowCandidatePar;
  Produces<o2::aod::HfLcParEs> rowCandidateParE;
  Produces<o2::aod::HfLcSels> rowCandidateSel;
  Produces<o2::aod::HfLcMls> rowCandidateMl;
  Produces<o2::aod::HfLcIds> rowCandidateId;
  Produces<o2::aod::HfLcMcs> rowCandidateMc;

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

  HfHelper hfHelper;
  SliceCache cache;
  static constexpr double mass{o2::constants::physics::MassLambdaCPlus};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 1 || aod::hf_sel_candidate_lc::isSelLcToPiKP >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));

  void init(InitContext const&)
  {
    std::array<bool, 8> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T, typename U>
  void fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, const U& prong2, int candFlag, double invMass,
                           double ct, double y, int8_t flagMc, int8_t origin, int8_t swapping, const std::vector<float>& mlScores)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.chi2PCA(),
        candidate.nProngsContributorsPV(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
        candidate.impactParameterNormalised2(),
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaPr(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaPr(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaPr(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaKa(),
        prong2.tpcNSigmaPi(),
        prong2.tpcNSigmaPr(),
        prong2.tofNSigmaPi(),
        prong2.tofNSigmaPr(),
        prong2.tpcTofNSigmaPi(),
        prong2.tpcTofNSigmaPr());
    }
    if (fillCandidateParE) {
      rowCandidateParE(
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.errorDecayLength(),
        candidate.errorDecayLengthXY(),
        candidate.rSecondaryVertex(),
        RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
        RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.pxProng2(),
        candidate.pyProng2(),
        candidate.pzProng2(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        candidate.errorImpactParameter2(),
        ct);
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
        candidate.prong0Id(),
        candidate.prong1Id(),
        candidate.prong2Id());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin,
        swapping);
    }
  }

  template <bool isMl, bool isMc, bool onlyBkg, bool onlySig, typename CollType, typename CandType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
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
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (isMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0, swapping = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (isMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          swapping = candidate.isCandidateSwapped();
          if constexpr (onlyBkg) {
            if (TESTBIT(std::abs(flagMcRec), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
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
            if (!TESTBIT(std::abs(flagMcRec), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
              continue;
            }
          }
        }
        auto prong0 = candidate.template prong0_as<TracksWPid>();
        auto prong1 = candidate.template prong1_as<TracksWPid>();
        auto prong2 = candidate.template prong2_as<TracksWPid>();
        double ct = hfHelper.ctLc(candidate);
        double y = hfHelper.yLc(candidate);
        float massLcToPKPi = hfHelper.invMassLcToPKPi(candidate);
        float massLcToPiKP = hfHelper.invMassLcToPiKP(candidate);
        std::vector<float> mlScoresLcToPKPi, mlScoresLcToPiKP;
        if constexpr (isMl) {
          std::copy(candidate.mlProbLcToPKPi().begin(), candidate.mlProbLcToPKPi().end(), std::back_inserter(mlScoresLcToPKPi));
          std::copy(candidate.mlProbLcToPiKP().begin(), candidate.mlProbLcToPiKP().end(), std::back_inserter(mlScoresLcToPiKP));
        }
        if (candidate.isSelLcToPKPi()) {
          fillTablesCandidate(candidate, prong0, prong1, prong2, 0, massLcToPKPi, ct, y, flagMcRec, origin, swapping, mlScoresLcToPKPi);
        }
        if (candidate.isSelLcToPiKP()) {
          fillTablesCandidate(candidate, prong0, prong1, prong2, 1, massLcToPiKP, ct, y, flagMcRec, origin, swapping, mlScoresLcToPiKP);
        }
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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processData, "Process data", true);

  void processMcSig(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcSig, "Process MC only for signals", false);

  void processMcBkg(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcBkg, "Process MC only for background", false);

  void processMcAll(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<false, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processDataMl, "Process data with ML", false);

  void processMcMlSig(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, true>(collisions, candidatesMcMlSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlSig, "Process MC with ML only for signals", false);

  void processMcMlBkg(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, true, false>(collisions, candidatesMcMlBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlBkg, "Process MC with ML only for background", false);

  void processMcMlAll(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<true, true, false, false>(collisions, candidatesMcMlAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlAll, "Process MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorLcToPKPi>(cfgc)};
}
