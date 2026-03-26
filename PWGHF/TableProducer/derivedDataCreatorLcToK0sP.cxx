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

/// \file derivedDataCreatorLcToK0sP.cxx
/// \brief Producer of derived tables of Lc candidates, collisions and MC particles
/// \note Based on treeCreatorLcToK0sP.cxx and derivedDataCreatorD0ToKPi.cxx
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, Marietta-Blau Institute
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"
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
using namespace o2::analysis::hf_derived;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorLcToK0sP {
  HfProducesDerivedData<
    o2::aod::HfLcToK0sPBases,
    o2::aod::HfLcToK0sPCollBases,
    o2::aod::HfLcToK0sPCollIds,
    o2::aod::HfLcToK0sPMcCollBases,
    o2::aod::HfLcToK0sPMcCollIds,
    o2::aod::HfLcToK0sPMcRCollIds,
    o2::aod::HfLcToK0sPPBases,
    o2::aod::HfLcToK0sPPIds>
    rowsCommon;
  // Candidates
  Produces<o2::aod::HfLcToK0sPPars> rowCandidatePar;
  Produces<o2::aod::HfLcToK0sPParEs> rowCandidateParE;
  Produces<o2::aod::HfLcToK0sPSels> rowCandidateSel;
  Produces<o2::aod::HfLcToK0sPMls> rowCandidateMl;
  Produces<o2::aod::HfLcToK0sPIds> rowCandidateId;
  Produces<o2::aod::HfLcToK0sPMcs> rowCandidateMc;

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
  static constexpr double Mass{o2::constants::physics::MassLambdaCPlus};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfCandCascadeMcRec, aod::HfSelLcToK0sP>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP, aod::HfMlLcToK0sP>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfCandCascadeMcRec, aod::HfSelLcToK0sP, aod::HfMlLcToK0sP>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandCascadeMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_mc_flag::flagMcMatchGen) == 1;

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == 1;
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) != 1;
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == 1;
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) != 1;

  void init(InitContext const&)
  {
    std::array<bool, 9> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll, doprocessMcGenOnly};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T, typename U>
  void fillTablesCandidate(const T& candidate, const U& bach, int candFlag, double invMass,
                           double ct, double ctV0, double y, int8_t flagMc, int8_t origin, const std::vector<float>& mlScores)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.chi2PCA(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.v0radius(),
        candidate.v0cosPA(),
        candidate.mLambda(),
        candidate.mAntiLambda(),
        candidate.mK0Short(),
        candidate.mGamma(),
        candidate.dcaV0daughters(),
        candidate.ptV0Pos(),
        candidate.ptV0Neg(),
        candidate.decayLengthV0(),
        candidate.dcanegtopv(),
        candidate.dcapostopv(),
        bach.tpcNSigmaPr(),
        bach.tofNSigmaPr(),
        ctV0);
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
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.pxpos(),
        candidate.pypos(),
        candidate.pzpos(),
        candidate.pxneg(),
        candidate.pyneg(),
        candidate.pzneg(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
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
        if constexpr (IsMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          if constexpr (OnlyBkg) {
            if (std::abs(flagMcRec) == 1) {
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
            if (std::abs(flagMcRec) != 1) {
              continue;
            }
          }
        } else {
          if (downSampleBkgFactor < 1.) {
            float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
            if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
              continue;
            }
          }
        }
        auto bach = candidate.template prong0_as<TracksWPid>(); // bachelor
        double const ct = HfHelper::ctLc(candidate);
        double const ctV0 = HfHelper::ctV0K0s(candidate);
        double const y = HfHelper::yLc(candidate);
        float const massLcToK0sP = HfHelper::invMassLcToK0sP(candidate);
        std::vector<float> mlScoresLcToK0sP;
        if constexpr (IsMl) {
          std::copy(candidate.mlProbLcToK0sP().begin(), candidate.mlProbLcToK0sP().end(), std::back_inserter(mlScoresLcToK0sP));
        }
        if (candidate.isSelLcToK0sP()) {
          fillTablesCandidate(candidate, bach, 0, massLcToK0sP, ct, ctV0, y, flagMcRec, origin, mlScoresLcToK0sP);
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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processData, "Process data", true);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcSig, "Process MC only for signals", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcBkg, "Process MC only for background", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processDataMl, "Process data with ML", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcMlSig, "Process MC with ML only for signals", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcMlBkg, "Process MC with ML only for background", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcMlAll, "Process MC with ML", false);

  void processMcGenOnly(TypeMcCollisions const& mcCollisions,
                        MatchedGenCandidatesMc const& mcParticles)
  {
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToK0sP, processMcGenOnly, "Process MC gen. only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorLcToK0sP>(cfgc)};
}
