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

/// \file derivedDataCreatorDstarToD0Pi.cxx
/// \brief Producer of derived tables of D*+ candidates, collisions and MC particles
/// \note Based on derivedDataCreatorLcToPKPi.cxx
///
/// \author Mingze Li <mingze.li@cern.cn>, CCNU

#include "PWGHF/Core/DecayChannels.h"
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
using namespace o2::analysis::hf_derived;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorDstarToD0Pi {
  HfProducesDerivedData<
    o2::aod::HfDstarBases,
    o2::aod::HfDstarCollBases,
    o2::aod::HfDstarCollIds,
    o2::aod::HfDstarMcCollBases,
    o2::aod::HfDstarMcCollIds,
    o2::aod::HfDstarMcRCollIds,
    o2::aod::HfDstarPBases,
    o2::aod::HfDstarPIds>
    rowsCommon;
  // Candidates
  Produces<o2::aod::HfDstarPars> rowCandidatePar;
  Produces<o2::aod::HfDstarSels> rowCandidateSel;
  Produces<o2::aod::HfDstarMls> rowCandidateMl;
  Produces<o2::aod::HfDstarIds> rowCandidateId;
  Produces<o2::aod::HfDstarMcs> rowCandidateMc;

  // Switches for filling tables
  HfConfigurableDerivedData confDerData;
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  SliceCache cache;
  static constexpr double mass{o2::constants::physics::MassDStar};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfCandDstarMcRec, aod::HfSelDstarToD0Pi>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfCandDstarMcRec, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandDstarMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Filter filterMcGenMatching = nabs(aod::hf_cand_dstar::flagMcMatchGen) == static_cast<int8_t>(hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_dstar::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_dstar::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_dstar::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_dstar::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);

  void init(InitContext const&)
  {
    std::array<bool, 9> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll, doprocessMcGenOnly};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T, typename U>
  void fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, const U& prongSoftPi, int candFlag, double invMass,
                           double y, int8_t flagMc, int8_t origin, const std::vector<float>& mlScores)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.chi2PCAD0(),
        candidate.cpaD0(),
        candidate.cpaXYD0(),
        candidate.decayLengthD0(),
        candidate.decayLengthXYD0(),
        candidate.decayLengthNormalisedD0(),
        candidate.decayLengthXYNormalisedD0(),
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.pxD0(),
        candidate.pyD0(),
        candidate.pzD0(),
        candidate.pxSoftPi(),
        candidate.pySoftPi(),
        candidate.pzSoftPi(),
        candidate.signSoftPi(),
        candidate.invMassD0(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impParamSoftPi(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
        candidate.normalisedImpParamSoftPi(),
        prong0.tpcNSigmaPi(),
        prong0.tofNSigmaPi(),
        prong0.tpcTofNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaKa(),
        prongSoftPi.tpcNSigmaPi(),
        prongSoftPi.tofNSigmaPi(),
        prongSoftPi.tpcTofNSigmaPi());
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
        candidate.prongPiId());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
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
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (isMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (isMl) {
          if (!TESTBIT(candidate.isSelDstarToD0Pi(), aod::SelectionStep::RecoMl)) {
            continue;
          }
        }
        if constexpr (isMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          if constexpr (onlyBkg) {
            if (std::abs(flagMcRec) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) {
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
            if (std::abs(flagMcRec) != hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) {
              continue;
            }
          }
        }
        auto prong0 = candidate.template prong0_as<TracksWPid>();
        auto prong1 = candidate.template prong1_as<TracksWPid>();
        auto prongSoftPi = candidate.template prongPi_as<TracksWPid>();
        double y = candidate.y(o2::constants::physics::MassDStar);
        double massDstar = candidate.invMassDstar();
        std::vector<float> mlScoresDstarToD0Pi;
        bool isD0 = prongSoftPi.sign() < 0;
        if constexpr (isMl) {
          std::copy(candidate.mlProbDstarToD0Pi().begin(), candidate.mlProbDstarToD0Pi().end(), std::back_inserter(mlScoresDstarToD0Pi));
        }
        fillTablesCandidate(candidate, prong0, prong1, prongSoftPi, isD0 ? 0 : 1, y, massDstar, flagMcRec, origin, mlScoresDstarToD0Pi);
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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processData, "Process data", true);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcSig, "Process MC only for signals", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcBkg, "Process MC only for background", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processDataMl, "Process data with ML", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcMlSig, "Process MC with ML only for signals", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcMlBkg, "Process MC with ML only for background", false);

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
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcMlAll, "Process MC with ML", false);

  void processMcGenOnly(TypeMcCollisions const& mcCollisions,
                        MatchedGenCandidatesMc const& mcParticles)
  {
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorDstarToD0Pi, processMcGenOnly, "Process MC gen. only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorDstarToD0Pi>(cfgc)};
}
