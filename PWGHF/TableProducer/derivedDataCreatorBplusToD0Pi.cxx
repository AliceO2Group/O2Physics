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

/// \file derivedDataCreatorBplusToD0Pi.cxx
/// \brief Producer of derived tables of B+ candidates, collisions and MC particles
/// \note Based on derivedDataCreatorLcToPKPi.cxx
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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"
#include "PWGHF/Utils/utilsPid.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;
using namespace o2::analysis::hf_derived;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorBplusToD0Pi {
  // Candidates
  Produces<o2::aod::HfBplusBases> rowCandidateBase;
  Produces<o2::aod::HfBplusPars> rowCandidatePar;
  Produces<o2::aod::HfBplusParD0s> rowCandidateParD0;
  Produces<o2::aod::HfBplusParEs> rowCandidateParE;
  Produces<o2::aod::HfBplusMls> rowCandidateMl;
  Produces<o2::aod::HfBplusMlD0s> rowCandidateMlD0;
  Produces<o2::aod::HfBplusIds> rowCandidateId;
  Produces<o2::aod::HfBplusMcs> rowCandidateMc;
  // Collisions
  Produces<o2::aod::HfCollBases> rowCollBase;
  Produces<o2::aod::HfCollIds> rowCollId;
  // MC collisions
  Produces<o2::aod::HfMcCollBases> rowMcCollBase;
  Produces<o2::aod::HfMcCollIds> rowMcCollId;
  Produces<o2::aod::HfMcRCollIds> rowMcRCollId;
  // MC particles
  Produces<o2::aod::HfBplusPBases> rowParticleBase;
  Produces<o2::aod::HfBplusPIds> rowParticleId;

  // Switches for filling tables
  Configurable<bool> fillCandidateBase{"fillCandidateBase", true, "Fill candidate base properties"};
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParD0{"fillCandidateParD0", true, "Fill D0 candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateMlD0{"fillCandidateMlD0", true, "Fill D0 candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  Configurable<bool> fillCollBase{"fillCollBase", true, "Fill collision base properties"};
  Configurable<bool> fillCollId{"fillCollId", true, "Fill original collision indices"};
  Configurable<bool> fillMcCollBase{"fillMcCollBase", true, "Fill MC collision base properties"};
  Configurable<bool> fillMcCollId{"fillMcCollId", true, "Fill original MC collision indices"};
  Configurable<bool> fillMcRCollId{"fillMcRCollId", true, "Fill indices of saved derived reconstructed collisions matched to saved derived MC collisions"};
  Configurable<bool> fillParticleBase{"fillParticleBase", true, "Fill MC particle properties"};
  Configurable<bool> fillParticleId{"fillParticleId", true, "Fill original MC indices"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;
  SliceCache cache;
  std::map<int, std::vector<int>> matchedCollisions; // indices of derived reconstructed collisions matched to the global indices of MC collisions
  std::map<int, bool> hasMcParticles;                // flags for MC collisions with HF particles

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfCandBplusMcRec, aod::HfSelBplusToD0Pi>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfMlBplusToD0Pi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfCandBplusMcRec, aod::HfSelBplusToD0Pi, aod::HfMlBplusToD0Pi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>;
  using TypeMcCollisions = aod::McCollisions;
  using THfCandDaughters = aod::HfCand2ProngWPid;
  using THfCandDaughtersMl = soa::Join<THfCandDaughters, aod::HfMlD0>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_bplus::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_bplus::DecayType::BplusToD0Pi));

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_bplus::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_bplus::DecayType::BplusToD0Pi));
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_bplus::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_bplus::DecayType::BplusToD0Pi));
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_bplus::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_bplus::DecayType::BplusToD0Pi));
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_bplus::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_bplus::DecayType::BplusToD0Pi));

  void init(InitContext const&)
  {
    std::array<bool, 8> doprocess{doprocessData, doprocessMcSig, doprocessMcBkg, doprocessMcAll, doprocessDataMl, doprocessMcMlSig, doprocessMcMlBkg, doprocessMcMlAll};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
  }

  template <bool isMC, typename T>
  // void fillTablesCollision(const T& collision, int isEventReject, int runNumber)
  void fillTablesCollision(const T& collision)
  {
    if (fillCollBase) {
      rowCollBase(
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0A(),
        collision.centFT0C(),
        collision.centFT0M(),
        collision.centFV0A(),
        collision.multZeqNTracksPV());
      // isEventReject,
      // runNumber);
    }
    if (fillCollId) {
      rowCollId(
        collision.globalIndex());
    }
    if constexpr (isMC) {
      if (fillMcRCollId && collision.has_mcCollision()) {
        // Save rowCollBase.lastIndex() at key collision.mcCollisionId()
        LOGF(debug, "Rec. collision %d: Filling derived-collision index %d for MC collision %d", collision.globalIndex(), rowCollBase.lastIndex(), collision.mcCollisionId());
        matchedCollisions[collision.mcCollisionId()].push_back(rowCollBase.lastIndex()); // [] inserts an empty element if it does not exist
      }
    }
  }

  template <typename T>
  void fillTablesMcCollision(const T& mcCollision)
  {
    if (fillMcCollBase) {
      rowMcCollBase(
        mcCollision.posX(),
        mcCollision.posY(),
        mcCollision.posZ());
    }
    if (fillMcCollId) {
      rowMcCollId(
        mcCollision.globalIndex());
    }
    if (fillMcRCollId) {
      // Fill the table with the vector of indices of derived reconstructed collisions matched to mcCollision.globalIndex()
      rowMcRCollId(
        matchedCollisions[mcCollision.globalIndex()]);
    }
  }

  template <typename T, typename U, typename V>
  void fillTablesCandidate(const T& candidate, const U& prongCharm, const V& prongBachelor, int candFlag, double invMass,
                           double ct, double y, int8_t flagMc, int8_t origin, float mlScore, const std::vector<float>& mlScoresCharm)
  {
    if (fillCandidateBase) {
      rowCandidateBase(
        rowCollBase.lastIndex(),
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        invMass,
        y);
    }
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
        candidate.ptProng1(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
        prongBachelor.tpcNSigmaPi(),
        prongBachelor.tofNSigmaPi(),
        prongBachelor.tpcTofNSigmaPi(),
        prongBachelor.tpcNSigmaKa(),
        prongBachelor.tofNSigmaKa(),
        prongBachelor.tpcTofNSigmaKa(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct());
    }
    if (fillCandidateParD0) {
      std::array<std::array<std::array<float, 3>, 2>, 2> sigmas{}; // PID nSigma [Expected][Hypothesis][TPC/TOF/TPC+TOF]
      if (candFlag == 0) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], prongCharm, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], prongCharm, 0, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], prongCharm, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], prongCharm, 1, Ka)
      } else if (candFlag == 1) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], prongCharm, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], prongCharm, 1, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], prongCharm, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], prongCharm, 0, Ka)
      }
      rowCandidateParD0(
        prongCharm.cpa(),
        prongCharm.decayLength(),
        prongCharm.impactParameter0(),
        prongCharm.impactParameter1(),
        prongCharm.impactParameterProduct(),
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
        hfHelper.cosThetaStarBplus(candidate),
        ct);
    }
    if (fillCandidateMl) {
      rowCandidateMl(
        mlScore);
    }
    if (fillCandidateMlD0) {
      rowCandidateMlD0(
        mlScoresCharm);
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        prongCharm.prong0Id(),
        prongCharm.prong1Id(),
        candidate.prong1Id());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
    }
  }

  template <typename T, typename U>
  void fillTablesParticle(const T& particle, U mass)
  {
    if (fillParticleBase) {
      rowParticleBase(
        rowMcCollBase.lastIndex(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecayPtEtaPhi::y(particle.pt(), particle.eta(), mass),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
    if (fillParticleId) {
      rowParticleId(
        particle.mcCollisionId(),
        particle.globalIndex());
    }
  }

  template <bool isMl, bool isMc, bool onlyBkg, bool onlySig, typename CollType, typename CandType, typename CandCharmType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
                         CandCharmType const&,
                         TracksWPid const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    if constexpr (isMc) {
      if (fillMcRCollId) {
        matchedCollisions.clear();
      }
    }
    auto sizeTableColl = collisions.size();
    reserveTable(rowCollBase, fillCollBase, sizeTableColl);
    reserveTable(rowCollId, fillCollId, sizeTableColl);
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candidatesThisColl = candidates->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache); // FIXME
      auto sizeTableCand = candidatesThisColl.size();
      LOGF(debug, "Rec. collision %d has %d candidates", thisCollId, sizeTableCand);
      // Skip collisions without HF candidates (and without HF particles in matched MC collisions if saving indices of reconstructed collisions matched to MC collisions)
      bool mcCollisionHasMcParticles{false};
      if constexpr (isMc) {
        mcCollisionHasMcParticles = fillMcRCollId && collision.has_mcCollision() && hasMcParticles[collision.mcCollisionId()];
        LOGF(debug, "Rec. collision %d has MC collision %d with MC particles? %s", thisCollId, collision.mcCollisionId(), mcCollisionHasMcParticles ? "yes" : "no");
      }
      if (sizeTableCand == 0 && (!fillMcRCollId || !mcCollisionHasMcParticles)) {
        LOGF(debug, "Skipping rec. collision %d", thisCollId);
        continue;
      }
      LOGF(debug, "Filling rec. collision %d at derived index %d", thisCollId, rowCollBase.lastIndex() + 1);
      // fillTablesCollision(collision, 0, collision.bc().runNumber());
      fillTablesCollision<isMc>(collision);

      // Fill candidate properties
      reserveTable(rowCandidateBase, fillCandidateBase, sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParD0, fillCandidateParD0, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateMlD0, fillCandidateMlD0, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (isMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (isMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          if constexpr (onlyBkg) {
            if (TESTBIT(std::abs(flagMcRec), aod::hf_cand_bplus::DecayType::BplusToD0Pi)) {
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
            if (!TESTBIT(std::abs(flagMcRec), aod::hf_cand_bplus::DecayType::BplusToD0Pi)) {
              continue;
            }
          }
        }
        auto prongCharm = candidate.template prong0_as<CandCharmType>();
        auto prongBachelor = candidate.template prong1_as<TracksWPid>();
        double ct = hfHelper.ctBplus(candidate);
        double y = hfHelper.yBplus(candidate);
        float massBplusToD0Pi = hfHelper.invMassBplusToD0Pi(candidate);
        float mlScoreBplusToD0Pi{-1.f};
        std::vector<float> mlScoresD0;
        bool isD0 = prongBachelor.sign() < 0; // D0 or D0bar
        if constexpr (isMl) {
          mlScoreBplusToD0Pi = candidate.mlProbBplusToD0Pi();
          if (isD0) {
            std::copy(prongCharm.mlProbD0().begin(), prongCharm.mlProbD0().end(), std::back_inserter(mlScoresD0));
          } else {
            std::copy(prongCharm.mlProbD0bar().begin(), prongCharm.mlProbD0bar().end(), std::back_inserter(mlScoresD0));
          }
        }
        // flag = 0 for D0 pi-, flag = 1 for D0bar pi+
        fillTablesCandidate(candidate, prongCharm, prongBachelor, isD0 ? 0 : 1, massBplusToD0Pi, ct, y, flagMcRec, origin, mlScoreBplusToD0Pi, mlScoresD0);
      }
    }
  }

  template <typename CollisionType, typename ParticleType>
  void preProcessMcCollisions(CollisionType const& mcCollisions,
                              ParticleType const& mcParticles)
  {
    if (!fillMcRCollId) {
      return;
    }
    hasMcParticles.clear();
    // Fill MC collision flags
    for (const auto& mcCollision : mcCollisions) {
      auto thisMcCollId = mcCollision.globalIndex();
      auto particlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, thisMcCollId);
      LOGF(debug, "MC collision %d has %d MC particles (preprocess)", thisMcCollId, particlesThisMcColl.size());
      hasMcParticles[thisMcCollId] = (particlesThisMcColl.size() > 0);
    }
  }

  template <typename CollisionType, typename ParticleType>
  void processMcParticles(CollisionType const& mcCollisions,
                          ParticleType const& mcParticles)
  {
    // Fill MC collision properties
    auto sizeTableMcColl = mcCollisions.size();
    reserveTable(rowMcCollBase, fillMcCollBase, sizeTableMcColl);
    reserveTable(rowMcCollId, fillMcCollId, sizeTableMcColl);
    reserveTable(rowMcRCollId, fillMcRCollId, sizeTableMcColl);
    for (const auto& mcCollision : mcCollisions) {
      auto thisMcCollId = mcCollision.globalIndex();
      auto particlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, thisMcCollId);
      auto sizeTablePart = particlesThisMcColl.size();
      LOGF(debug, "MC collision %d has %d MC particles", thisMcCollId, sizeTablePart);
      // Skip MC collisions without HF particles (and without HF candidates in matched reconstructed collisions if saving indices of reconstructed collisions matched to MC collisions)
      LOGF(debug, "MC collision %d has %d saved derived rec. collisions", thisMcCollId, matchedCollisions[thisMcCollId].size());
      if (sizeTablePart == 0 && (!fillMcRCollId || matchedCollisions[thisMcCollId].empty())) {
        LOGF(debug, "Skipping MC collision %d", thisMcCollId);
        continue;
      }
      LOGF(debug, "Filling MC collision %d at derived index %d", thisMcCollId, rowMcCollBase.lastIndex() + 1);
      fillTablesMcCollision(mcCollision);

      // Fill MC particle properties
      reserveTable(rowParticleBase, fillParticleBase, sizeTablePart);
      reserveTable(rowParticleId, fillParticleId, sizeTablePart);
      for (const auto& particle : particlesThisMcColl) {
        fillTablesParticle(particle, o2::constants::physics::MassBPlus);
      }
    }
  }

  void processData(CollisionsWCentMult const& collisions,
                   SelectedCandidates const&,
                   THfCandDaughters const& candidatesDaughters,
                   TracksWPid const& tracks,
                   aod::BCs const& bcs)
  {
    processCandidates<false, false, false, false>(collisions, candidatesAll, candidatesDaughters, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processData, "Process data", true);

  void processMcSig(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughters const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, false, true>(collisions, candidatesMcSig, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcSig, "Process MC only for signals", false);

  void processMcBkg(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughters const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, true, false>(collisions, candidatesMcBkg, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcBkg, "Process MC only for background", false);

  void processMcAll(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    THfCandDaughters const& candidatesDaughters,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, false, false>(collisions, candidatesMcAll, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcAll, "Process MC", false);

  // ML versions

  void processDataMl(CollisionsWCentMult const& collisions,
                     SelectedCandidatesMl const&,
                     THfCandDaughtersMl const& candidatesDaughters,
                     TracksWPid const& tracks,
                     aod::BCs const& bcs)
  {
    processCandidates<true, false, false, false>(collisions, candidatesMlAll, candidatesDaughters, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processDataMl, "Process data with ML", false);

  void processMcMlSig(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, false, true>(collisions, candidatesMcMlSig, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcMlSig, "Process MC with ML only for signals", false);

  void processMcMlBkg(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, true, false>(collisions, candidatesMcMlBkg, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcMlBkg, "Process MC with ML only for background", false);

  void processMcMlAll(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      THfCandDaughtersMl const& candidatesDaughters,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, false, false>(collisions, candidatesMcMlAll, candidatesDaughters, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorBplusToD0Pi, processMcMlAll, "Process MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorBplusToD0Pi>(cfgc)};
}
