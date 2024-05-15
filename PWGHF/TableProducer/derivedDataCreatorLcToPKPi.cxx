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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorLcToPKPi {
  // Candidates
  Produces<o2::aod::Hf3PBases> rowCandidateBase;
  Produces<o2::aod::Hf3PPars> rowCandidatePar;
  Produces<o2::aod::Hf3PParEs> rowCandidateParE;
  Produces<o2::aod::Hf3PSels> rowCandidateSel;
  Produces<o2::aod::Hf3PMls> rowCandidateMl;
  Produces<o2::aod::Hf3PIds> rowCandidateId;
  Produces<o2::aod::Hf3PMcs> rowCandidateMc;
  // Collisions
  Produces<o2::aod::Hf3PCollBases> rowCollBase;
  Produces<o2::aod::Hf3PCollIds> rowCollId;
  // MC collisions
  Produces<o2::aod::Hf3PMcCollBases> rowMcCollBase;
  Produces<o2::aod::Hf3PMcCollIds> rowMcCollId;
  // MC particles
  Produces<o2::aod::Hf3PPBases> rowParticleBase;
  Produces<o2::aod::Hf3PPIds> rowParticleId;

  // Switches for filling tables
  Configurable<bool> fillCandidateBase{"fillCandidateBase", true, "Fill candidate base properties"};
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  Configurable<bool> fillCollBase{"fillCollBase", true, "Fill collision base properties"};
  Configurable<bool> fillCollId{"fillCollId", true, "Fill original collision indices"};
  Configurable<bool> fillMcCollBase{"fillMcCollBase", true, "Fill MC collision base properties"};
  Configurable<bool> fillMcCollId{"fillMcCollId", true, "Fill indices of saved derived reconstructed collisions matched to saved derived MC collisions"};
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
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using TypeMcCollisions = aod::McCollisions;

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
  }

  template <typename T>
  void reserveTable(T& table, const Configurable<bool>& enabled, const uint64_t size)
  {
    if (enabled.value) {
      table.reserve(size);
    }
  };

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
      if (fillMcCollId && collision.has_mcCollision()) {
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
      // Fill the table with the vector of indices of derived reconstructed collisions matched to mcCollision.globalIndex()
      rowMcCollId(
        matchedCollisions[mcCollision.globalIndex()]);
    }
  }

  template <typename T, typename U>
  void fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, const U& prong2, int candFlag, double invMass,
                           double ct, int8_t flagMc, int8_t origin, int8_t swapping, const std::vector<float>& mlScores)
  {
    if (fillCandidateBase) {
      rowCandidateBase(
        rowCollBase.lastIndex(),
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        invMass);
    }
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
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
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

  template <typename T, typename U>
  void fillTablesParticle(const T& particle, U /*mass*/)
  {
    if (fillParticleBase) {
      rowParticleBase(
        rowMcCollBase.lastIndex(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
    if (fillParticleId) {
      rowParticleId(
        particle.mcCollisionId(),
        particle.globalIndex());
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
      if (fillMcCollId) {
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
        mcCollisionHasMcParticles = fillMcCollId && collision.has_mcCollision() && hasMcParticles[collision.mcCollisionId()];
        LOGF(debug, "Rec. collision %d has MC collision %d with MC particles? %s", thisCollId, collision.mcCollisionId(), mcCollisionHasMcParticles ? "yes" : "no");
      }
      if (sizeTableCand == 0 && (!fillMcCollId || !mcCollisionHasMcParticles)) {
        LOGF(debug, "Skipping rec. collision %d", thisCollId);
        continue;
      }
      LOGF(debug, "Filling rec. collision %d at derived index %d", thisCollId, rowCollBase.lastIndex() + 1);
      // fillTablesCollision(collision, 0, collision.bc().runNumber());
      fillTablesCollision<isMc>(collision);

      // Fill candidate properties
      reserveTable(rowCandidateBase, fillCandidateBase, sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
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
              float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
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
        float massLcToPKPi = hfHelper.invMassLcToPKPi(candidate);
        float massLcToPiKP = hfHelper.invMassLcToPiKP(candidate);
        std::vector<float> mlScoresLcToPKPi, mlScoresLcToPiKP;
        if constexpr (isMl) {
          std::copy(candidate.mlProbLcToPKPi().begin(), candidate.mlProbLcToPKPi().end(), std::back_inserter(mlScoresLcToPKPi));
          std::copy(candidate.mlProbLcToPiKP().begin(), candidate.mlProbLcToPiKP().end(), std::back_inserter(mlScoresLcToPiKP));
        }
        if (candidate.isSelLcToPKPi()) {
          fillTablesCandidate(candidate, prong0, prong1, prong2, 0, massLcToPKPi, ct, flagMcRec, origin, swapping, mlScoresLcToPKPi);
        }
        if (candidate.isSelLcToPiKP()) {
          fillTablesCandidate(candidate, prong0, prong1, prong2, 1, massLcToPiKP, ct, flagMcRec, origin, swapping, mlScoresLcToPiKP);
        }
      }
    }
  }

  template <typename CollisionType, typename ParticleType>
  void preProcessMcCollisions(CollisionType const& mcCollisions,
                              ParticleType const& mcParticles)
  {
    if (!fillMcCollId) {
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
    for (const auto& mcCollision : mcCollisions) {
      auto thisMcCollId = mcCollision.globalIndex();
      auto particlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, thisMcCollId);
      auto sizeTablePart = particlesThisMcColl.size();
      LOGF(debug, "MC collision %d has %d MC particles", thisMcCollId, sizeTablePart);
      // Skip MC collisions without HF particles (and without HF candidates in matched reconstructed collisions if saving indices of reconstructed collisions matched to MC collisions)
      LOGF(debug, "MC collision %d has %d saved derived rec. collisions", thisMcCollId, matchedCollisions[thisMcCollId].size());
      if (sizeTablePart == 0 && (!fillMcCollId || matchedCollisions[thisMcCollId].empty())) {
        LOGF(debug, "Skipping MC collision %d", thisMcCollId);
        continue;
      }
      LOGF(debug, "Filling MC collision %d at derived index %d", thisMcCollId, rowMcCollBase.lastIndex() + 1);
      fillTablesMcCollision(mcCollision);

      // Fill MC particle properties
      reserveTable(rowParticleBase, fillParticleBase, sizeTablePart);
      reserveTable(rowParticleId, fillParticleId, sizeTablePart);
      for (const auto& particle : particlesThisMcColl) {
        fillTablesParticle(particle, o2::constants::physics::MassLambdaCPlus);
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
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcSig, "Process MC only for signals", false);

  void processMcBkg(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcBkg, "Process MC only for background", false);

  void processMcAll(CollisionsWMcCentMult const& collisions,
                    SelectedCandidatesMc const&,
                    TypeMcCollisions const& mcCollisions,
                    MatchedGenCandidatesMc const& mcParticles,
                    TracksWPid const& tracks,
                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<false, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
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
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, false, true>(collisions, candidatesMcMlSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlSig, "Process MC with ML only for signals", false);

  void processMcMlBkg(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, true, false>(collisions, candidatesMcMlBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlBkg, "Process MC with ML only for background", false);

  void processMcMlAll(CollisionsWMcCentMult const& collisions,
                      SelectedCandidatesMcMl const&,
                      TypeMcCollisions const& mcCollisions,
                      MatchedGenCandidatesMc const& mcParticles,
                      TracksWPid const& tracks,
                      aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<true, true, false, false>(collisions, candidatesMcMlAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcMlAll, "Process MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorLcToPKPi>(cfgc)};
}
