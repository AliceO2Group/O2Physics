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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorLcToPKPi {
  Produces<o2::aod::Hf3PBases> rowCandidateBase;
  Produces<o2::aod::Hf3PPars> rowCandidatePar;
  Produces<o2::aod::Hf3PParEs> rowCandidateParE;
  Produces<o2::aod::Hf3PSels> rowCandidateSel;
  Produces<o2::aod::Hf3PIds> rowCandidateId;
  Produces<o2::aod::Hf3PMcs> rowCandidateMc;
  Produces<o2::aod::Hf3PCollBases> rowCollBase;
  Produces<o2::aod::Hf3PCollIds> rowCollId;
  Produces<o2::aod::Hf3PPBases> rowParticleBase;
  Produces<o2::aod::Hf3PPIds> rowParticleId;

  // Switches for filling tables
  Configurable<bool> fillCandidateBase{"fillCandidateBase", true, "Fill candidate base properties"};
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill candidate indices"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  Configurable<bool> fillCollBase{"fillCollBase", true, "Fill collision base properties"};
  Configurable<bool> fillCollId{"fillCollId", true, "Fill collision indices"};
  Configurable<bool> fillParticleBase{"fillParticleBase", true, "Fill particle properties"};
  Configurable<bool> fillParticleId{"fillParticleId", true, "Fill particle indices"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;
  SliceCache cache;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 1 || aod::hf_sel_candidate_lc::isSelLcToPiKP >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_lc::isSelLcToPKPi >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi));

  void init(InitContext const&)
  {
    std::array<bool, 4> doprocess{doprocessDataWithDCAFitterN, doprocessMcWithDCAFitterSig, doprocessMcWithDCAFitterBkg, doprocessMcWithDCAFitterAll};
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

  template <typename T>
  void fillTablesCollision(const T& collision, int isEventReject, int runNumber)
  {
    if (fillCollBase) {
      rowCollBase(
        collision.numContrib(),
        isEventReject,
        runNumber);
    }
    if (fillCollId) {
      rowCollId(
        collision.globalIndex());
    }
  }

  template <typename T, typename U>
  auto fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, int candFlag, double invMass,
                           double ct, int8_t flagMc, int8_t origin)
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
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaKa(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaKa(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        candidate.cpa(),
        candidate.cpaXY());
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
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        ct);
    }
    if (fillCandidateSel) {
      rowCandidateSel(
        BIT(candFlag));
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        candidate.prong0Id(),
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

  template <bool isMc, bool onlyBkg, bool onlySig, typename CandType>
  void processCandidates(aod::Collisions const& collisions,
                         Partition<CandType>& candidates,
                         TracksWPid const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    auto sizeTableColl = collisions.size();
    reserveTable(rowCollBase, fillCollBase, sizeTableColl);
    reserveTable(rowCollId, fillCollId, sizeTableColl);
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candidatesThisColl = candidates->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache); // FIXME
      auto sizeTableCand = candidatesThisColl.size();
      // skip collisions without HF candidates
      if (sizeTableCand == 0) {
        continue;
      }
      fillTablesCollision(collision, 0, collision.bc().runNumber());

      // Fill candidate properties
      reserveTable(rowCandidateBase, fillCandidateBase, sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (isMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec, origin;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (isMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
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
        } else {
          flagMcRec = 0;
          origin = 0;
        }
        auto prong0 = candidate.template prong0_as<TracksWPid>();
        auto prong1 = candidate.template prong1_as<TracksWPid>();
        double ct = hfHelper.ctLc(candidate);
        float massLcToPKPi = hfHelper.invMassLcToPKPi(candidate);
        float massLcToPiKP = hfHelper.invMassLcToPiKP(candidate);
        if (candidate.isSelLcToPKPi()) {
          fillTablesCandidate(candidate, prong0, prong1, 0, massLcToPKPi, ct, flagMcRec, origin);
        }
        if (candidate.isSelLcToPiKP()) {
          fillTablesCandidate(candidate, prong0, prong1, 1, massLcToPiKP, ct, flagMcRec, origin);
        }
      }
    }
  }

  template <typename ParticleType>
  void processMcParticles(ParticleType const& mcParticles)
  {
    // Fill MC particle properties
    auto sizeTablePart = mcParticles.size();
    reserveTable(rowParticleBase, fillParticleBase, sizeTablePart);
    reserveTable(rowParticleId, fillParticleId, sizeTablePart);
    for (const auto& particle : mcParticles) {
      if (!TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      fillTablesParticle(particle, o2::constants::physics::MassLambdaCPlus);
    }
  }

  void processDataWithDCAFitterN(aod::Collisions const& collisions,
                                 SelectedCandidates const&,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<false, false, false>(collisions, candidatesAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processDataWithDCAFitterN, "Process data", true);

  void processMcWithDCAFitterSig(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcWithDCAFitterSig, "Process MC only for signals", false);

  void processMcWithDCAFitterBkg(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcWithDCAFitterBkg, "Process MC only for background", false);

  void processMcWithDCAFitterAll(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorLcToPKPi, processMcWithDCAFitterAll, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorLcToPKPi>(cfgc)};
}
