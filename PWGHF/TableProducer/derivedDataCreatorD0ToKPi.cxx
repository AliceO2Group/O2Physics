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

/// \file derivedDataCreatorD0ToKPi.cxx
/// \brief Producer of derived tables of D0 candidates, collisions and MC particles
/// \note Based on treeCreatorD0ToKPi.cxx
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
struct HfDerivedDataCreatorD0ToKPi {
  Produces<o2::aod::HfD0Bases> rowCandidateBase;
  Produces<o2::aod::HfD0Pars> rowCandidatePar;
  Produces<o2::aod::HfD0ParEs> rowCandidateParE;
  Produces<o2::aod::HfD0Sels> rowCandidateSel;
  Produces<o2::aod::HfD0Ids> rowCandidateId;
  Produces<o2::aod::HfD0Mcs> rowCandidateMc;
  Produces<o2::aod::HfD0CollBases> rowCollBase;
  Produces<o2::aod::HfD0CollIds> rowCollId;
  Produces<o2::aod::HfD0PBases> rowParticleBase;
  Produces<o2::aod::HfD0PIds> rowParticleId;

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

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPiExt, aod::TracksPidKaExt>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF, aod::HfSelD0>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcKf = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_2prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesKf> candidatesKfPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcKf> candidatesMcKfPerCollision = aod::hf_cand::collisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesKf> candidatesKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcKf> candidatesMcKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> candidatesMcKfSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> candidatesMcKfBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  void init(InitContext const&)
  {
    std::array<bool, 8> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterSig, doprocessMcWithDCAFitterBkg, doprocessMcWithDCAFitterAll, doprocessMcWithKFParticleSig, doprocessMcWithKFParticleBkg, doprocessMcWithKFParticleAll};
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
  auto fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, int candFlag, double invMass, double cosThetaStar, double topoChi2,
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
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct());
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
        topoChi2,
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
        cosThetaStar,
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

  template <int reconstructionType, bool isMc, bool onlyBkg, bool onlySig, typename CandType>
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
            if (TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
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
            if (!TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
              continue;
            }
          }
        } else {
          flagMcRec = 0;
          origin = 0;
        }
        auto prong0 = candidate.template prong0_as<TracksWPid>();
        auto prong1 = candidate.template prong1_as<TracksWPid>();
        double ctD = hfHelper.ctD0(candidate);
        float massD0, massD0bar;
        float topolChi2PerNdf = -999.;
        if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
          massD0 = candidate.kfGeoMassD0();
          massD0bar = candidate.kfGeoMassD0bar();
          topolChi2PerNdf = candidate.kfTopolChi2OverNdf();
        } else {
          massD0 = hfHelper.invMassD0ToPiK(candidate);
          massD0bar = hfHelper.invMassD0barToKPi(candidate);
        }
        if (candidate.isSelD0()) {
          fillTablesCandidate(candidate, prong0, prong1, 0, massD0, hfHelper.cosThetaStarD0(candidate), topolChi2PerNdf, ctD, flagMcRec, origin);
        }
        if (candidate.isSelD0bar()) {
          fillTablesCandidate(candidate, prong0, prong1, 1, massD0bar, hfHelper.cosThetaStarD0bar(candidate), topolChi2PerNdf, ctD, flagMcRec, origin);
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
      if (!TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      fillTablesParticle(particle, o2::constants::physics::MassD0);
    }
  }

  void processDataWithDCAFitterN(aod::Collisions const& collisions,
                                 SelectedCandidates const&,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, false, false>(collisions, candidatesAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithKFParticle(aod::Collisions const& collisions,
                                 SelectedCandidatesKf const&,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, false, false>(collisions, candidatesKfAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  void processMcWithDCAFitterSig(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterSig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterBkg(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterAll(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const&,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithKFParticleSig(aod::Collisions const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  TracksWPid const& tracks,
                                  aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, true>(collisions, candidatesMcKfSig, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleSig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleBkg(aod::Collisions const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  TracksWPid const& tracks,
                                  aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false>(collisions, candidatesMcKfBkg, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleAll(aod::Collisions const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  TracksWPid const& tracks,
                                  aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, false>(collisions, candidatesMcKfAll, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorD0ToKPi>(cfgc)};
}
