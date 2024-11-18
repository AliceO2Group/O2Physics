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
struct HfDerivedDataCreatorD0ToKPi {
  // Candidates
  Produces<o2::aod::HfD0Bases> rowCandidateBase;
  Produces<o2::aod::HfD0Pars> rowCandidatePar;
  Produces<o2::aod::HfD0ParEs> rowCandidateParE;
  Produces<o2::aod::HfD0Sels> rowCandidateSel;
  Produces<o2::aod::HfD0Mls> rowCandidateMl;
  Produces<o2::aod::HfD0Ids> rowCandidateId;
  Produces<o2::aod::HfD0Mcs> rowCandidateMc;
  // Collisions
  Produces<o2::aod::HfD0CollBases> rowCollBase;
  Produces<o2::aod::HfD0CollIds> rowCollId;
  // MC collisions
  Produces<o2::aod::HfD0McCollBases> rowMcCollBase;
  Produces<o2::aod::HfD0McCollIds> rowMcCollId;
  Produces<o2::aod::HfD0McRCollIds> rowMcRCollId;
  // MC particles
  Produces<o2::aod::HfD0PBases> rowParticleBase;
  Produces<o2::aod::HfD0PIds> rowParticleId;

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
  // using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>>;
  using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcKf = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesKfMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesMcKfMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
  using TypeMcCollisions = aod::McCollisions;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_2prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesKf> candidatesKfPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcKf> candidatesMcKfPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesKfMl> candidatesKfMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcKfMl> candidatesMcKfMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesKf> candidatesKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcKf> candidatesMcKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesKfMl> candidatesKfMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> candidatesMcKfSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> candidatesMcKfBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  void init(InitContext const&)
  {
    std::array<bool, 16> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterSig, doprocessMcWithDCAFitterBkg, doprocessMcWithDCAFitterAll, doprocessMcWithKFParticleSig, doprocessMcWithKFParticleBkg, doprocessMcWithKFParticleAll,
                                   doprocessDataWithDCAFitterNMl, doprocessDataWithKFParticleMl, doprocessMcWithDCAFitterMlSig, doprocessMcWithDCAFitterMlBkg, doprocessMcWithDCAFitterMlAll, doprocessMcWithKFParticleMlSig, doprocessMcWithKFParticleMlBkg, doprocessMcWithKFParticleMlAll};
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

  template <typename T>
  void fillTablesCandidate(const T& candidate, int candFlag, double invMass, double cosThetaStar, double topoChi2,
                           double ct, double y, int8_t flagMc, int8_t origin, const std::vector<float>& mlScores)
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
      std::array<std::array<std::array<float, 3>, 2>, 2> sigmas{}; // PID nSigma [Expected][Hypothesis][TPC/TOF/TPC+TOF]
      if (candFlag == 0) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], candidate, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], candidate, 0, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], candidate, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], candidate, 1, Ka)
      } else if (candFlag == 1) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], candidate, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], candidate, 1, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], candidate, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], candidate, 0, Ka)
      }
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
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][2],
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct());
    }
    if (fillCandidateParE) {
      rowCandidateParE(
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
    if (fillCandidateMl) {
      rowCandidateMl(
        mlScores);
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

  template <aod::hf_cand::VertexerType reconstructionType, bool isMl, bool isMc, bool onlyBkg, bool onlySig, typename CollType, typename CandType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
                         aod::Tracks const&,
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
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
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
            if (TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
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
            if (!TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
              continue;
            }
          }
        } else {
          if (downSampleBkgFactor < 1.) {
            float pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
            if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
              continue;
            }
          }
        }

        double ct = hfHelper.ctD0(candidate);
        double y = hfHelper.yD0(candidate);
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
        std::vector<float> mlScoresD0, mlScoresD0bar;
        if constexpr (isMl) {
          std::copy(candidate.mlProbD0().begin(), candidate.mlProbD0().end(), std::back_inserter(mlScoresD0));
          std::copy(candidate.mlProbD0bar().begin(), candidate.mlProbD0bar().end(), std::back_inserter(mlScoresD0bar));
        }
        if (candidate.isSelD0()) {
          fillTablesCandidate(candidate, 0, massD0, hfHelper.cosThetaStarD0(candidate), topolChi2PerNdf, ct, y, flagMcRec, origin, mlScoresD0);
        }
        if (candidate.isSelD0bar()) {
          fillTablesCandidate(candidate, 1, massD0bar, hfHelper.cosThetaStarD0bar(candidate), topolChi2PerNdf, ct, y, flagMcRec, origin, mlScoresD0bar);
        }
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
        fillTablesParticle(particle, o2::constants::physics::MassD0);
      }
    }
  }

  void processDataWithDCAFitterN(CollisionsWCentMult const& collisions,
                                 SelectedCandidates const&,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, false, false, false>(collisions, candidatesAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithKFParticle(CollisionsWCentMult const& collisions,
                                 SelectedCandidatesKf const&,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, false, false, false>(collisions, candidatesKfAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  void processMcWithDCAFitterSig(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterSig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterBkg(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterAll(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithKFParticleSig(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, false, true>(collisions, candidatesMcKfSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleSig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleBkg(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, true, false>(collisions, candidatesMcKfBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleAll(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, false, false>(collisions, candidatesMcKfAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);

  // ML versions

  void processDataWithDCAFitterNMl(CollisionsWCentMult const& collisions,
                                   SelectedCandidatesMl const&,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterNMl, "Process data with DCAFitterN and ML", false);

  void processDataWithKFParticleMl(CollisionsWCentMult const& collisions,
                                   SelectedCandidatesKfMl const&,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, false, false>(collisions, candidatesKfMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticleMl, "Process data with KFParticle and ML", false);

  void processMcWithDCAFitterMlSig(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false, true>(collisions, candidatesMcMlSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlSig, "Process MC with DCAFitterN and ML only for signals", false);

  void processMcWithDCAFitterMlBkg(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, true, false>(collisions, candidatesMcMlBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlBkg, "Process MC with DCAFitterN and ML only for background", false);

  void processMcWithDCAFitterMlAll(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false, false>(collisions, candidatesMcMlAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlAll, "Process MC with DCAFitterN and ML", false);

  void processMcWithKFParticleMlSig(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false, true>(collisions, candidatesMcKfMlSig, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlSig, "Process MC with KFParticle and ML only for signals", false);

  void processMcWithKFParticleMlBkg(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, true, false>(collisions, candidatesMcKfMlBkg, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlBkg, "Process MC with KFParticle and ML only for background", false);

  void processMcWithKFParticleMlAll(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    preProcessMcCollisions(mcCollisions, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false, false>(collisions, candidatesMcKfMlAll, tracks, bcs);
    processMcParticles(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlAll, "Process MC with KFParticle and ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorD0ToKPi>(cfgc)};
}
