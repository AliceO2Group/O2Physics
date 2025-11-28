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

/// \file utilsCharmReso.h
/// \brief utility functions for charm-hadron resonance derived data creators
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#ifndef PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_
#define PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Utils/utilsMcMatching.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>

#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace o2::analysis {
  namespace hf_charm_reso {

    // event types
    enum Event : uint8_t {
      Processed = 0,
      NoDV0Selected,
      DV0Selected,
      kNEvent
    };

    enum BachelorType : uint8_t {
      K0s = 0,
      Lambda,
      AntiLambda,
      Track
    };

    enum DType : uint8_t {
      Dplus = 1,
      Dstar,
      D0
    };

    enum PairingType : uint8_t {
      V0Only,
      TrackOnly,
      V0AndTrack
    };

    enum D0Sel : uint8_t {
      SelectedD0 = 0,
      SelectedD0Bar
    };

    // Helper structs to pass D and V0 informations

    struct candidateV0 {
      std::array<float, 3> pos;
      std::array<float, 3> mom;
      std::array<float, 3> momPos;
      std::array<float, 3> momNeg;
      float pT;
      float cosPA;
      float dcaV0ToPv;
      float dcaDau;
      float alpha;
      float eta;
      float radius;
      float mK0Short;
      float mLambda;
      uint8_t v0Type;
    };

    struct varContainer {
      float invMassD;
      float ptD;
      float invMassD0;
      float invMassD0Bar;
      float invMassReso;
      float ptReso;
      int8_t signD;
      std::array<float, 3> pVectorProng0;
      std::array<float, 3> pVectorProng1;
      std::array<float, 3> pVectorProng2;
    };

    /// Basic track quality selections for V0 daughters
    /// \param Tr is a track
    /// \param dDaughtersIds are the IDs of the D meson daughter tracks
    /// \param cfgV0Cuts are the cuts to be applied to the V0
    /// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
    template <typename Tr, typename Cuts>
    bool selectV0Daughter(Tr const& track, const std::array<int, 3>& dDaughtersIds, const Cuts& cfgV0Cuts, bool rejectPairsWithCommonDaughter)
    {
      // acceptance selection
      if (std::abs(track.eta()) > cfgV0Cuts.etaMaxDau.value) {
        return false;
      }
      // Tpc Refit
      if (!(track.hasTPC())) {
        return false;
      }
      // track quality selection
      if (track.itsNCls() < cfgV0Cuts.trackNclusItsCut.value ||
          track.tpcNClsFound() < cfgV0Cuts.trackNCrossedRowsTpc.value ||
          track.tpcNClsCrossedRows() < cfgV0Cuts.trackNCrossedRowsTpc.value ||
          track.tpcNClsCrossedRows() < cfgV0Cuts.trackFracMaxindableTpcCls.value * track.tpcNClsFindable() ||
          track.tpcNClsShared() > cfgV0Cuts.trackNsharedClusTpc.value) {
        return false;
      }
      // rejection of tracks that share a daughter with the D meson
      if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
        return false;
      }
      return true;
    }

    /// Utility to find which v0 daughter carries the largest fraction of the mother longitudinal momentum
    /// \param momV0 is the momentum of the V0
    /// \param momDau0 is the momentum of first daughter
    /// \param momDau1 is the momentum of second daughter
    /// \return alphaAP
    float alphaAP(std::array<float, 3> const& momV0, std::array<float, 3> const& momDau0, std::array<float, 3> const& momDau1)
    {
      float const momTot = std::sqrt(std::pow(momV0[0], 2.) + std::pow(momV0[1], 2.) + std::pow(momV0[2], 2.));
      float const lQlPos = (momDau0[0] * momV0[0] + momDau0[1] * momV0[1] + momDau0[2] * momV0[2]) / momTot;
      float const lQlNeg = (momDau1[0] * momV0[0] + momDau1[1] * momV0[1] + momDau1[2] * momV0[2]) / momTot;
      return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
    }

    /// Utility to find DCA of V0 to Primary vertex
    /// \param x is the x-coordinate
    /// \param y is the y-coordinate
    /// \param z is the z-coordinate
    /// \param px is the x-component of the momentum
    /// \param py is the y-component of the momentum
    /// \param pz is the z-component of the momentum
    /// \param pvX is the x-coordinate of the PV
    /// \param pvY is the y-coordinate of the PV
    /// \param pvZ is the z-coordinate of the PV
    /// \return the DCA
    float calculateDCAStraightToPV(float x, float y, float z, float px, float py, float pz, float pvX, float pvY, float pvZ)
    {
      return std::sqrt((std::pow((pvY - y) * pz - (pvZ - z) * py, 2) + std::pow((pvX - x) * pz - (pvZ - z) * px, 2) + std::pow((pvX - x) * py - (pvY - y) * px, 2)) / (px * px + py * py + pz * pz));
    }

    /// Basic selection of V0 candidates
    /// \param collision is the current collision
    /// \param dauTracks are the v0 daughter tracks
    /// \param dDaughtersIds are the IDs of the D meson daughter tracks
    /// \param fitter is the DCAFitter object
    /// \param cfgV0Cuts are the cuts to be applied to the V0
    /// \param v0 is the V0 candidate
    /// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
    /// \return a bitmap with mass hypotesis if passes all cuts
    template <typename Coll, typename Tr, typename Cuts>
    bool buildAndSelectV0(const Coll& collision, const std::array<int, 3>& dDaughtersIds, const std::array<Tr, 2>& dauTracks, const Cuts& cfgV0Cuts, o2::vertexing::DCAFitterN<2>& fitter, candidateV0& v0, bool rejectPairsWithCommonDaughter)
    {
      const auto& trackPos = dauTracks[0];
      const auto& trackNeg = dauTracks[1];
      // single-tracks selection
      if (!selectV0Daughter(trackPos, dDaughtersIds, cfgV0Cuts, rejectPairsWithCommonDaughter) || !selectV0Daughter(trackNeg, dDaughtersIds, cfgV0Cuts, rejectPairsWithCommonDaughter)) {
        return false;
      }
      // daughters DCA to V0's collision primary vertex
      std::array<float, 2> dcaInfo{};
      auto trackPosPar = getTrackPar(trackPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPosPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto trackPosDcaXY = dcaInfo[0];
      auto trackNegPar = getTrackPar(trackNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackNegPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto trackNegDcaXY = dcaInfo[0];
      if (std::fabs(trackPosDcaXY) < cfgV0Cuts.dcaMaxDauToPv.value || std::fabs(trackNegDcaXY) < cfgV0Cuts.dcaMaxDauToPv.value) {
        return false;
      }
      // vertex reconstruction
      auto trackPosCov = getTrackParCov(trackPos);
      auto trackNegCov = getTrackParCov(trackNeg);
      int nCand = 0;
      try {
        nCand = fitter.process(trackPosCov, trackNegCov);
      } catch (...) {
        return false;
      }
      if (nCand == 0) {
        return false;
      }
      // compute candidate momentum from tracks propagated to decay vertex
      auto& trackPosProp = fitter.getTrack(0);
      auto& trackNegProp = fitter.getTrack(1);
      trackPosProp.getPxPyPzGlo(v0.momPos);
      trackNegProp.getPxPyPzGlo(v0.momNeg);

      v0.mom = RecoDecay::pVec(v0.momPos, v0.momNeg);

      v0.pT = std::hypot(v0.mom[0], v0.mom[1]);
      // topological selections:
      // v0 eta
      v0.eta = RecoDecay::eta(v0.mom);
      if (std::abs(v0.eta) > cfgV0Cuts.etaMax.value) {
        return false;
      }
      // daughters DCA
      v0.dcaDau = std::sqrt(fitter.getChi2AtPCACandidate());
      if (v0.dcaDau > cfgV0Cuts.dcaDau.value) {
        return false;
      }
      // v0 radius
      const auto& vtx = fitter.getPCACandidate();
      v0.radius = std::hypot(vtx[0], vtx[1]);
      if (v0.radius < cfgV0Cuts.radiusMin.value) {
        return false;
      }
      std::copy(vtx.begin(), vtx.end(), v0.pos.begin());

      // v0 DCA to primary vertex
      v0.dcaV0ToPv = calculateDCAStraightToPV(
        vtx[0], vtx[1], vtx[2],
        v0.momPos[0] + v0.momNeg[0],
        v0.momPos[1] + v0.momNeg[1],
        v0.momPos[2] + v0.momNeg[2],
        collision.posX(), collision.posY(), collision.posZ());
      if (std::abs(v0.dcaV0ToPv) > cfgV0Cuts.dcaPv.value) {
        return false;
      }
      // v0 cosine of pointing angle
      std::array<float, 3> const primVtx = {collision.posX(), collision.posY(), collision.posZ()};
      v0.cosPA = RecoDecay::cpa(primVtx, vtx, v0.mom);
      if (v0.cosPA < cfgV0Cuts.cosPa.value) {
        return false;
      }
      // distinguish between K0s, and Lambda hypotesys
      v0.v0Type = {BIT(BachelorType::K0s) | BIT(BachelorType::Lambda) | BIT(BachelorType::AntiLambda)};
      // for lambda hypotesys define if its lambda or anti-lambda
      v0.alpha = alphaAP(v0.mom, v0.momPos, v0.momNeg);
      bool const matter = v0.alpha > 0;
      CLRBIT(v0.v0Type, matter ? BachelorType::AntiLambda : BachelorType::Lambda);
      auto massPos = matter ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged;
      auto massNeg = matter ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
      // mass hypotesis
      v0.mLambda = RecoDecay::m(std::array{v0.momPos, v0.momNeg}, std::array{massPos, massNeg});
      v0.mK0Short = RecoDecay::m(std::array{v0.momPos, v0.momNeg}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
      if (std::fabs(v0.mK0Short - o2::constants::physics::MassK0) > cfgV0Cuts.deltaMassK0s.value) {
        CLRBIT(v0.v0Type, BachelorType::K0s);
      }
      if (std::fabs(v0.mLambda - o2::constants::physics::MassLambda0) > cfgV0Cuts.deltaMassLambda.value) {
        CLRBIT(v0.v0Type, BachelorType::Lambda);
        CLRBIT(v0.v0Type, BachelorType::AntiLambda);
      }
      // PID
      if (TESTBIT(v0.v0Type, BachelorType::K0s)) {
        if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value) ||
            (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value)) {
          CLRBIT(v0.v0Type, BachelorType::K0s);
        }
      }
      if (TESTBIT(v0.v0Type, BachelorType::Lambda)) {
        if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc.value) ||
            (trackPos.hasTOF() && std::fabs(trackPos.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr.value) ||
            (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value)) {
          CLRBIT(v0.v0Type, BachelorType::Lambda);
        }
      }
      if (TESTBIT(v0.v0Type, BachelorType::AntiLambda)) {
        if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value) ||
            (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc.value) ||
            (trackNeg.hasTOF() && std::fabs(trackNeg.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr.value)) {
          CLRBIT(v0.v0Type, BachelorType::AntiLambda);
        }
      }
      if (v0.v0Type == 0) {
        return false;
      }
      return true;
    }

    /// Basic selection of tracks
    /// \param track is the track
    /// \param dDaughtersIds are the IDs of the D meson daughter tracks
    /// \param cfgV0Cuts are the cuts to be applied to the track
    /// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
    /// \return true if passes all cuts
    template <typename Tr, typename Cuts>
    bool isTrackSelected(const Tr& track, const std::array<int, 3>& dDaughtersIds, const Cuts& cfgSingleTrackCuts, bool rejectPairsWithCommonDaughter)
    {
      if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
        return false;
      }
      switch (cfgSingleTrackCuts.setTrackSelections.value) {
        case 1:
          if (!track.isGlobalTrackWoDCA()) {
            return false;
          }
          break;
        case 2:
          if (!track.isGlobalTrack()) {
            return false;
          }
          break;
        case 3:
          if (!track.isQualityTrackITS()) {
            return false;
          }
          break;
      }
      if (track.pt() < cfgSingleTrackCuts.minPt.value) {
        return false;
      }
      if (std::abs(track.eta()) > cfgSingleTrackCuts.maxEta.value) {
        return false;
      }
      if (!track.hasTPC()) {
        return false;
      }
      bool const isPion = std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value;
      bool const isKaon = std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value;
      bool const isProton = std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value;
      return (isPion || isKaon || isProton); // we keep the track if is it compatible with at least one of the PID hypotheses selected
    }

    /// Matching of V0 candidates to MC truth
    /// \param particlesMc is the table of MC particles
    /// \param arrDaughtersV0 is the array of V0 daughter tracks
    /// \return the MC matching flag for the V0
    template <typename PParticles, typename TrIU>
    int8_t getMatchingFlagV0(PParticles const& particlesMc, const std::array<TrIU, 2>& arrDaughtersV0)
    {
      int8_t signV0{0};
      int indexRec{-1};
      int flagV0{0};
      indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersV0, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2);
      if (indexRec > -1) {
        flagV0 = hf_decay::hf_cand_reso::PartialMatchMc::K0Matched;
      } else {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersV0, kLambda0, std::array{+kProton, -kPiPlus}, true, &signV0, 2);
        if (indexRec > -1) {
          flagV0 = signV0 * hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched;
        }
      }
      return flagV0; // Placeholder, should return the actual flag based on matching logic
    }

    /// Matching of V0 candidates to MC truth
    /// \param bachTrack is the track
    /// \return the MC matching flag for the track
    template <typename Tr>
    int8_t getMatchingFlagTrack(Tr const& bachTrack)
    {
      auto particle = bachTrack.mcParticle();
      auto pdgCode = std::abs(particle.pdgCode());
      if (pdgCode == kPiPlus) {
        return hf_decay::hf_cand_reso::PartialMatchMc::PionMatched;
      }
      if (pdgCode == kKPlus) {
        return hf_decay::hf_cand_reso::PartialMatchMc::KaonMatched;
      }
      if (pdgCode == kProton) {
        return hf_decay::hf_cand_reso::PartialMatchMc::ProtonMatched;
      }
      return 0;
    }

    /// Matching of V0 candidates to MC truth
    /// \param particlesMc is the table of MC particles
    /// \param indexRec is the index of the MC particle associated to the reconstructed canditate
    /// \param pdg is the O2DatabasePDG service
    /// \return the generated invariant mass
    template <typename PParticles>
    float computeInvMassGen(PParticles const& particlesMc, int indexRec, o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg)
    {
      if (indexRec < 0) {
        return -1.f;
      }
      auto particleReso = particlesMc.iteratorAt(indexRec);
      auto dau1 = particlesMc.iteratorAt(particleReso.daughtersIds().front());
      auto dau2 = particlesMc.iteratorAt(particleReso.daughtersIds().back());
      std::array<std::array<float, 3>, 2> pArr = {{{dau1.px(), dau1.py(), dau1.pz()}, {dau2.px(), dau2.py(), dau2.pz()}}};
      std::array<float, 2> mArr = {static_cast<float>(pdg->Mass(dau1.pdgCode())), static_cast<float>(pdg->Mass(dau2.pdgCode()))};
      return static_cast<float>(RecoDecay::m(pArr, mArr));
    }

    /// Function for filling MC reco information of DV0 candidates in the tables
    /// \tparam dType is the D meson type (Dstar, Dplus or D0)
    /// \param particlesMc is the table with MC particles
    /// \param candCharmBach is the D meson candidate
    /// \param bachelorV0 is the V0 candidate
    /// \param tracks is the table with tracks
    /// \param indexHfCandCharm is the index of the charm-hadron bachelor in the reduced table
    /// \param indexCandV0TrBach is the index of the v0 bachelor in the reduced table
    /// \param pdg is the O2DatabasePDG service
    /// \param registry is the histogram registry
    /// \param rowMcRecReduced is the table to be filled
    template <uint8_t DType, typename PParticles, typename CCand, typename BBachV0, typename Tr, typename Table>
    void fillMcRecoInfoDV0(PParticles const& particlesMc,
                          CCand const& candCharmBach,
                          BBachV0 const& bachelorV0,
                          Tr const& tracks,
                          int64_t& indexHfCandCharm,
                          int64_t& indexCandV0Bach,
                          o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                          o2::framework::HistogramRegistry& registry,
                          Table& rowMcRecReduced)
    {
      std::vector<typename Tr::iterator> vecDaughtersReso{};
      int8_t sign{0}, nKinkedTracks{0}, origin{0}, flagCharmBach{0}, flagCharmBachInterm{0}, flagV0{0}, flagReso{0};
      int indexRec{-1}, debugMcRec{0};
      float ptGen{-1.f}, invMassGen{-1.f};
      if constexpr (DType == DType::Dstar) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prongPiId()));
        // Check if D* is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DstarMatched);
          origin = candCharmBach.originMcRec();
        }
        // Check if D0 is matched
        flagCharmBachInterm = candCharmBach.flagMcMatchRecD0();
        if (flagCharmBachInterm != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
        }
        // Check if V0 is matched
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
        flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]});
        if (flagV0 != 0) {
          SETBIT(debugMcRec, std::abs(flagV0));
        }
        // If both D* and K0s are matched, try to match resonance
        if (flagCharmBach != 0 && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
          std::array<int, 5> const pdgCodesDaughters = {+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus};
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
        } else if (flagCharmBachInterm != 0 && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
          std::array<int, 4> const pdgCodesDaughters = {+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus};
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[3], vecDaughtersReso[4]};
          // Peaking background of D0K0s <- Ds* with spurious soft pion
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
              break;
            }
          }
        }
        // No physical channel expected in D*Lambda
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      } else if constexpr (DType == DType::Dplus) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong2Id()));
        // Check if D+ is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DplusMatched);
          origin = candCharmBach.originMcRec();
        }
        // Check if V0 is matched
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
        flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]});
        if (flagV0 != 0) {
          SETBIT(debugMcRec, std::abs(flagV0));
        }
        // If both D+ and K0s are matched, try to match resonance
        if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
          auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], +kPiPlus, -kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusK0s) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
          // Partial matching of Dsj -> D*K0s -> (D+ pi0) (K0s) with missing neutral
          if (indexRec < 0) {
            for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
              indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
              if (indexRec > -1) {
                flagReso = sign * decayChannelFlag;
                SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
                break;
              }
            }
          }

        } else if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && std::abs(flagV0) == hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched) {
          // Peaking background of D+Lambda <- Ds* with spurious soft pion
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
          auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], +kProton, -kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusLambda) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
        }
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      } else if constexpr (DType == DType::D0) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        // Check if D0 is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
          origin = candCharmBach.originMcRec();
        }
        // Check if V0 is matched
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
        flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[2], vecDaughtersReso[3]});
        if (flagV0 != 0) {
          SETBIT(debugMcRec, std::abs(flagV0));
        }
        // No physical channel expected in D0 K0s
        // If both D0 and Lambda are matched, try to match resonance
        if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && std::abs(flagV0) == hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3]};
          auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kProton, -kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Lambda) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
        }
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      }
      registry.fill(HIST("hMCRecDebug"), debugMcRec);
      if (indexRec > -1) {
        registry.fill(HIST("hMCRecCounter"), flagReso);
        registry.fill(HIST("hMCRecOrigin"), origin);
        registry.fill(HIST("hMCRecMassGen"), invMassGen);
      }
      if (flagCharmBach != 0) {
        registry.fill(HIST("hMCRecCharmDau"), flagCharmBach);
      }
    }

    // Function for filling MC reco information of D Track candidates in the tables
    /// \tparam dType is the D meson type (Dstar, Dplus or D0)
    /// \param particlesMc is the table with MC particles
    /// \param candCharmBach is the D meson candidate
    /// \param bachelorTrack is the bachelor track
    /// \param tracks is the table with tracks
    /// \param indexHfCandCharm is the index of the charm-hadron bachelor in the reduced table
    /// \param indexCandTrBach is the index of the v0 bachelor in the reduced table
    /// \param pdg is the O2DatabasePDG service
    /// \param registry is the histogram registry
    /// \param rowMcRecReduced is the table to be filled
    template <uint8_t DType, typename PParticles, typename CCand, typename BBachTr, typename Tr, typename Table>
    void fillMcRecoInfoDTrack(PParticles const& particlesMc,
                              CCand const& candCharmBach,
                              BBachTr const& bachelorTrack,
                              Tr const& tracks,
                              const int64_t indexHfCandCharm,
                              const int64_t indexCandTrBach,
                              o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                              o2::framework::HistogramRegistry& registry,
                              Table& rowMcRecReduced)
    {
      std::vector<typename Tr::iterator> vecDaughtersReso{};
      int8_t sign{0}, nKinkedTracks{0}, origin{0}, flagCharmBach{0}, flagCharmBachInterm{0}, flagTrack{0}, flagReso{0};
      int indexRec{-1};
      uint16_t debugMcRec{0};
      float ptGen{-1.f}, invMassGen{-1.f};
      if constexpr (DType == DType::Dstar) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prongPiId()));
        // Check if D* is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DstarMatched);
          origin = candCharmBach.originMcRec();
        }
        // Check if D0 is matched
        flagCharmBachInterm = candCharmBach.flagMcMatchRecD0();
        if (flagCharmBachInterm != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
        }
        // Check if Track is matched
        flagTrack = getMatchingFlagTrack(bachelorTrack);
        if (flagTrack != 0) {
          SETBIT(debugMcRec, flagTrack);
        }
        // If both D* and Track are matched, try to match resonance
        if (flagCharmBach != 0 && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], bachelorTrack};
          auto pdgCodesDaughters = std::array{+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
        }
        // No channels in D*K+ or D*Pr
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      } else if constexpr (DType == DType::Dplus) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong2Id()));
        // Check if D+ is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DplusMatched);
          origin = candCharmBach.originMcRec();
        }
        // Check if Track is matched
        flagTrack = getMatchingFlagTrack(bachelorTrack);
        if (flagTrack != 0) {
          SETBIT(debugMcRec, flagTrack);
        }
        // If both D+ and Track are matched, try to match resonance
        if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], bachelorTrack};
          auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], -kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusPi) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
          // Partial matching of Dj -> D*Pi -> (D+ pi0) (pi) with missing neutral
          if (indexRec < 0) {
            for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
              indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
              if (indexRec > -1) {
                flagReso = sign * decayChannelFlag;
                SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
                break;
              }
            }
          }
        }
        // No channels in D+K+ or D+Pr
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      } else if constexpr (DType == DType::D0) {
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
        vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
        // Check if D0 is matched
        flagCharmBach = candCharmBach.flagMcMatchRec();
        flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
        if (flagCharmBach != 0) {
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
          origin = candCharmBach.originMcRec();
        }
        flagTrack = getMatchingFlagTrack(bachelorTrack);
        if (flagTrack != 0) {
          SETBIT(debugMcRec, flagTrack);
        }
        if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], bachelorTrack};
          auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kPiPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Pi) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
          // Partial matching of Dj -> D*Pi -> (D0 pi) (pi) with missing pion
          if (indexRec < 0) {
            for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
              indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
              if (indexRec > -1) {
                flagReso = sign * decayChannelFlag;
                SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
                break;
              }
            }
          }
        } else if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::KaonMatched) {
          auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], bachelorTrack};
          auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
          auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kKPlus};
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Kplus) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
            if (indexRec > -1) {
              flagReso = sign * decayChannelFlag;
              break;
            }
          }
        }
        if (indexRec > -1) {
          auto particleReso = particlesMc.iteratorAt(indexRec);
          ptGen = particleReso.pt();
          invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
        }
        rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                        flagReso, flagCharmBach,
                        flagCharmBachInterm, debugMcRec,
                        origin, ptGen, invMassGen,
                        nKinkedTracks);
      }
      registry.fill(HIST("hMCRecDebug"), debugMcRec);
      if (indexRec > -1) {
        registry.fill(HIST("hMCRecCounter"), flagReso);
        registry.fill(HIST("hMCRecOrigin"), origin);
        registry.fill(HIST("hMCRecMassGen"), invMassGen);
      }
      if (flagCharmBach != 0) {
        registry.fill(HIST("hMCRecCharmDau"), flagCharmBach);
      }
    } // fillMcRecoInfoDTrack

    // Function for derived data creation
    /// \tparam dType is the D meson type (Dstar, Dplus or D0)
    /// \param collision is the current collision
    /// \param candsD are the D meson candidates in the current collision
    /// \param bachelorV0s are the V0 candidates in the current collision
    /// \param bachelorTrks are the track ids in the current collision
    /// \param tracks is the track table
    /// \param tracksIU is the trackIU table
    /// \param particlesMc is the MC particle table
    /// \param hfRejMap is the event rejection map from the HF event selection util
    /// \param pdg is the O2DatabasePDG service
    /// \param registry is the histogram registry
    /// \param matCorr is the material correction type to be used in the track propagation
    /// \param fitter is the DCAFitter object
    /// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
    /// \param rowCollisionReduced is the collision reduced table to be filled
    /// \param rowCandDmesReduced is the D-meson reduced table to be filled
    /// \param rowCandV0Reduced is the V0 reduced table to be filled
    /// \param rowTrkReduced is the track reduced table to be filled
    /// \param rowMcRecV0Reduced is the MC reco D-V0 reduced table to be filled
    /// \param rowMcRecTrkReduced is the MC reco D-track reduced table to be filled
    /// \param rowCandDmesMlReduced is the ML reduced table to be filled
    template <bool WithMl, bool DoMc, uint8_t DType, uint8_t PairingType, typename Coll, typename CCands, typename Tr, typename TrIU, typename PParticles, typename BBachV0s, typename BBachTracks, typename DmesCuts, typename TrkCuts, typename V0Cuts, typename QaConfig, typename TableCollRed, typename TableCandDRed, typename TableCandV0Red, typename TableTrkRed, typename TableMcRecV0Red, typename TableMcRecTrkRed, typename TableCandDMlRed>
    void runDataCreation(Coll const& collision,
                         CCands const& candsD,
                         BBachV0s const& bachelorV0s,
                         BBachTracks const& bachelorTrks,
                         Tr const& tracks,
                         TrIU const& tracksIU,
                         PParticles const& particlesMc,
                         o2::hf_evsel::HfCollisionRejectionMask hfRejMap,
                         const float bz,
                         o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                         o2::framework::HistogramRegistry& registry,
                         o2::base::Propagator::MatCorrType const& matCorr,
                         o2::vertexing::DCAFitterN<2>& fitter,
                         DmesCuts const& cfgDmesCuts,
                         TrkCuts const& cfgSingleTrackCuts,
                         V0Cuts const& cfgV0Cuts,
                         QaConfig const& cfgQaPlots,
                         bool rejectPairsWithCommonDaughter,
                         TableCollRed& rowCollisionReduced,
                         TableCandDRed& rowCandDmesReduced,
                         TableCandV0Red& rowCandV0Reduced,
                         TableTrkRed& rowTrkReduced,
                         TableMcRecV0Red& rowMcRecV0Reduced,
                         TableMcRecTrkRed& rowMcRecTrkReduced,
                         TableCandDMlRed& rowCandDmesMlReduced)
    {
      int const indexHfReducedCollision = rowCollisionReduced.lastIndex() + 1;
      // std::map where the key is the V0.globalIndex() and
      // the value is the V0 index in the table of the selected v0s
      std::map<int64_t, int64_t> selectedV0s;
      std::map<int64_t, int64_t> selectedTracks;
      bool fillHfReducedCollision = false;
      constexpr bool DoTracks = PairingType == PairingType::TrackOnly || PairingType == PairingType::V0AndTrack;
      constexpr bool DoV0s = PairingType == PairingType::V0Only || PairingType == PairingType::V0AndTrack;
      // loop on D candidates
      for (const auto& candD : candsD) {
        // initialize variables depending on D meson type
        bool fillHfCandD = false;
        std::array<float, 3> secondaryVertexD{};
        std::array<int, 3> prongIdsD{};
        std::array<float, 6> bdtScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
        std::vector<std::decay_t<typename TrIU::iterator>> charmHadDauTracks{};
        varContainer varUtils;
        varUtils.ptD = candD.pt();
        if constexpr (DType == DType::Dstar) {
          varUtils.signD = candD.signSoftPi();
          if (varUtils.signD > 0) {
            varUtils.invMassD = candD.invMassDstar();
            varUtils.invMassD0 = candD.invMassD0();
          } else {
            varUtils.invMassD = candD.invMassAntiDstar();
            varUtils.invMassD0 = candD.invMassD0Bar();
          }
          secondaryVertexD[0] = candD.xSecondaryVertexD0();
          secondaryVertexD[1] = candD.ySecondaryVertexD0();
          secondaryVertexD[2] = candD.zSecondaryVertexD0();
          prongIdsD[0] = candD.prong0Id();
          prongIdsD[1] = candD.prong1Id();
          prongIdsD[2] = candD.prongPiId();
          varUtils.pVectorProng0 = candD.pVectorProng0();
          varUtils.pVectorProng1 = candD.pVectorProng1();
          varUtils.pVectorProng2 = candD.pVecSoftPi();
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
          if constexpr (WithMl) {
            std::copy(candD.mlProbDstarToD0Pi().begin(), candD.mlProbDstarToD0Pi().end(), bdtScores.begin());
          }
          registry.fill(HIST("hMassVsPtDstarAll"), varUtils.ptD, varUtils.invMassD - varUtils.invMassD0);
        } else if constexpr (DType == DType::Dplus) {
          auto prong0 = tracksIU.rawIteratorAt(candD.prong0Id());
          varUtils.invMassD = HfHelper::invMassDplusToPiKPi(candD);
          secondaryVertexD[0] = candD.xSecondaryVertex();
          secondaryVertexD[1] = candD.ySecondaryVertex();
          secondaryVertexD[2] = candD.zSecondaryVertex();
          prongIdsD[0] = candD.prong0Id();
          prongIdsD[1] = candD.prong1Id();
          prongIdsD[2] = candD.prong2Id();
          varUtils.signD = prong0.sign();
          varUtils.pVectorProng0 = candD.pVectorProng0();
          varUtils.pVectorProng1 = candD.pVectorProng1();
          varUtils.pVectorProng2 = candD.pVectorProng2();
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong2Id()));
          if constexpr (WithMl) {
            std::copy(candD.mlProbDplusToPiKPi().begin(), candD.mlProbDplusToPiKPi().end(), bdtScores.begin());
          }
          registry.fill(HIST("hMassVsPtDplusAll"), varUtils.ptD, varUtils.invMassD);
        } else if constexpr (DType == DType::D0) {
          varUtils.invMassD0 = HfHelper::invMassD0ToPiK(candD);
          varUtils.invMassD0Bar = HfHelper::invMassD0barToKPi(candD);
          secondaryVertexD[0] = candD.xSecondaryVertex();
          secondaryVertexD[1] = candD.ySecondaryVertex();
          secondaryVertexD[2] = candD.zSecondaryVertex();
          prongIdsD[0] = candD.prong0Id();
          prongIdsD[1] = candD.prong1Id();
          prongIdsD[2] = -1; // D0 does not have a third prong
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
          charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
          varUtils.pVectorProng0 = candD.pVectorProng0();
          varUtils.pVectorProng1 = candD.pVectorProng1();
          varUtils.pVectorProng2 = {0.f, 0.f, 0.f}; // D0 does not have a third prong
          if constexpr (WithMl) {
            std::copy(candD.mlProbD0().begin(), candD.mlProbD0().end(), bdtScores.begin());
            std::copy(candD.mlProbD0bar().begin(), candD.mlProbD0bar().end(), bdtScores.begin() + 3);
          }
          if (candD.isSelD0() >= cfgDmesCuts.selectionFlagD0.value) {
            registry.fill(HIST("hMassVsPtD0All"), varUtils.ptD, varUtils.invMassD0);
          }
          if (candD.isSelD0bar() >= cfgDmesCuts.selectionFlagD0Bar.value) {
            registry.fill(HIST("hMassVsPtD0BarAll"), varUtils.ptD, varUtils.invMassD0Bar);
          }
        } // end of dType switch

        // Get single track variables
        float chi2TpcDauMax = -1.f;
        int nItsClsDauMin = 8, nTpcCrossRowsDauMin = 200;
        float chi2TpcSoftPi = -1.f;
        int nItsClsSoftPi = 8, nTpcCrossRowsSoftPi = 200;
        for (const auto& charmHadTrack : charmHadDauTracks) {
          if (charmHadTrack.itsNCls() < nItsClsDauMin) {
            nItsClsDauMin = charmHadTrack.itsNCls();
          }
          if (charmHadTrack.tpcNClsCrossedRows() < nTpcCrossRowsDauMin) {
            nTpcCrossRowsDauMin = charmHadTrack.tpcNClsCrossedRows();
          }
          if (charmHadTrack.tpcChi2NCl() > chi2TpcDauMax) {
            chi2TpcDauMax = charmHadTrack.tpcChi2NCl();
          }
        }
        if constexpr (DType == DType::Dstar) {
          auto softPi = tracksIU.rawIteratorAt(candD.prongPiId());
          nItsClsSoftPi = softPi.itsNCls();
          nTpcCrossRowsSoftPi = softPi.tpcNClsCrossedRows();
          chi2TpcSoftPi = softPi.tpcChi2NCl();
          charmHadDauTracks.push_back(softPi);
        }
        // Loop on the bachelor V0s
        if constexpr (DoV0s) {
          for (const auto& v0 : bachelorV0s) {
            auto trackPos = tracksIU.rawIteratorAt(v0.posTrackId());
            auto trackNeg = tracksIU.rawIteratorAt(v0.negTrackId());
            // Apply selsection
            auto v0DauTracks = std::array{trackPos, trackNeg};
            candidateV0 candV0;
            if (!buildAndSelectV0(collision, prongIdsD, v0DauTracks, cfgV0Cuts, fitter, candV0, rejectPairsWithCommonDaughter)) {
              continue;
            }
            // Get single track variables
            float chi2TpcDauV0Max = -1.f;
            int nItsClsDauV0Min = 8, nTpcCrossRowsDauV0Min = 200;
            for (const auto& v0Track : v0DauTracks) {
              if (v0Track.itsNCls() < nItsClsDauV0Min) {
                nItsClsDauV0Min = v0Track.itsNCls();
              }
              if (v0Track.tpcNClsCrossedRows() < nTpcCrossRowsDauV0Min) {
                nTpcCrossRowsDauV0Min = v0Track.tpcNClsCrossedRows();
              }
              if (v0Track.tpcChi2NCl() > chi2TpcDauV0Max) {
                chi2TpcDauV0Max = v0Track.tpcChi2NCl();
              }
            }
            // propagate V0 to primary vertex (if enabled)
            if (cfgV0Cuts.propagateV0toPV.value) {
              std::array<float, 3> const pVecV0Orig = {candV0.mom[0], candV0.mom[1], candV0.mom[2]};
              std::array<float, 2> dcaInfo{};
              auto trackParK0 = o2::track::TrackPar(candV0.pos, pVecV0Orig, 0, true);
              trackParK0.setPID(o2::track::PID::K0);
              trackParK0.setAbsCharge(0);
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParK0, 2.f, matCorr, &dcaInfo);
              getPxPyPz(trackParK0, candV0.mom);
            }
            // compute resonance invariant mass and filling of QA histograms
            if (TESTBIT(candV0.v0Type, BachelorType::K0s)) {
              registry.fill(HIST("hMassVsPtK0s"), candV0.pT, candV0.mK0Short);
              if constexpr (DType == DType::Dstar) {
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
                if (varUtils.signD > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
                }
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                    varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value &&
                    candV0.mK0Short > cfgQaPlots.cutMassK0sMin.value &&
                    candV0.mK0Short < cfgQaPlots.cutMassK0sMax.value)) {
                  registry.fill(HIST("hMassDstarK0s"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              } else if constexpr (DType == DType::Dplus) {
                varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                    varUtils.invMassD < cfgQaPlots.cutMassDMax.value &&
                    candV0.mK0Short > cfgQaPlots.cutMassK0sMin.value &&
                    candV0.mK0Short < cfgQaPlots.cutMassK0sMax.value)) {
                  registry.fill(HIST("hMassDplusK0s"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              } // end of dType switch
            } // matched with K0s
            bool const isLambda = TESTBIT(candV0.v0Type, BachelorType::Lambda);
            bool const isAntiLambda = TESTBIT(candV0.v0Type, BachelorType::AntiLambda);
            if (isLambda || isAntiLambda) {
              registry.fill(HIST("hMassVsPtLambda"), candV0.pT, candV0.mLambda);
              if constexpr (DType == DType::Dstar) {
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
                if (varUtils.signD > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
                }
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                    varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value &&
                    candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                    candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
                  registry.fill(HIST("hMassDstarLambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              } else if constexpr (DType == DType::Dplus) {
                varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                    varUtils.invMassD < cfgQaPlots.cutMassDMax.value &&
                    candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                    candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
                  registry.fill(HIST("hMassDplusLambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              } else if constexpr (DType == DType::D0) {
                if (isLambda) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassLambda});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassLambda});
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, candV0.mom));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                      (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value)) &&
                    candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                    candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
                  if (isLambda) {
                    registry.fill(HIST("hMassD0Lambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
                  } else {
                    registry.fill(HIST("hMassD0Lambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
                  }
                }
              } // end of dType switch
            } // matched with Lambda or AntiLambda
            // fill V0 table
            // if information on V0 already stored, go to next V0
            if (!selectedV0s.count(v0.globalIndex())) {
              rowCandV0Reduced(trackPos.globalIndex(), trackNeg.globalIndex(),
                               indexHfReducedCollision,
                               candV0.pos[0], candV0.pos[1], candV0.pos[2],
                               candV0.momPos[0], candV0.momPos[1], candV0.momPos[2],
                               candV0.momNeg[0], candV0.momNeg[1], candV0.momNeg[2],
                               candV0.cosPA,
                               candV0.dcaV0ToPv,
                               nItsClsDauV0Min, nTpcCrossRowsDauV0Min, chi2TpcDauV0Max,
                               candV0.v0Type);
              selectedV0s[v0.globalIndex()] = rowCandV0Reduced.lastIndex();
            }
            fillHfCandD = true;
            // Optional filling of MC Rec table, for now only implemented for Ds1->D*K0s and Ds2*->D+K0s
            if constexpr (DoMc) {
              auto indexHfCandCharm = rowCandDmesReduced.lastIndex() + 1;
              fillMcRecoInfoDV0<DType>(particlesMc, candD, v0, tracksIU, indexHfCandCharm, selectedV0s[v0.globalIndex()], pdg, registry, rowMcRecV0Reduced);
            }
          } // end of loop on V0 candidates
        } // end of do V0s
        // Loop on the bachelor tracks
        if constexpr (DoTracks) {
          for (const auto& trackIndex : bachelorTrks) {
            auto track = tracks.rawIteratorAt(trackIndex.trackId());
            if (!isTrackSelected(track, prongIdsD, cfgSingleTrackCuts, rejectPairsWithCommonDaughter)) {
              continue;
            }
            // if the track has been reassociated, re-propagate it to PV (minor difference)
            auto trackParCovTrack = getTrackParCov(track);
            std::array<float, 2> dcaTrack{track.dcaXY(), track.dcaZ()};
            std::array<float, 3> pVecTrack = track.pVector();
            if (track.collisionId() != collision.globalIndex()) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovTrack, 2.f, matCorr, &dcaTrack);
              getPxPyPz(trackParCovTrack, pVecTrack);
            }
            registry.fill(HIST("hdEdxVsP"), track.p(), track.tpcSignal());
            // compute invariant mass and filling of QA histograms
            if constexpr (DType == DType::Dstar) {
              // D* pi
              if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
                if (varUtils.signD > 0 && track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                } else if (varUtils.signD < 0 && track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                    varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
                  registry.fill(HIST("hMassDstarPi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
              if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
                if (varUtils.signD > 0 && track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
                } else if (varUtils.signD < 0 && track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                    varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
                  registry.fill(HIST("hMassDstarK"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
              // D* p
              if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
                if (varUtils.signD > 0 && track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                } else if (varUtils.signD < 0 && track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                    varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
                  registry.fill(HIST("hMassDstarProton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
            } else if constexpr (DType == DType::Dplus) {
              // D+ pi
              if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
                if (varUtils.signD * track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                    varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
                  registry.fill(HIST("hMassDplusPi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
              // D+ K
              if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
                if (varUtils.signD * track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                    varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
                  registry.fill(HIST("hMassDplusK"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
              // D+ pr
              if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
                if (varUtils.signD * track.sign() < 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                } else {
                  varUtils.invMassReso = -1.f; // invalid case
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                    varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
                  registry.fill(HIST("hMassDplusProton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
                }
              }
            } else if constexpr (DType == DType::D0) {
              // D0 pi
              if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
                if (track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged});
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                    (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
                  if (track.sign() > 0) {
                    registry.fill(HIST("hMassD0Pi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
                  } else {
                    registry.fill(HIST("hMassD0Pi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
                  }
                }
              }
              // D0 K
              if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
                if (track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                    (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
                  if (track.sign() > 0) {
                    registry.fill(HIST("hMassD0K"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
                  } else {
                    registry.fill(HIST("hMassD0K"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
                  }
                }
              }
              // D0 p
              if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
                if (track.sign() > 0) {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
                } else {
                  varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
                }
                varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
                if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                    ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                    (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                      varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
                  if (track.sign() > 0) {
                    registry.fill(HIST("hMassD0Proton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
                  } else {
                    registry.fill(HIST("hMassD0Proton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
                  }
                }
              }
            } // end of DType switch
            // fill track table
            if (!selectedTracks.count(track.globalIndex())) {
              rowTrkReduced(track.globalIndex(),
                            indexHfReducedCollision,
                            track.px(), track.py(), track.pz(), track.sign(),
                            track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                            track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                            track.hasTOF(), track.hasTPC(), track.itsNCls(), track.tpcNClsCrossedRows(), track.tpcChi2NCl());
              selectedTracks[track.globalIndex()] = rowTrkReduced.lastIndex();
            }
            fillHfCandD = true;
            if constexpr (DoMc) {
              auto indexHfCandCharm = rowCandDmesReduced.lastIndex() + 1;
              fillMcRecoInfoDTrack<DType>(particlesMc, candD, track, tracks, indexHfCandCharm, selectedTracks[track.globalIndex()], pdg, registry, rowMcRecTrkReduced);
            }
          } // end of loop on bachelor tracks
        } // end of do tracks
        // fill D candidate table
        if (fillHfCandD) { // fill candDplus table only once per D candidate, only if at least one V0 is found
          if constexpr (DType == DType::Dplus) {
            rowCandDmesReduced(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                               indexHfReducedCollision,
                               secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                               candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                               candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                               varUtils.pVectorProng2[0], varUtils.pVectorProng2[1], varUtils.pVectorProng2[2],
                               nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax, varUtils.signD);
            if constexpr (WithMl) {
              rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
            }
          } else if constexpr (DType == DType::D0) {
            uint8_t selFlagD0 = {BIT(D0Sel::SelectedD0) | BIT(D0Sel::SelectedD0Bar)};
            if (candD.isSelD0() < cfgDmesCuts.selectionFlagD0.value) {
              CLRBIT(selFlagD0, D0Sel::SelectedD0);
            }
            if (candD.isSelD0bar() < cfgDmesCuts.selectionFlagD0Bar.value) {
              CLRBIT(selFlagD0, D0Sel::SelectedD0Bar);
            }
            rowCandDmesReduced(prongIdsD[0], prongIdsD[1],
                               indexHfReducedCollision,
                               secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                               candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                               candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                               nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax,
                               selFlagD0);
            if constexpr (WithMl) {
              rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
            }
          } else if constexpr (DType == DType::Dstar) {
            rowCandDmesReduced(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                               indexHfReducedCollision,
                               secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                               candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                               candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                               varUtils.pVectorProng2[0], varUtils.pVectorProng2[1], varUtils.pVectorProng2[2],
                               nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax,
                               nItsClsSoftPi, nTpcCrossRowsSoftPi, chi2TpcSoftPi,
                               varUtils.signD);
            if constexpr (WithMl) {
              rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
            }
          }
          fillHfReducedCollision = true;
          if constexpr (DType == DType::Dstar) {
            registry.fill(HIST("hMassVsPtDstarPaired"), candD.pt(), varUtils.invMassD - varUtils.invMassD0);
          } else if constexpr (DType == DType::Dplus) {
            registry.fill(HIST("hMassVsPtDplusPaired"), candD.pt(), varUtils.invMassD);
          } else if constexpr (DType == DType::D0) {
            if (candD.isSelD0() >= cfgDmesCuts.selectionFlagD0.value) {
              registry.fill(HIST("hMassVsPtD0Paired"), varUtils.ptD, varUtils.invMassD0);
            }
            if (candD.isSelD0bar() >= cfgDmesCuts.selectionFlagD0Bar.value) {
              registry.fill(HIST("hMassVsPtD0BarPaired"), varUtils.ptD, varUtils.invMassD0Bar);
            }
          }
        }
      } // candsD loop
      registry.fill(HIST("hEvents"), 1 + Event::Processed);
      if (!fillHfReducedCollision) {
        registry.fill(HIST("hEvents"), 1 + Event::NoDV0Selected);
        return;
      }
      registry.fill(HIST("hEvents"), 1 + Event::DV0Selected);
      // fill collision table if it contains a DPi pair a minima
      rowCollisionReduced(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
    } // end of runDataCreation function


    // Function for derived data creation
    /// \tparam dType is the D meson type (Dstar, Dplus or D0)
    /// \param mcParticles is the MC particle table
    /// \param collInfos is the reco collision table with MC info
    /// \param mcCollisions is the MC collision table
    /// \param hfEvSelMc is the HF event selection util object from MC
    /// \param rejectCollisionsWithBadEvSel is the flag to not store collisions rejected by event selection
    /// \param registry is the histogram registry
    /// \param pdg is the O2DatabasePDG service
    /// \param rowHfResoMcGenReduced is the MC gen reduced table
    template <uint8_t DType, uint8_t PairingType, typename McParticles, typename McParticlesPerMcColl, typename CCs, typename CollPerMcColl, typename McCollisions, typename TableMcGenRed, typename BCsInfo>
    void runMcGen(McParticles const& mcParticles,
                  McParticlesPerMcColl const& mcParticlesPerMcCollision,
                  CCs const& collInfos,
                  CollPerMcColl const& colPerMcCollision,
                  McCollisions const& mcCollisions,
                  o2::hf_evsel::HfEventSelectionMc& hfEvSelMc,
                  bool rejectCollisionsWithBadEvSel,
                  o2::framework::HistogramRegistry& registry,
                  o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                  TableMcGenRed& rowHfResoMcGenReduced,
                  BCsInfo const&)
    {
      bool const doV0s = (PairingType == PairingType::V0Only || PairingType == PairingType::V0AndTrack);
      bool const doTracks = (PairingType == PairingType::TrackOnly || PairingType == PairingType::V0AndTrack);
      for (const auto& mcCollision : mcCollisions) {
        // Slice the particles table to get the particles for the current MC collision
        const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
        // Slice the collisions table to get the collision info for the current MC collision
        float centrality{-1.f};
        o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
        int const nSplitColl = 0;
        const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::None>(mcCollision, collSlice, centrality);
        hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::None>(mcCollision, rejectionMask, nSplitColl);
        if (rejectCollisionsWithBadEvSel && rejectionMask != 0) {
          // at least one event selection not satisfied --> reject all gen particles from this collision
          continue;
        }
        for (const auto& particle : mcParticlesPerMcColl) {
          int8_t sign{0};
          int8_t flag{0};
          int8_t signD{0};
          int8_t signBach{0};
          int8_t origin{0};
          bool matchedReso{false}, matchedD{false}, matchedV0Tr{false};
          std::vector<int> idxBhadMothers{};
          if constexpr (DType == DType::Dstar) {
            if (doV0s) {
              // D* K0s
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDStar), +kK0}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
                  matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signBach, 2);
                  break;
                }
              }
            }
            if (doTracks && !matchedReso) {
              // D*+ pi-
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDStar), -static_cast<int>(kPiPlus)}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  matchedV0Tr = true;
                  break;
                }
              }
            }
            if (matchedReso && matchedV0Tr) {
              auto candDstarMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
              matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candDstarMC, o2::constants::physics::Pdg::kDStar, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kPiPlus)}, true, &signD, 1);
              if (matchedD) {
                auto candD0MC = mcParticles.rawIteratorAt(candDstarMC.daughtersIds().front());
                matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candD0MC, o2::constants::physics::Pdg::kD0, std::array{-kKPlus, +kPiPlus}, true, &signD, 2);
              }
            }
          } else if constexpr (DType == DType::Dplus) {
            if (doV0s) {
              // D+ K0s
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusK0s) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), +kK0}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
                  matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signBach, 2);
                  break;
                }
              }
              if (!matchedReso) {
                // D+ lambda
                for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusLambda) {
                  matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), +kLambda0}, true, &sign, 1);
                  if (matchedReso) {
                    flag = sign * decayChannelFlag;
                    auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
                    matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kLambda0, std::array{+kProton, -kPiPlus}, true, &signBach, 1);
                    break;
                  }
                }
              }
            }
            if (doTracks && !matchedReso) {
              // D+ pi-
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusPi) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), -static_cast<int>(kPiPlus)}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  matchedV0Tr = true;
                  break;
                }
              }
            }
            if (matchedReso && matchedV0Tr) {
              auto candDplusMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
              matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candDplusMC, o2::constants::physics::Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signD, 2);
            }
          } else if constexpr (DType == DType::D0) {
            if (doV0s) {
              // D0 Lambda
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Lambda) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +kLambda0}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
                  matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kLambda0, std::array{+kProton, -kPiPlus}, true, &signBach, 1);
                  break;
                }
              }
            }
            if (doTracks && !matchedReso) {
              // D0 pi+
              for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Pi) {
                matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kPiPlus)}, true, &sign, 1);
                if (matchedReso) {
                  flag = sign * decayChannelFlag;
                  matchedV0Tr = true;
                  break;
                }
              }
              // D0 K+
              if (!matchedReso) {
                for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Kplus) {
                  matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kKPlus)}, true, &sign, 1);
                  if (matchedReso) {
                    flag = sign * decayChannelFlag;
                    matchedV0Tr = true;
                    break;
                  }
                }
              }
            }
            if (matchedReso && matchedV0Tr) {
              auto candD0MC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
              matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candD0MC, o2::constants::physics::Pdg::kD0, std::array{-kKPlus, +kPiPlus}, true, &signD, 2);
            }
          }
          if (matchedReso && matchedD && matchedV0Tr) {
            origin = RecoDecay::getCharmHadronOrigin(mcParticlesPerMcColl, particle, false, &idxBhadMothers);
            registry.fill(HIST("hMCGenOrigin"), origin);
            auto ptParticle = particle.pt();
            auto invMassGen = computeInvMassGen(mcParticles, particle.globalIndex(), pdg);
            auto yParticle = RecoDecay::y(particle.pVector(), invMassGen);
            auto etaParticle = particle.eta();

            std::array<float, 2> ptProngs{};
            std::array<float, 2> yProngs{};
            std::array<float, 2> etaProngs{};
            int counter = 0;
            for (const auto& daught : particle.template daughters_as<McParticles>()) {
              ptProngs[counter] = daught.pt();
              etaProngs[counter] = daught.eta();
              yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
              counter++;
            }
            registry.fill(HIST("hMCGenCounter"), flag, ptParticle);
            rowHfResoMcGenReduced(flag, origin, ptParticle, yParticle, etaParticle,
                                  ptProngs[0], yProngs[0], etaProngs[0],
                                  ptProngs[1], yProngs[1], etaProngs[1],
                                  invMassGen, rejectionMask);
          }
        }
      }
    }
  } // namespace hf_charm_reso
} // namespace o2::analysis

#endif // PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_
