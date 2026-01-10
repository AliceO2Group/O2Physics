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

/// \file cksspinalignder.cxx
/// \brief Analysis task for ChargedKStar spin alignment
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/LFCKSSpinalignmentTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/BinningPolicy.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TMath.h>

#include <fairlogger/Logger.h>

#include <cmath> // for std::fabs
#include <deque>
// #include <iostream>
#include <algorithm>
#include <iterator>
#include <random>
#include <set> // <<< CHANGED: for dedup sets
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map> // <<< CHANGED: for seenMap
#include <unordered_set>
#include <utility>
#include <vector>

// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct cksspinalignder {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // event sel/////////
  Configurable<float> centMin{"centMin", 0, "Minimum Centrality"};
  Configurable<float> centMax{"centMax", 80, "Maximum Centrality"};
  // V0 selection ////////////
  Configurable<float> cosPA{"cosPA", 0.995, "Cosine Pointing Angle"};
  Configurable<float> radiusMin{"radiusMin", 1.2, "Minimum V0 radius"};
  Configurable<float> radiusMax{"radiusMax", 100, "Maximum V0 radius"};
  Configurable<float> dcaPion{"dcaPion", 0.1, "DCA Pion"};
  Configurable<float> dcaDaughters{"dcaDaughters", 1.0, "DCA between daughters"};
  Configurable<float> lifetimeMax{"lifetimeMax", 20, "Maximum V0 lifetime"};
  Configurable<float> ptMin{"ptMin", 0.5, "V0 Pt minimum"};
  Configurable<float> ptMax{"ptMax", 10.0, "V0 Pt maximum"};
  Configurable<float> v0eta{"v0eta", 0.8, "Eta cut on lambda"};
  // pion sel/////////
  Configurable<bool> usePID{"usePID", false, "Enable PID selection using stored bitmask"};
  // Choose which PID bit to require (1..4). If pidExclusive=true -> require ONLY that PID.
  Configurable<int> pidChoice{"pidChoice", 2, "PID choice: 1,2,3,4"};
  Configurable<bool> pidExclusive{"pidExclusive", false, "If true require ONLY the chosen PID (mask == bit). If false require (mask & bit)!=0"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {8, 0.0, 80}, "Mixing bins - centrality"};
  Configurable<float> ptMix{"ptMix", 0.2f, "ME: pT bin width"};
  Configurable<float> etaMix{"etaMix", 0.2f, "ME: eta bin width"};
  Configurable<float> phiMix{"phiMix", 0.3f, "ME: phi bin width"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {50, 1.09, 1.14}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.0}, "#it{p}_{T}"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {20, -1.0, 1.0}, "cos#it{#theta *}"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0.0, 80.0}, "Centrality"};

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{80, 0, 80}});
    histos.add("hKShortMass", "hKShortMass", kTH1F, {{100, 0.45, 0.55}});

    histos.add("hSparsesame", "hSparsesame", kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisSA, configThnAxisCentrality});
    histos.add("hSparsemix", "hSparsemix", kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisSA, configThnAxisCentrality});

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  template <typename T>
  bool selectionV0(T const& candidate)
  {
    auto kshortPt = std::sqrt(candidate.kShortPx() * candidate.kShortPx() + candidate.kShortPy() * candidate.kShortPy());
    auto px = candidate.kShortPx();
    auto py = candidate.kShortPy();
    auto pz = candidate.kShortPz();
    auto p = std::sqrt(px * px + py * py + pz * pz);
    auto kshortEta = 0.5 * std::log((p + pz) / (p - pz));

    if (std::abs(kshortEta) > v0eta) {
      return false;
    }
    if (candidate.v0Cospa() < cosPA) {
      return false;
    }
    if (candidate.v0Radius() > radiusMax) {
      return false;
    }
    if (candidate.v0Radius() < radiusMin) {
      return false;
    }
    if (candidate.dcaBetweenDaughter() > dcaDaughters) {
      return false;
    }
    if (std::abs(candidate.dcaPositive()) < dcaPion && std::abs(candidate.dcaNegative()) < dcaPion) {
      return false;
    }
    if (candidate.v0Lifetime() > lifetimeMax) {
      return false;
    }
    if (kshortPt < ptMin) {
      return false;
    }
    if (kshortPt > ptMax) {
      return false;
    }
    return true;
  }

  // ----- PID via stored bitmask (no selectionPID function needed) -----
  // Assumes your PionTracks table has a uint8_t column named pionPidMask().
  // Bit mapping: PID1=1, PID2=2, PID3=4, PID4=8.
  static constexpr uint8_t kPID1 = 1u << 0;
  static constexpr uint8_t kPID2 = 1u << 1;
  static constexpr uint8_t kPID3 = 1u << 2;
  static constexpr uint8_t kPID4 = 1u << 3;

  uint8_t requiredPidBit() const
  {
    const int choice = pidChoice.value;
    switch (choice) {
      case 1:
        return kPID1;
      case 2:
        return kPID2;
      case 3:
        return kPID3;
      case 4:
        return kPID4;
      default:
        LOGF(warn, "pidChoice=%d is invalid. Using PID2.", choice);
        return kPID2;
    }
  }

  template <typename TPion>
  bool passSelectedPID(const TPion& pion) const
  {
    if (!usePID.value) {
      return true;
    }
    const uint8_t bit = requiredPidBit();
    const uint8_t mask = pion.pionPidMask();

    if (pidExclusive.value) {
      return mask == bit; // ONLY chosen PID
    }
    return (mask & bit) != 0; // chosen PID (maybe others too)
  }
  // ---------------------------------------------------------------

  std::tuple<float, float, float> computePtEtaPhi(float px, float py, float pz)
  {
    float pt = std::sqrt(px * px + py * py);
    float p = std::sqrt(px * px + py * py + pz * pz);
    float eta = (p != std::abs(pz)) ? 0.5 * std::log((p + pz) / (p - pz)) : 0.0f; // avoid division by zero
    float phi = RecoDecay::constrainAngle(std::atan2(py, px));
    return {pt, eta, phi};
  }

  ROOT::Math::PtEtaPhiMVector kshort, pion, chkstar;
  ROOT::Math::PtEtaPhiMVector kshortmix, pionmix, chkstarmix;
  ROOT::Math::PxPyPzMVector fourVecDauCM, fourVecDauCMmix;
  ROOT::Math::XYZVector threeVecDauCM, eventplaneVecNorm, threeVecDauCMmix, eventplaneVecNormmix;

  Filter centralityFilter = (nabs(aod::kshortpionevent::cent) < centMax && nabs(aod::kshortpionevent::cent) > centMin);

  using EventCandidates = soa::Filtered<aod::KShortpionEvents>;

  void processSameData(EventCandidates::iterator const& collision, aod::KShortTracks const& V0s, aod::PionTracks const& piontracks)
  {
    auto centrality = collision.cent();
    histos.fill(HIST("hCentrality"), centrality);

    auto psiFT0C = collision.psiFT0C();

    for (const auto& v0 : V0s) {
      if (!selectionV0(v0)) {
        continue;
      }
      auto [kshortPt, kshortEta, kshortPhi] = computePtEtaPhi(v0.kShortPx(), v0.kShortPy(), v0.kShortPz());
      kshort = ROOT::Math::PtEtaPhiMVector(kshortPt, kshortEta, kshortPhi, v0.kShortMass());
      histos.fill(HIST("hKShortMass"), kshort.M());

      for (const auto& piontrack : piontracks) {

        // PID selection via stored mask
        if (!passSelectedPID(piontrack)) {
          continue;
        }
        auto [pionPt, pionEta, pionPhi] = computePtEtaPhi(piontrack.pionBachPx(), piontrack.pionBachPy(), piontrack.pionBachPz());
        pion = ROOT::Math::PtEtaPhiMVector(pionPt, pionEta, pionPhi, o2::constants::physics::MassPionCharged);

        if (piontrack.pionBachIndex() == v0.pionIndex1() || piontrack.pionBachIndex() == v0.pionIndex2())
          continue; // checking if bachelor pion is khort daughter or not -> skip further processing if such is the case

        chkstar = kshort + pion;

        ROOT::Math::Boost boost{chkstar.BoostToCM()};
        fourVecDauCM = boost(kshort);
        threeVecDauCM = fourVecDauCM.Vect();
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
        auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());

        histos.fill(HIST("hSparsesame"), chkstar.M(), chkstar.Pt(), cosThetaStar, centrality);
      }
    }
  }
  PROCESS_SWITCH(cksspinalignder, processSameData, "Process same data", true);

  // Processing Event Mixing
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::kshortpionevent::Posz, aod::kshortpionevent::Cent>;
  BinningType colBinning{{cfgVtxBins, cfgMultBins}, true};
  Preslice<aod::KShortTracks> tracksPerCollisionV0 = aod::kshortpionpair::kshortpioneventId;
  Preslice<aod::PionTracks> tracksPerCollisionBach = aod::kshortpionpair::kshortpioneventId;

  void processMixedData(EventCandidates const& collisions, aod::KShortTracks const& V0s, aod::PionTracks const& piontracks)
  {
    for (const auto& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index()) {
        continue;
      }
      auto centrality = collision1.cent();
      auto psiFT0Cmix = collision1.psiFT0C();

      auto groupV0 = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      auto groupPion = piontracks.sliceBy(tracksPerCollisionBach, collision2.index());
      for (const auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupV0, groupPion))) {
        if (!selectionV0(t1))
          continue;
        auto [kshortPtmix, kshortEtamix, kshortPhimix] = computePtEtaPhi(t1.kShortPx(), t1.kShortPy(), t1.kShortPz());
        kshortmix = ROOT::Math::PtEtaPhiMVector(kshortPtmix, kshortEtamix, kshortPhimix, t1.kShortMass());

        // PID selection via stored mask (on the bachelor pion from the mixed event)
        if (!passSelectedPID(t2)) {
          continue;
        }
        auto [pionPtmix, pionEtamix, pionPhimix] = computePtEtaPhi(t2.pionBachPx(), t2.pionBachPy(), t2.pionBachPz());
        pionmix = ROOT::Math::PtEtaPhiMVector(pionPtmix, pionEtamix, pionPhimix, o2::constants::physics::MassPionCharged);

        chkstarmix = kshortmix + pionmix;

        ROOT::Math::Boost boost{chkstarmix.BoostToCM()};
        fourVecDauCMmix = boost(kshortmix);
        threeVecDauCMmix = fourVecDauCMmix.Vect();
        eventplaneVecNormmix = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0Cmix), -std::cos(2.0 * psiFT0Cmix), 0);
        auto cosThetaStarmix = eventplaneVecNormmix.Dot(threeVecDauCMmix) / std::sqrt(threeVecDauCMmix.Mag2()) / std::sqrt(eventplaneVecNormmix.Mag2());

        histos.fill(HIST("hSparsemix"), chkstarmix.M(), chkstarmix.Pt(), cosThetaStarmix, centrality);
      }
    }
  }
  PROCESS_SWITCH(cksspinalignder, processMixedData, "Process mixed data", true);

  // =====================================================================================
  // MEV2-style mixing for Charged K* : KEEP K0s (from event E) FIXED, MIX bachelor pion
  // pion is taken from other events E' but matched in (event class + pion pt/eta/phi)
  // =====================================================================================

  // One status bin (unused dimension, kept for compatibility with linearKey)
  static constexpr int N_STATUS = 1;

  // ---------- Binner for bachelor-pion kinematics ----------
  struct MixBinnerPi {
    float ptMin, ptMax, ptStep;
    float etaMin, etaMax, etaStep;
    float phiMin, phiMax, phiStep;

    // Dummy "mass" axis (we keep 1 bin to keep the same 6D infrastructure)
    static constexpr float mMin = 0.f;
    static constexpr float mMax = 1.f;
    static constexpr int nM_ = 1;
    static constexpr float mStep = (mMax - mMin) / nM_;

    int nPt_{1}, nEta_{1}, nPhi_{1};

    MixBinnerPi(float ptMin_, float ptMax_, float ptStep_,
                float etaAbsMax_, float etaStep_,
                float phiStep_)
      : ptMin(ptMin_), ptMax(ptMax_), ptStep(ptStep_), etaMin(-etaAbsMax_), etaMax(+etaAbsMax_), etaStep(etaStep_), phiMin(0.f), phiMax(static_cast<float>(2.0 * TMath::Pi())), phiStep(phiStep_)
    {
      ptStep = (ptStep > 0.f ? ptStep : 0.2f);
      etaStep = (etaStep > 0.f ? etaStep : 0.2f);
      phiStep = (phiStep > 0.f ? phiStep : 0.2f);

      nPt_ = std::max(1, static_cast<int>(std::floor((ptMax - ptMin) / ptStep + 0.5f)));
      nEta_ = std::max(1, static_cast<int>(std::floor((etaMax - etaMin) / etaStep + 0.5f)));
      nPhi_ = std::max(1, static_cast<int>(std::ceil((phiMax - phiMin) / phiStep)));
    }

    inline int nPt() const { return nPt_; }
    inline int nEta() const { return nEta_; }
    inline int nPhi() const { return nPhi_; }
    inline int nM() const { return nM_; }

    inline int binFromValue(float v, float vmin, float step, int nBins) const
    {
      if (!std::isfinite(v))
        return -1;
      const float x = (v - vmin) / step;
      int b = static_cast<int>(std::floor(x + 1e-6f));
      if (b < 0)
        return -1;
      if (b >= nBins)
        b = nBins - 1;
      return b;
    }

    inline int ptBin(float pt) const { return binFromValue(pt, ptMin, ptStep, nPt_); }
    inline int etaBin(float eta) const { return binFromValue(eta, etaMin, etaStep, nEta_); }
    inline int phiBin(float phi) const { return binFromValue(phi, phiMin, phiStep, nPhi_); }
    inline int massBin(float m) const { return binFromValue(m, mMin, mStep, nM_); } // always 0 if finite
  };

  // Buffer candidate: one pion stored in a bin
  struct BufferCandPi {
    int64_t collisionIdx; // to enforce different event
    int64_t rowIndex;     // row id in PionTracks (for iteratorAt)
    uint16_t ptBin, etaBin, phiBin, mBin;
  };

  struct MatchRef {
    int64_t collisionIdx;
    int64_t rowIndex;
  };

  // Key: (colBin, stat, pt, eta, phi, m)
  static inline size_t linearKey(int colBin, int statBin,
                                 int ptBin, int etaBin, int phiBin, int mBin,
                                 int nStat, int nPt, int nEta, int nPhi, int nM)
  {
    return ((((((static_cast<size_t>(colBin) * nStat + statBin) * nPt + ptBin) * nEta + etaBin) * nPhi + phiBin) * nM + mBin));
  }

  // =====================================================================================
  // Mixed-event method:
  // PASS 1: build buffer of pions by (event-class bin) + (pi pt/eta/phi bins)
  // PASS 2: loop same-event (K0, piSame) pairs; replace piSame by piX from buffer
  //         keep K0 from current event fixed; use psiFT0C of current event.
  // =====================================================================================
  void processMixedDataMEV2_Pion(EventCandidates const& collisions,
                                 aod::KShortTracks const& V0s,
                                 aod::PionTracks const& piontracks)
  {
    static thread_local std::mt19937 rng(0xBADC0FFE);

    // Build pion binner from your configurables
    // Use your own pt range config for bachelor pions if you have it; otherwise reuse ptMin/ptMax
    MixBinnerPi mb{
      ptMin, ptMax, /*ptStep*/ ptMix.value,
      /*|eta|max*/ v0eta.value, /*etaStep*/ etaMix.value,
      /*phiStep*/ phiMix.value};

    const int nCol = colBinning.getAllBinsCount(); // event-class bins (vz, cent)
    const int nStat = N_STATUS;                    // 1
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCandPi>> buffer(nKeys);

    // ---------------- PASS 1: fill 6D buffer with bachelor pions ----------------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      if (colBin < 0)
        continue;

      auto slicePi = piontracks.sliceBy(tracksPerCollisionBach, col.index());

      for (auto const& pi : slicePi) {

        // apply configurable PID selection stored in bitmask (your passSelectedPID)
        if (!passSelectedPID(pi))
          continue;

        auto [pt, eta, phi] = computePtEtaPhi(pi.pionBachPx(), pi.pionBachPy(), pi.pionBachPz());
        phi = RecoDecay::constrainAngle(phi, 0.0F);

        const int ptB = mb.ptBin(pt);
        const int etaB = mb.etaBin(eta);
        const int phiB = mb.phiBin(phi);
        const int mB = 0; // dummy axis has 1 bin
        if (ptB < 0 || etaB < 0 || phiB < 0)
          continue;

        const int stat = 0;
        const size_t key = linearKey(colBin, stat, ptB, etaB, phiB, mB, nStat, nPt, nEta, nPhi, nM);

        const int64_t row = static_cast<int64_t>(pi.globalIndex());

        buffer[key].push_back(BufferCandPi{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = row,
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
      }
    }

    // ---------------- PASS 2: same-event (K0,piSame), replace ONLY pion by piX ----------------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col1.posz(), col1.cent()));
      if (colBin < 0)
        continue;

      auto poolK0 = V0s.sliceBy(tracksPerCollisionV0, col1.index());
      auto poolPi = piontracks.sliceBy(tracksPerCollisionBach, col1.index());
      if (poolK0.size() == 0 || poolPi.size() == 0)
        continue;

      const float centrality = col1.cent();
      const float psiFT0Cmix = col1.psiFT0C();

      // Loop over K0s from the *current* event (fixed object)
      for (auto const& k0 : poolK0) {
        if (!selectionV0(k0))
          continue;

        // Build K0 4-vector once
        auto [kpt, keta, kphi] = computePtEtaPhi(k0.kShortPx(), k0.kShortPy(), k0.kShortPz());
        kphi = RecoDecay::constrainAngle(kphi, 0.0F);
        kshortmix = ROOT::Math::PtEtaPhiMVector(kpt, keta, kphi, k0.kShortMass());

        // Loop over same-event pions to define the "same-event pair" base
        for (auto const& piSame : poolPi) {
          if (!passSelectedPID(piSame))
            continue;

          // avoid self-combination: bachelor pion can't be a K0 daughter
          if (piSame.pionBachIndex() == k0.pionIndex1() || piSame.pionBachIndex() == k0.pionIndex2()) {
            continue;
          }

          // Compute the kinematic bin of THIS pion (this defines where to pick mixed pions)
          auto [ppt, peta, pphi] = computePtEtaPhi(piSame.pionBachPx(), piSame.pionBachPy(), piSame.pionBachPz());
          pphi = RecoDecay::constrainAngle(pphi, 0.0F);

          const int ptB = mb.ptBin(ppt);
          const int etaB = mb.etaBin(peta);
          const int phiB = mb.phiBin(pphi);
          const int mB = 0;
          if (ptB < 0 || etaB < 0 || phiB < 0)
            continue;

          const int stat = 0;
          const size_t key = linearKey(colBin, stat, ptB, etaB, phiB, mB, nStat, nPt, nEta, nPhi, nM);
          auto const& binVec = buffer[key];
          if (binVec.empty())
            continue;

          // collect partners from same bin but different collision
          std::vector<MatchRef> matches;
          matches.reserve(binVec.size());
          const int64_t curColIdx = static_cast<int64_t>(col1.index());

          for (auto const& bc : binVec) {
            if (bc.collisionIdx == curColIdx)
              continue;
            matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
          }
          if (matches.empty())
            continue;

          // random unique sampling (cap = nEvtMixing)
          const int cap = nEvtMixing.value;
          if (cap > 0 && cap < static_cast<int>(matches.size())) {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(matches.size()) - 1);
            std::unordered_set<int> chosen;
            chosen.reserve(static_cast<size_t>(cap) * 2);
            while (static_cast<int>(chosen.size()) < cap)
              chosen.insert(dist(rng));

            std::vector<MatchRef> subset;
            subset.reserve(cap);
            for (int idx : chosen)
              subset.push_back(matches[idx]);
            matches.swap(subset);
          } else {
            std::shuffle(matches.begin(), matches.end(), rng);
          }

          // const float wBase = 1.0f / static_cast<float>(matches.size());
          const float wBase = 1.0;

          // Replace pion by mixed pion piX, keep k0 fixed
          for (auto const& m : matches) {
            auto piX = piontracks.iteratorAt(m.rowIndex); // must match PASS 1 rowIndex convention

            // safety: PID selection again
            if (!passSelectedPID(piX))
              continue;

            // OPTIONAL keep: also avoid piX being a K0 daughter (same k0 from current event)
            if (piX.pionBachIndex() == k0.pionIndex1() || piX.pionBachIndex() == k0.pionIndex2()) {
              continue;
            }

            auto [xpt, xeta, xphi] = computePtEtaPhi(piX.pionBachPx(), piX.pionBachPy(), piX.pionBachPz());
            xphi = RecoDecay::constrainAngle(xphi, 0.0F);
            pionmix = ROOT::Math::PtEtaPhiMVector(xpt, xeta, xphi, o2::constants::physics::MassPionCharged);

            chkstarmix = kshortmix + pionmix;

            ROOT::Math::Boost boost{chkstarmix.BoostToCM()};
            fourVecDauCMmix = boost(kshortmix);
            threeVecDauCMmix = fourVecDauCMmix.Vect();

            eventplaneVecNormmix = ROOT::Math::XYZVector(std::sin(2.f * psiFT0Cmix),
                                                         -std::cos(2.f * psiFT0Cmix), 0.f);

            const float cosThetaStarmix =
              eventplaneVecNormmix.Dot(threeVecDauCMmix) /
              std::sqrt(threeVecDauCMmix.Mag2()) /
              std::sqrt(eventplaneVecNormmix.Mag2());

            histos.fill(HIST("hSparsemix"),
                        chkstarmix.M(), chkstarmix.Pt(), cosThetaStarmix, centrality, wBase);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(cksspinalignder, processMixedDataMEV2_Pion, "Process mixed data MEV4 (keep K0 fixed, replace pion with mixed pion)", true);

  /*
// ---- Configuration knobs (already in your task or add similarly) ----
// ptMix: step for K0 pT binning, etaMix step for eta, phiMix step for phi
// configThnAxisInvMass etc are your histogram axes (not used for binning here)

// We can keep one status bin for now
static constexpr int N_STATUS = 1;

// Binner for K0 kinematics
struct MixBinnerK0 {
  float ptMin, ptMax, ptStep;
  float etaMin, etaMax, etaStep;
  float phiMin, phiMax, phiStep;

  // K0 mass binning (you can widen if you want)
  static constexpr float mMin = 0.45f;
  static constexpr float mMax = 0.55f; // exclusive-ish
  static constexpr int   nM_  = 1;
  static constexpr float mStep = (mMax - mMin) / nM_;

  int nPt_, nEta_, nPhi_;

  MixBinnerK0(float ptMin_, float ptMax_, float ptStep_,
              float etaAbsMax_, float etaStep_,
              float phiStep_)
    : ptMin(ptMin_), ptMax(ptMax_), ptStep(ptStep_),
      etaMin(-etaAbsMax_), etaMax(+etaAbsMax_), etaStep(etaStep_),
      phiMin(0.f), phiMax(static_cast<float>(2.0 * TMath::Pi())), phiStep(phiStep_)
  {
    ptStep  = (ptStep  > 0.f ? ptStep  : 0.1f);
    etaStep = (etaStep > 0.f ? etaStep : 0.1f);
    phiStep = (phiStep > 0.f ? phiStep : 0.1f);

    nPt_  = std::max(1, static_cast<int>(std::floor((ptMax - ptMin) / ptStep + 0.5f)));
    nEta_ = std::max(1, static_cast<int>(std::floor((etaMax - etaMin) / etaStep + 0.5f)));
    nPhi_ = std::max(1, static_cast<int>(std::ceil((phiMax - phiMin) / phiStep)));
  }

  inline int nPt()  const { return nPt_; }
  inline int nEta() const { return nEta_; }
  inline int nPhi() const { return nPhi_; }
  inline int nM()   const { return nM_; }

  inline int binFromValue(float v, float vmin, float step, int nBins) const
  {
    if (!std::isfinite(v)) return -1;
    const float x = (v - vmin) / step;
    int b = static_cast<int>(std::floor(x + 1e-6f));
    if (b < 0) return -1;
    if (b >= nBins) b = nBins - 1;
    return b;
  }

  inline int ptBin(float pt)   const { return binFromValue(pt,  ptMin,  ptStep,  nPt_); }
  inline int etaBin(float eta) const { return binFromValue(eta, etaMin, etaStep, nEta_); }
  inline int phiBin(float phi) const { return binFromValue(phi, phiMin, phiStep, nPhi_); }
  inline int massBin(float m)  const { return binFromValue(m,  mMin,   mStep,   nM_); }
};

  struct BufferCandK0 {
    int64_t collisionIdx; // to enforce different event
    int64_t rowIndex;     // row id in KShortTracks (for iteratorAt)
    uint16_t ptBin, etaBin, phiBin, mBin;
  };

  struct MatchRef {
    int64_t collisionIdx;
    int64_t rowIndex;
  };

  // Key: (colBin, stat, pt, eta, phi, m)
  static inline size_t linearKey(int colBin, int statBin,
         int ptBin, int etaBin, int phiBin, int mBin,
         int nStat, int nPt, int nEta, int nPhi, int nM)
  {
    return ((((((static_cast<size_t>(colBin) * nStat + statBin) * nPt + ptBin) * nEta + etaBin) * nPhi + phiBin) * nM + mBin));
  }

  void processMixedDataMEV4(EventCandidates const& collisions,
          aod::KShortTracks const& V0s,
          aod::PionTracks const& piontracks)
  {
    static thread_local std::mt19937 rng(0xBADC0FFE);

    // Build K0 binner from your configurables
    MixBinnerK0 mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0eta.value, etaMix.value,
      phiMix.value
    };

    const int nCol  = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt   = mb.nPt();
    const int nEta  = mb.nEta();
    const int nPhi  = mb.nPhi();
    const int nM    = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCandK0>> buffer(nKeys);

    // ---------------- PASS 1: fill 6D buffer with K0s ----------------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      if (colBin < 0) continue;

      auto sliceK0 = V0s.sliceBy(tracksPerCollisionV0, col.index());

      for (auto const& k0 : sliceK0) {
  if (!selectionV0(k0)) continue;

  auto [kpt, keta, kphi] = computePtEtaPhi(k0.kShortPx(), k0.kShortPy(), k0.kShortPz());
  kphi = RecoDecay::constrainAngle(kphi, 0.0F);

  const int ptB  = mb.ptBin(kpt);
  const int etaB = mb.etaBin(keta);
  const int phiB = mb.phiBin(kphi);
  const int mB   = mb.massBin(k0.kShortMass());
  if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) continue;

  const int stat = 0;
  const size_t key = linearKey(colBin, stat, ptB, etaB, phiB, mB, nStat, nPt, nEta, nPhi, nM);

  // IMPORTANT: rowIndex must match iteratorAt accessor below.
  const int64_t row = static_cast<int64_t>(k0.globalIndex()); // if this fails -> use k0.index() consistently

  buffer[key].push_back(BufferCandK0{
      .collisionIdx = static_cast<int64_t>(col.index()),
      .rowIndex     = row,
      .ptBin        = static_cast<uint16_t>(ptB),
      .etaBin       = static_cast<uint16_t>(etaB),
      .phiBin       = static_cast<uint16_t>(phiB),
      .mBin         = static_cast<uint16_t>(mB)
    });
      }
    }

    // ---------------- PASS 2: build same-event (K0,π) pairs, replace only K0 by K0X ----------------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col1.posz(), col1.cent()));
      if (colBin < 0) continue;

      auto poolK0 = V0s.sliceBy(tracksPerCollisionV0, col1.index());
      auto poolPi = piontracks.sliceBy(tracksPerCollisionBach, col1.index());

      if (poolK0.size() == 0 || poolPi.size() == 0) continue; // no .empty() in ASoA slice

      const float centrality = col1.cent();
      const float psiFT0Cmix = col1.psiFT0C();

      // same-event pair loop
      for (auto const& k0 : poolK0) {
  if (!selectionV0(k0)) continue;

  auto [kpt, keta, kphi] = computePtEtaPhi(k0.kShortPx(), k0.kShortPy(), k0.kShortPz());
  kphi = RecoDecay::constrainAngle(kphi, 0.0F);

  const int ptB  = mb.ptBin(kpt);
  const int etaB = mb.etaBin(keta);
  const int phiB = mb.phiBin(kphi);
  const int mB   = mb.massBin(k0.kShortMass());
  if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) continue;

  const int stat = 0;
  const size_t key = linearKey(colBin, stat, ptB, etaB, phiB, mB, nStat, nPt, nEta, nPhi, nM);
  auto const& binVec = buffer[key];
  if (binVec.empty()) continue;

  // pre-collect eligible K0 partners from other events (same bin, different collision)
  std::vector<MatchRef> matches;
  matches.reserve(binVec.size());
  const int64_t curColIdx = static_cast<int64_t>(col1.index());

  for (auto const& bc : binVec) {
    if (bc.collisionIdx == curColIdx) continue;
    matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
  }
  if (matches.empty()) continue;

  // choose random unique subset of K0 partners
  const int cap = nEvtMixing.value;
  if (cap > 0 && cap < static_cast<int>(matches.size())) {
    std::uniform_int_distribution<int> dist(0, static_cast<int>(matches.size()) - 1);
    std::unordered_set<int> chosen;
    chosen.reserve(static_cast<size_t>(cap) * 2);
    while (static_cast<int>(chosen.size()) < cap) chosen.insert(dist(rng));

    std::vector<MatchRef> subset;
    subset.reserve(cap);
    for (int idx : chosen) subset.push_back(matches[idx]);
    matches.swap(subset);
  } else {
    std::shuffle(matches.begin(), matches.end(), rng);
  }

  const float wBase = 1.0f / static_cast<float>(matches.size());

  // Now keep π from SAME event fixed (this is your requirement)
  for (auto const& piSame : poolPi) {
    if (!passSelectedPID(piSame)) continue;

    // avoid self-combination if you stored K0 daughter indices
    if (piSame.pionBachIndex() == k0.pionIndex1() || piSame.pionBachIndex() == k0.pionIndex2()) {
      continue;
    }

    auto [ppt, peta, pphi] = computePtEtaPhi(piSame.pionBachPx(), piSame.pionBachPy(), piSame.pionBachPz());
    pphi = RecoDecay::constrainAngle(pphi, 0.0F);
    pionmix = ROOT::Math::PtEtaPhiMVector(ppt, peta, pphi, o2::constants::physics::MassPionCharged);

    // replace K0 by mixed K0X from other event, but keep piSame
    for (auto const& m : matches) {
      auto k0X = V0s.iteratorAt(m.rowIndex); // if globalIndex fails -> use index consistently in PASS1

      if (!selectionV0(k0X)) continue;

      auto [xpt, xeta, xphi] = computePtEtaPhi(k0X.kShortPx(), k0X.kShortPy(), k0X.kShortPz());
      xphi = RecoDecay::constrainAngle(xphi, 0.0F);
      kshortmix = ROOT::Math::PtEtaPhiMVector(xpt, xeta, xphi, k0X.kShortMass());

      chkstarmix = kshortmix + pionmix;

      ROOT::Math::Boost boost{chkstarmix.BoostToCM()};
      fourVecDauCMmix = boost(kshortmix);
      threeVecDauCMmix = fourVecDauCMmix.Vect();

      eventplaneVecNormmix = ROOT::Math::XYZVector(std::sin(2.f * psiFT0Cmix),
               -std::cos(2.f * psiFT0Cmix), 0.f);

      const float cosThetaStarmix =
        eventplaneVecNormmix.Dot(threeVecDauCMmix) /
        std::sqrt(threeVecDauCMmix.Mag2()) /
        std::sqrt(eventplaneVecNormmix.Mag2());

      histos.fill(HIST("hSparsemix"),
      chkstarmix.M(), chkstarmix.Pt(), cosThetaStarmix, centrality, wBase);
    }
  }
      }
    }
  }

  PROCESS_SWITCH(cksspinalignder, processMixedDataMEV4,
     "Process mixed data MEV4 (replace K0, keep same-event pion)", true);
  */
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cksspinalignder>(cfgc)};
}
