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
#include <iterator>
#include <set> // <<< CHANGED: for dedup sets
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map> // <<< CHANGED: for seenMap
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
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "N sigma TPC cut for bachelor pions"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "N sigma TOF cut for bachelor pions"};
  Configurable<bool> usePID{"usePID", false, "Flag for using PID selection"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {8, 0.0, 80}, "Mixing bins - centrality"};

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

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    auto px = candidate.pionBachPx();
    auto py = candidate.pionBachPy();
    auto pz = candidate.pionBachPz();
    auto p = std::sqrt(px * px + py * py + pz * pz);
    float lowmom = 0.5;
    if (p < lowmom) {
      if (!candidate.pionBachTOFHit() && std::abs(candidate.pionBachTPC()) < nsigmaCutTPC) {
        return true;
      } else if (candidate.pionBachTOFHit() && std::sqrt(candidate.pionBachTPC() * candidate.pionBachTPC() + candidate.pionBachTOF() * candidate.pionBachTOF()) < nsigmaCutTOF) {
        return true;
      }
    } else if (candidate.pionBachTOFHit() && std::sqrt(candidate.pionBachTPC() * candidate.pionBachTPC() + candidate.pionBachTOF() * candidate.pionBachTOF()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

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
        auto [pionPt, pionEta, pionPhi] = computePtEtaPhi(piontrack.pionBachPx(), piontrack.pionBachPy(), piontrack.pionBachPz());
        pion = ROOT::Math::PtEtaPhiMVector(pionPt, pionEta, pionPhi, o2::constants::physics::MassPionCharged);

        if (piontrack.pionBachIndex() == v0.pionIndex1() || piontrack.pionBachIndex() == v0.pionIndex2())
          continue; // checking if bachelor pion is khort daughter or not -> skip further processing if such is the case

        if (usePID && !selectionPID(piontrack))
          continue; // checking PID

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

        auto [pionPtmix, pionEtamix, pionPhimix] = computePtEtaPhi(t2.pionBachPx(), t2.pionBachPy(), t2.pionBachPz());
        if (usePID && !selectionPID(t2))
          continue; // checking PID
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
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cksspinalignder>(cfgc)};
}
