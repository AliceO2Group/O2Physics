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
///
/// \author Junlee Kim, (junlee.kim@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/ReducedLambdaLambdaTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TMath.h>

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lambdalambdatable {

  // Produce derived tables
  Produces<aod::RedLLEvents> redLLEvents;
  Produces<aod::LLTracks> llTrack;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  Configurable<bool> cfgUseGlobalTrack{"cfgUseGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.2f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.2f, "DCAz range for tracks"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 50, "Number of TPC cluster"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 50, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};

  Configurable<float> cfgMinLambdaMass{"cfgMinLambdaMass", 1.105, "Minimum lambda mass"};
  Configurable<float> cfgMaxLambdaMass{"cfgMaxLambdaMass", 1.125, "Maximum lambda mass"};

  ConfigurableAxis massAxis{"massAxis", {200, 2.1, 2.3}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 20, 50, 100}, "Centrality interval"};
  ConfigurableAxis vertexAxis{"vertexAxis", {10, -10, 10}, "vertex axis for mixing"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr>>;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  SliceCache cache;

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  float centrality;

  void init(o2::framework::InitContext&)
  {
    histos.add("hEventstat", "", {HistType::kTH1F, {{3, 0, 3}}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (cfgUseGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int pid)
  {
    if (pid == 0) {
      if (std::abs(candidate.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi) {
        return false;
      }
    } else if (pid == 2) {
      if (std::abs(candidate.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr) {
        return false;
      }
    }
    return true;
  }

  template <typename TCollision, typename V0>
  bool selectionV0(TCollision const& collision, V0 const& candidate)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (std::abs(candidate.dcapostopv()) < cfgDCAPosToPVMin)
      return false;
    if (std::abs(candidate.dcanegtopv()) < cfgDCANegToPVMin)
      return false;
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.yLambda() < cfgV0EtaMin)
      return false;
    if (candidate.yLambda() > cfgV0EtaMax)
      return false;
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime)
      return false;

    return true;
  }

  template <typename T>
  bool selectionV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  ROOT::Math::PxPyPzMVector DauVec1, DauVec2, LLMesonMother, LLVectorDummy, LLd1dummy, LLd2dummy;

  void processLLReducedTable(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    bool keepEventLL = false;
    int numberLambda = 0;
    auto currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    centrality = collision.centFT0M();

    std::vector<int> LLdId = {};

    std::vector<int64_t> LLdd1Index = {};
    std::vector<int64_t> LLdd2Index = {};

    std::vector<float> LLdd1TPC = {};
    std::vector<float> LLdd2TPC = {};

    std::vector<float> LLdx = {};
    std::vector<float> LLdy = {};
    std::vector<float> LLdz = {};

    std::vector<ROOT::Math::PtEtaPhiMVector> llresonance;

    histos.fill(HIST("hEventstat"), 0.5);
    if (!(collision.sel8() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)))
      return;
    histos.fill(HIST("hEventstat"), 1.5);

    for (auto& v0 : V0s) {
      auto postrack_v0 = v0.template posTrack_as<TrackCandidates>();
      auto negtrack_v0 = v0.template negTrack_as<TrackCandidates>();

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (selectionV0Daughter(postrack_v0, 0) && selectionV0Daughter(negtrack_v0, 1))
        LambdaTag = 1;

      if (selectionV0Daughter(negtrack_v0, 0) && selectionV0Daughter(postrack_v0, 1))
        aLambdaTag = 1;

      if (LambdaTag == aLambdaTag)
        continue;

      if (!selectionV0(collision, v0))
        continue;

      if (LambdaTag) {
        if (v0.mLambda() < cfgMinLambdaMass || v0.mLambda() > cfgMaxLambdaMass) {
          continue;
        }
        DauVec1 = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        DauVec2 = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
        LLdId.push_back(3122);
      } else if (aLambdaTag) {
        if (v0.mAntiLambda() < cfgMinLambdaMass || v0.mAntiLambda() > cfgMaxLambdaMass) {
          continue;
        }
        DauVec1 = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        DauVec2 = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        LLdId.push_back(-3122);
      }
      numberLambda++;

      LLdx.push_back(v0.x());
      LLdy.push_back(v0.y());
      LLdz.push_back(v0.z());

      LLMesonMother = DauVec1 + DauVec2;

      ROOT::Math::PtEtaPhiMVector temp3(LLMesonMother.Pt(), LLMesonMother.Eta(), LLMesonMother.Phi(), LLMesonMother.M());
      llresonance.push_back(temp3);

      if (LambdaTag) {
        LLdd1TPC.push_back(postrack_v0.tpcNSigmaPr());
        LLdd2TPC.push_back(negtrack_v0.tpcNSigmaPi());
      } else if (aLambdaTag) {
        LLdd1TPC.push_back(postrack_v0.tpcNSigmaPi());
        LLdd2TPC.push_back(negtrack_v0.tpcNSigmaPr());
      }

      LLdd1Index.push_back(postrack_v0.globalIndex());
      LLdd2Index.push_back(negtrack_v0.globalIndex());
    } // select collision

    if (numberLambda < 2)
      return;

    keepEventLL = true;

    if (keepEventLL) {
      histos.fill(HIST("hEventstat"), 2.5);
      /////////// Fill collision table///////////////
      redLLEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), centrality, numberLambda);
      auto indexEvent = redLLEvents.lastIndex();
      //// Fill track table for LL//////////////////
      for (auto if1 = llresonance.begin(); if1 != llresonance.end(); ++if1) {
        auto i5 = std::distance(llresonance.begin(), if1);
        LLVectorDummy = llresonance.at(i5);
        llTrack(indexEvent, LLdId.at(i5), LLVectorDummy.Px(), LLVectorDummy.Py(), LLVectorDummy.Pz(), LLdx.at(i5), LLdy.at(i5), LLdz.at(i5), LLVectorDummy.M(), LLdd1TPC.at(i5), LLdd2TPC.at(i5), LLdd1Index.at(i5), LLdd2Index.at(i5));
      }
    }
  } // process
  PROCESS_SWITCH(lambdalambdatable, processLLReducedTable, "Process table creation for double ll", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<lambdalambdatable>(cfg)};
}
