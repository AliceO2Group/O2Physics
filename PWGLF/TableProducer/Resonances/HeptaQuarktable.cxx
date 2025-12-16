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

/// \file heptaquarktable.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Junlee Kim, (junlee.kim@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/ReducedHeptaQuarkTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
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
#include <Math/Vector3D.h>
#include <TMath.h>

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct heptaquarktable {

  // Produce derived tables
  Produces<aod::RedHQEvents> redHQEvents;
  Produces<aod::HQTracks> hqTrack;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  Configurable<bool> cfgUseGlobalTrack{"cfgUseGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cfgNsigmaTPCKa{"cfgNsigmaTPCKa", 3.0, "Value of the TPC Nsigma cut for kaon"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};

  Configurable<float> cfgMinPhiMass{"cfgMinPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> cfgMaxPhiMass{"cfgMaxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<bool> cfgDeepAngleFlag{"cfgDeepAngleFlag", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 3, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 3, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};

  Configurable<float> cfgMinLambdaMass{"cfgMinLambdaMass", 1.105, "Minimum lambda mass"};
  Configurable<float> cfgMaxLambdaMass{"cfgMaxLambdaMass", 1.125, "Maximum lambda mass"};

  ConfigurableAxis massAxis{"massAxis", {200, 3.0, 3.4}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 20, 50, 100}, "Centrality interval"};
  ConfigurableAxis vertexAxis{"vertexAxis", {10, -10, 10}, "vertex axis for mixing"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>>;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;
  double massKa = o2::constants::physics::MassKPlus;

  float centrality;

  void init(o2::framework::InitContext&)
  {
    histos.add("hEventstat", "", {HistType::kTH1F, {{3, 0, 3}}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (cfgUseGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
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
    } else if (pid == 1) {
      if (std::abs(candidate.tpcNSigmaKa()) > cfgNsigmaTPCKa) {
        return false;
      }
    } else if (pid == 2) {
      if (std::abs(candidate.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr) {
        return false;
      }
    }
    return true;
  }

  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgDeepAngleFlag && angle < cfgDeepAngle) {
      return false;
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

  ROOT::Math::PxPyPzMVector DauVec1, DauVec2, HQMesonMother, HQVectorDummy, HQd1dummy, HQd2dummy;
  ROOT::Math::XYZVector HQPosVectorDummy;

  void processHQReducedTable(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    o2::aod::ITSResponse itsResponse;
    bool keepEventDoubleHQ = false;
    int numberPhi = 0;
    int numberLambda = 0;
    auto currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    centrality = collision.centFT0M();

    std::vector<int> HQId = {};

    std::vector<int64_t> HQd1Index = {};
    std::vector<int64_t> HQd2Index = {};

    std::vector<float> HQd1Charge = {};
    std::vector<float> HQd2Charge = {};

    std::vector<float> HQd1TPC = {};
    std::vector<float> HQd2TPC = {};

    std::vector<int> HQd1TOFHit = {};
    std::vector<int> HQd2TOFHit = {};

    std::vector<float> HQd1TOF = {};
    std::vector<float> HQd2TOF = {};

    std::vector<ROOT::Math::PtEtaPhiMVector> hqresonance, hqresonanced1, hqresonanced2;
    std::vector<ROOT::Math::XYZVector> hqresonancePosition;

    histos.fill(HIST("hEventstat"), 0.5);
    if (!(collision.sel8() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)))
      return;
    histos.fill(HIST("hEventstat"), 1.5);

    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1))
        continue;

      if (!selectionPID(track1, 1))
        continue;

      if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 3.0))
        continue;

      /*
              qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track1.tpcNSigmaKa(), track1.pt());
              if (track1.hasTOF()) {
                qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track1.tofNSigmaKa(), track1.pt());
              }
      */
      auto track1ID = track1.globalIndex();
      for (auto track2 : negThisColl) {
        if (!selectionTrack(track2))
          continue;

        if (!selectionPID(track2, 1))
          continue;

        if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 3.0))
          continue;

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue;

        if (!selectionPair(track1, track2))
          continue;

        DauVec1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        DauVec2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        HQMesonMother = DauVec1 + DauVec2;
        if (!(HQMesonMother.M() > cfgMinPhiMass && HQMesonMother.M() < cfgMaxPhiMass))
          continue;

        numberPhi++;
        ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massKa);
        ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(), massKa);
        ROOT::Math::PtEtaPhiMVector temp3(HQMesonMother.pt(), HQMesonMother.eta(), HQMesonMother.phi(), HQMesonMother.M());

        hqresonanced1.push_back(temp1);
        hqresonanced2.push_back(temp2);
        hqresonance.push_back(temp3);

        ROOT::Math::XYZVector temppos(0, 0, 0);
        hqresonancePosition.push_back(temppos);

        HQId.push_back(333);

        HQd1Index.push_back(track1.globalIndex());
        HQd2Index.push_back(track2.globalIndex());

        HQd1Charge.push_back(track1.sign());
        HQd2Charge.push_back(track2.sign());

        HQd1TPC.push_back(track1.tpcNSigmaKa());
        HQd2TPC.push_back(track2.tpcNSigmaKa());

        auto d1TOFHit = -1;
        auto d2TOFHit = -1;
        auto d1TOF = -999.0;
        auto d2TOF = -999.0;

        if (track1.hasTOF()) {
          d1TOFHit = 1;
          d1TOF = track1.tofNSigmaKa();
        }
        if (track2.hasTOF()) {
          d2TOFHit = 1;
          d2TOF = track2.tofNSigmaKa();
        }

        HQd1TOFHit.push_back(d1TOFHit);
        HQd2TOFHit.push_back(d2TOFHit);

        HQd1TOF.push_back(d1TOF);
        HQd2TOF.push_back(d2TOF);
      }
    }
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

        HQId.push_back(3122);
      } else if (aLambdaTag) {
        if (v0.mAntiLambda() < cfgMinLambdaMass || v0.mAntiLambda() > cfgMaxLambdaMass) {
          continue;
        }
        DauVec1 = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        DauVec2 = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        HQId.push_back(-3122);
      }
      numberLambda++;

      HQMesonMother = DauVec1 + DauVec2;

      ROOT::Math::PtEtaPhiMVector temp1(DauVec1.Pt(), DauVec1.Eta(), DauVec1.Phi(), DauVec1.M());
      ROOT::Math::PtEtaPhiMVector temp2(DauVec2.Pt(), DauVec2.Eta(), DauVec2.Phi(), DauVec2.M());
      ROOT::Math::PtEtaPhiMVector temp3(HQMesonMother.Pt(), HQMesonMother.Eta(), HQMesonMother.Phi(), HQMesonMother.M());

      hqresonanced1.push_back(temp1);
      hqresonanced2.push_back(temp2);
      hqresonance.push_back(temp3);

      ROOT::Math::XYZVector temppos(v0.x(), v0.y(), v0.z());
      hqresonancePosition.push_back(temppos);

      HQd1Index.push_back(postrack_v0.globalIndex());
      HQd2Index.push_back(negtrack_v0.globalIndex());

      HQd1Charge.push_back(postrack_v0.sign());
      HQd2Charge.push_back(negtrack_v0.sign());

      if (LambdaTag) {
        HQd1TPC.push_back(postrack_v0.tpcNSigmaPr());
        HQd2TPC.push_back(negtrack_v0.tpcNSigmaPi());
      } else if (aLambdaTag) {
        HQd1TPC.push_back(postrack_v0.tpcNSigmaPi());
        HQd2TPC.push_back(negtrack_v0.tpcNSigmaPr());
      }

      auto d1TOFHit = -1;
      auto d2TOFHit = -1;
      auto d1TOF = -999.0;
      auto d2TOF = -999.0;

      if (postrack_v0.hasTOF()) {
        d1TOFHit = 1;
        d1TOF = postrack_v0.tofNSigmaPr();
      }
      if (negtrack_v0.hasTOF()) {
        d2TOFHit = 1;
        d2TOF = negtrack_v0.tofNSigmaPr();
      } ////// TOF with PV assumption to be corrected
      HQd1TOFHit.push_back(d1TOFHit);
      HQd2TOFHit.push_back(d2TOFHit);

      HQd1TOF.push_back(d1TOF);
      HQd2TOF.push_back(d2TOF);
    } // select collision
    if (numberPhi < 2 || numberLambda < 1)
      return;

    keepEventDoubleHQ = true;

    if (keepEventDoubleHQ && numberPhi > 1 && numberLambda > 0 && (hqresonance.size() == hqresonanced1.size()) && (hqresonance.size() == hqresonanced2.size())) {
      histos.fill(HIST("hEventstat"), 2.5);
      /////////// Fill collision table///////////////
      redHQEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), centrality, numberPhi, numberLambda);
      auto indexEvent = redHQEvents.lastIndex();
      //// Fill track table for HQ//////////////////
      for (auto if1 = hqresonance.begin(); if1 != hqresonance.end(); ++if1) {
        auto i5 = std::distance(hqresonance.begin(), if1);
        HQPosVectorDummy = hqresonancePosition.at(i5);
        HQVectorDummy = hqresonance.at(i5);
        HQd1dummy = hqresonanced1.at(i5);
        HQd2dummy = hqresonanced2.at(i5);
        hqTrack(indexEvent, HQId.at(i5), HQVectorDummy.Px(), HQVectorDummy.Py(), HQVectorDummy.Pz(),
                HQd1dummy.Px(), HQd1dummy.Py(), HQd1dummy.Pz(), HQd2dummy.Px(), HQd2dummy.Py(), HQd2dummy.Pz(),
                HQPosVectorDummy.X(), HQPosVectorDummy.Y(), HQPosVectorDummy.Z(), HQVectorDummy.M(),
                HQd1Index.at(i5), HQd2Index.at(i5),
                HQd1Charge.at(i5), HQd2Charge.at(i5), HQd1TPC.at(i5), HQd2TPC.at(i5),
                HQd1TOFHit.at(i5), HQd2TOFHit.at(i5), HQd1TOF.at(i5), HQd2TOF.at(i5));
      }
    }
  } // process
  PROCESS_SWITCH(heptaquarktable, processHQReducedTable, "Process table creation for double hq", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<heptaquarktable>(cfg)};
}
