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
/// \brief Femtodream Tutorial 5
/// \author Luca Barioglio, Anton Riedel

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGCF/FemtoDream/FemtoDreamMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

// STEP 2
// Example task illustrating how to mix elements of different partitions and different events + process switches

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions,
                               aod::EvSels,
                               aod::Mults>;
using MyTracks = soa::Join<aod::FullTracks,
                           aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                           aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe>;
using MyCollision = MyCollisions::iterator;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

struct CFTutorialTask5 {
  SliceCache cache;
  Preslice<o2::aod::MyTracks> perCol = o2::aod::track::collisionId;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Defining configurables
  Configurable<float> ConfZvtxCut{"ConfZvtxCut", 10, "Z vtx cut"};
  Configurable<float> ConfEtaCut{"ConfEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> ConfMaxPtCut{"ConfMaxPtCut", 3.0, "Max Pt cut"};
  Configurable<float> ConfMinPtCut{"ConfMinPtCut", 0.5, "Min Pt cut"};
  Configurable<float> ConfMinNSigmaTPCCut{"ConfMinNSigmaTPCCut", 3., "N-sigma TPC cut"};

  // Defining filters
  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZvtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < ConfEtaCut) && (aod::track::pt > ConfMinPtCut) && (aod::track::pt < ConfMaxPtCut);

  // Applying filters
  using MyFilteredCollisions = soa::Filtered<o2::aod::MyCollisions>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredTracks = soa::Filtered<o2::aod::MyTracks>;

  Partition<MyFilteredTracks> positive = aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negative = aod::track::signed1Pt < 0.f;

  // TODO
  // define axis for event mixing, i.e. mix events with similar multiplicity and z-vertex position

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx", ";Z (cm)", kTH1F, {vtxZAxis});

    histos.add("hChargePos", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hPPos", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEtaPos", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPtPos", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPCPos", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{proton}", kTH2F, {{35, 0.5, 4.}, {100, -5., 5.}});
    histos.add("hSEPos", ";#k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});
    histos.add("hMEPos", ";#k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});

    histos.add("hChargeNeg", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hPNeg", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEtaNeg", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPtNeg", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPCNeg", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{antiproton}", kTH2F, {{35, 0.5, 4.}, {100, -5., 5.}});
    histos.add("hSENeg", ";#k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});
    histos.add("hMENeg", ";#k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});
  }

  // TODO
  // rename function to processSame event and activiate this function with a process switch (see TOD) below)
  // void processSame(MyFilteredCollision const& coll, MyFilteredTracks const& tracks)
  void process(MyFilteredCollision const& coll, MyFilteredTracks const& tracks)
  {
    auto groupPositive = positive->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    auto groupNegative = negative->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    histos.fill(HIST("hZvtx"), coll.posZ());

    for (auto track : groupPositive) {
      histos.fill(HIST("hChargePos"), track.sign());
      histos.fill(HIST("hPPos"), track.p());
      histos.fill(HIST("hPtPos"), track.pt());
      histos.fill(HIST("hEtaPos"), track.eta());
      histos.fill(HIST("hNsigmaTPCPos"), track.tpcInnerParam(), track.tpcNSigmaPr());
    }

    for (auto track : groupNegative) {
      histos.fill(HIST("hChargeNeg"), track.sign());
      histos.fill(HIST("hPNeg"), track.p());
      histos.fill(HIST("hPtNeg"), track.pt());
      histos.fill(HIST("hEtaNeg"), track.eta());
      histos.fill(HIST("hNsigmaTPCNeg"), track.tpcInnerParam(), track.tpcNSigmaPr());
    }

    float kstar = 0.;
    float mp = constants::physics::MassProton;

    for (auto& [p0, p1] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupPositive, groupPositive))) {
      if (fabs(p0.tpcNSigmaPr()) > 3. || fabs(p1.tpcNSigmaPr() > 3.)) {
        continue;
      }
      kstar = FemtoDreamMath::getkstar(p0, mp, p1, mp);
      histos.fill(HIST("hSEPos"), kstar);
    }

    for (auto& [p0, p1] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupNegative, groupNegative))) {
      if (fabs(p0.tpcNSigmaPr()) > 3. || fabs(p1.tpcNSigmaPr() > 3.)) {
        continue;
      }
      kstar = FemtoDreamMath::getkstar(p0, mp, p1, mp);
      histos.fill(HIST("hSENeg"), kstar);
    }
  }
  // TODO
  // implement process switch

  // void processMixed(MyFilteredCollisions const& colls, MyFilteredTracks const& tracks)
  // {
  //
  // TODO
  // implement process mixing using the configurable axis

  // reuse code from processSame for mixing tracks
  // auto groupPositive = positive->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
  // auto groupNegative = negative->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
  //
  // float kstar = 0.;
  // float mp = constants::physics::MassProton;
  //
  // for (auto& [p0, p1] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupPositive, groupPositive))) {
  //   if (fabs(p0.tpcNSigmaPr()) > 3. || fabs(p1.tpcNSigmaPr() > 3.)) {
  //     continue;
  //   }
  //   kstar = FemtoDreamMath::getkstar(p0, mp, p1, mp);
  //   histos.fill(HIST("hMEPos"), kstar);
  // }
  //
  // for (auto& [p0, p1] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupNegative, groupNegative))) {
  //   if (fabs(p0.tpcNSigmaPr()) > 3. || fabs(p1.tpcNSigmaPr() > 3.)) {
  //     continue;
  //   }
  //   kstar = FemtoDreamMath::getkstar(p0, mp, p1, mp);
  //   histos.fill(HIST("hMENeg"), kstar);
  // }

  // TODO
  // implement process switch
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask5>(cfgc)};
  return workflow;
}
