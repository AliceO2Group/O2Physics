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

/// \author Luca Barioglio

// O2 includes
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
  Preslice<o2::aod::MyTracks> perCol = aod::track::collisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Defining configurables
  Configurable<float> ConfZvtxCut{"ConfZvtxCut", 10, "Z vtx cut"};
  Configurable<float> ConfEtaCut{"ConfEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> ConfMaxPtCut{"ConfMaxPtCut", 3.0, "Max Pt cut"};
  Configurable<float> ConfMinPtCut{"ConfMinPtCut", 0.5, "Min Pt cut"};
  Configurable<float> ConfMinNSigmaTPCCut{"ConfMinNSigmaTPCCut", 3., "N-sigma TPC cut"};
  Configurable<float> ConfChargeCut{"ConfChargeCut", 0., "N-sigma TPC cut"};

  // Defining filters
  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZvtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < ConfEtaCut) && (aod::track::pt > ConfMinPtCut) && (aod::track::pt < ConfMaxPtCut);

  // Applying filters
  using MyFilteredCollisions = soa::Filtered<o2::aod::MyCollisions>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredTracks = soa::Filtered<o2::aod::MyTracks>;

  Partition<MyFilteredTracks> positive = aod::track::signed1Pt > ConfChargeCut;
  Partition<MyFilteredTracks> negative = aod::track::signed1Pt < ConfChargeCut;

  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0A>;

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
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEta", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPCP", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{35, 0.5, 4.}, {100, -5., 5.}});
    histos.add("hChargePos", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hChargeNeg", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hInvariantMass", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});
    histos.add("hInvariantMassMixed", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});
    histos.add("hInvariantMassMixedInterface", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});
  }

  void processSame(MyFilteredCollision const& coll, MyFilteredTracks const&)
  {
    auto groupPositive = positive->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    auto groupNegative = negative->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    histos.fill(HIST("hZvtx"), coll.posZ());

    for (auto track : groupPositive) {
      histos.fill(HIST("hChargePos"), track.sign());
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hNsigmaTPCP"), track.p(), track.tpcNSigmaPi());
    }

    for (auto track : groupNegative) {
      histos.fill(HIST("hChargeNeg"), track.sign());
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hNsigmaTPCP"), track.p(), track.tpcNSigmaPi());
    }

    for (auto& [pos, neg] : combinations(soa::CombinationsFullIndexPolicy(groupPositive, groupNegative))) {
      if (fabs(pos.tpcNSigmaPi()) > 3 or fabs(neg.tpcNSigmaPi()) > 3) {
        continue;
      }
      TLorentzVector posVec;
      posVec.SetPtEtaPhiM(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassPionCharged);
      TLorentzVector negVec;
      negVec.SetPtEtaPhiM(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassPionCharged);

      TLorentzVector sumVec(posVec);
      sumVec += negVec;
      histos.fill(HIST("hInvariantMass"), sumVec.M());
    }
  }
  PROCESS_SWITCH(CFTutorialTask5, processSame, "Enable processing same event", true);

  void processMixed(MyFilteredCollisions const& colls, MyFilteredTracks const&)
  {
    BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, colls, colls)) {
      auto groupPositive = positive->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto groupNegative = negative->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      for (auto& [pos, neg] : combinations(soa::CombinationsFullIndexPolicy(groupPositive, groupNegative))) {
        if (fabs(pos.tpcNSigmaPi()) > 3 or fabs(neg.tpcNSigmaPi()) > 3) {
          continue;
        }
        TLorentzVector posVec;
        posVec.SetPtEtaPhiM(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassPionCharged);
        TLorentzVector negVec;
        negVec.SetPtEtaPhiM(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassPionCharged);

        TLorentzVector sumVec(posVec);
        sumVec += negVec;
        histos.fill(HIST("hInvariantMassMixed"), sumVec.M());
      }
    }
  }
  PROCESS_SWITCH(CFTutorialTask5, processMixed, "Enable processing mixed event", true);

  void processMixedEventInterface(MyFilteredCollisions& colls, MyFilteredTracks& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};
    SameKindPair<MyFilteredCollisions, MyFilteredTracks, BinningType> pair{colBinning, 5, -1, colls, tracksTuple, &cache};
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      Partition<MyFilteredTracks> groupPositive = aod::track::signed1Pt > ConfChargeCut;
      groupPositive.bindTable(tracks1);
      Partition<MyFilteredTracks> groupNegative = aod::track::signed1Pt < ConfChargeCut;
      groupNegative.bindTable(tracks2);

      for (auto& [pos, neg] : combinations(soa::CombinationsFullIndexPolicy(groupPositive, groupNegative))) {
        if (fabs(pos.tpcNSigmaPi()) > 3 or fabs(neg.tpcNSigmaPi()) > 3) {
          continue;
        }
        TLorentzVector posVec;
        posVec.SetPtEtaPhiM(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassPionCharged);
        TLorentzVector negVec;
        negVec.SetPtEtaPhiM(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassPionCharged);

        TLorentzVector sumVec(posVec);
        sumVec += negVec;
        histos.fill(HIST("hInvariantMassMixedInterface"), sumVec.M());
      }
    }
  }
  PROCESS_SWITCH(CFTutorialTask5, processMixedEventInterface, "Enable processing mixed event with standard mixing interface", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask5>(cfgc)};
  return workflow;
}
