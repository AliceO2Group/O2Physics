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
/// \brief this is a starting point for the Resonances tutorial
/// \author sourav kundu
/// \since 02/11/2023

#include <Framework/Configurable.h>
#include <TLorentzVector.h>
#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TMath.h>
#include <fairlogger/Logger.h>
#include <iostream>
#include <iterator>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/ReducedF1ProtonTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct f1protoncorrelation {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // PID selection
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  // PID selection
  Configurable<bool> fillRotation{"fillRotation", 1, "Fill rotation"};
  Configurable<int> strategyPIDPion{"strategyPIDPion", 0, "PID strategy Pion"};
  Configurable<int> strategyPIDKaon{"strategyPIDKaon", 0, "PID strategy Kaon"};
  Configurable<float> maxKKS0Mass{"maxKKS0Mass", 1.025, "Maximum kaon kshort mass"};
  Configurable<float> maxMomentumPion{"maxMomentumPion", 4.0, "Maximum momentum Pion"};
  Configurable<float> maxMomentumKaon{"maxMomentumKaon", 4.0, "Maximum momentum Kaon"};
  Configurable<float> momentumTOFPion{"momentumTOFPion", 0.8, "Pion momentum TOF"};
  Configurable<float> momentumTOFKaon{"momentumTOFKaon", 0.8, "Kaon momentum TOF"};
  Configurable<float> momentumTOFProton{"momentumTOFProton", 0.7, "Proton momentum TOF"};
  Configurable<float> lowPtF1{"lowPtF1", 1.0, "PT cut F1"};
  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 1, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 30.0, 40.0, 50.0, 60.0, 80.0, 200.0}, "Mixing bins - number of contributor"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hNsigmaProtonTPCSE", "Nsigma Proton TPC distribution same event", kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 1.0f}});
    histos.add("hNsigmaProtonTPCME", "Nsigma Proton TPC distribution mixed event", kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 1.0f}});
    histos.add("h2SameEventPtCorrelation", "Pt correlation of F1 and proton", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {100, 0.0, 10.0}});

    histos.add("h2SameEventInvariantMassUnlike_mass", "Unlike Sign Invariant mass of f1 same event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});
    histos.add("h2SameEventInvariantMassLike_mass", "Like Sign Invariant mass of f1 same event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});
    histos.add("h2SameEventInvariantMassRot_mass", "Rotational Invariant mass of f1 same event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});

    histos.add("h2MixEventInvariantMassUnlike_mass", "Unlike Sign Invariant mass of f1 mix event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});
    histos.add("h2MixEventInvariantMassLike_mass", "Like Sign Invariant mass of f1 mix event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});
    histos.add("h2MixEventInvariantMassRot_mass", "Rotational Sign Invariant mass of f1 mix event", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {800, 1.0, 1.8}});
  }

  // get kstar
  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getkstar(const TLorentzVector part1,
                 const TLorentzVector part2)
  {
    // const TLorentzVector trackSum = part1 + part2;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    // TLorentzVector PartOneCMS(part1);
    // TLorentzVector PartTwoCMS(part2);
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    // const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  TLorentzVector F1, Proton, F1ProtonPair, Pion, Kaon, Kshort;
  TLorentzVector F1Rot, PionRot, KaonKshortPair;
  // Process the data in same event
  void process(aod::RedF1PEvents::iterator const& /*collision*/, aod::F1Tracks const& f1tracks, aod::ProtonTracks const& protontracks)
  {
    for (auto f1track : f1tracks) {
      if (f1track.f1MassKaonKshort() > maxKKS0Mass) {
        continue;
      }
      F1.SetXYZM(f1track.f1Px(), f1track.f1Py(), f1track.f1Pz(), f1track.f1Mass());
      Pion.SetXYZM(f1track.f1d1Px(), f1track.f1d1Py(), f1track.f1d1Pz(), 0.139);
      Kaon.SetXYZM(f1track.f1d2Px(), f1track.f1d2Py(), f1track.f1d2Pz(), 0.493);
      Kshort.SetXYZM(f1track.f1d3Px(), f1track.f1d3Py(), f1track.f1d3Pz(), 0.497);
      KaonKshortPair = Kaon + Kshort;
      if (Pion.P() > maxMomentumPion || Kaon.P() > maxMomentumKaon) {
        continue;
      }
      if (strategyPIDPion == 1 && Pion.P() > momentumTOFPion && f1track.f1d1TOFHit() != 1) {
        continue;
      }
      if (strategyPIDKaon == 1 && Kaon.P() > momentumTOFKaon && f1track.f1d2TOFHit() != 1) {
        continue;
      }
      for (auto protontrack : protontracks) {
        Proton.SetXYZM(protontrack.protonPx(), protontrack.protonPy(), protontrack.protonPz(), 0.938);
        if (Proton.P() < momentumTOFProton && TMath::Abs(protontrack.protonNsigmaTPC()) > 3) {
          continue;
        }
        if (Proton.P() >= momentumTOFProton && protontrack.protonTOFHit() != 1 && TMath::Abs(protontrack.protonNsigmaTOF()) > 3) {
          continue;
        }
        if ((f1track.f1PionIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KaonIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KshortPositiveIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KshortNegativeIndex() == protontrack.f1ProtonIndex())) {
          continue;
        }
        auto relative_momentum = getkstar(F1, Proton);
        histos.fill(HIST("h2SameEventPtCorrelation"), relative_momentum, F1.Pt(), Proton.Pt());
        if (f1track.f1SignalStat() == 1) {
          histos.fill(HIST("h2SameEventInvariantMassUnlike_mass"), relative_momentum, F1.Pt(), F1.M()); // F1 sign = 1 unlike, F1 sign = -1 like
          histos.fill(HIST("hNsigmaProtonTPCSE"), protontrack.protonNsigmaTPC(), relative_momentum);
        }
        if (f1track.f1SignalStat() == -1) {
          histos.fill(HIST("h2SameEventInvariantMassLike_mass"), relative_momentum, F1.Pt(), F1.M());
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < 9; nrotbkg++) {
            auto anglestart = 5.0 * TMath::Pi() / 6.0;
            auto angleend = 7.0 * TMath::Pi() / 6.0;
            auto anglestep = (angleend - anglestart) / (1.0 * (9.0 - 1.0));
            auto rotangle = anglestart + nrotbkg * anglestep;
            auto rotPionPx = f1track.f1d1Px() * std::cos(rotangle) - f1track.f1d1Py() * std::sin(rotangle);
            auto rotPionPy = f1track.f1d1Px() * std::sin(rotangle) + f1track.f1d1Py() * std::cos(rotangle);
            PionRot.SetXYZM(rotPionPx, rotPionPy, f1track.f1d1Pz(), 0.139);
            F1Rot = PionRot + KaonKshortPair;
            auto relative_momentum_rot = getkstar(F1Rot, Proton);
            if (f1track.f1SignalStat() == 1) {
              histos.fill(HIST("h2SameEventInvariantMassRot_mass"), relative_momentum_rot, F1Rot.Pt(), F1Rot.M());
            }
          }
        }
      }
    }
  }

  // Processing Event Mixing
  // using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  // for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
  // Pair<aod::RedF1PEvents, aod::F1Tracks, aod::ProtonTracks, BinningType> pairs{colBinning, nEvtMixing, -1, &cache}; // -1 is the number of the bin to skip
  //
  // tracks1 is an aod::Tracks table of f1tracks belonging to collision collision1 (aod::Collision::iterator)
  // tracks2 is an aod::Tracks table of protontracks belonging to collision collision2 (aod::Collision::iterator)
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  BinningType colBinning{{CfgVtxBins, CfgMultBins}, true};
  Preslice<aod::F1Tracks> tracksPerCollisionPresliceF1 = aod::f1protondaughter::redF1PEventId;
  Preslice<aod::ProtonTracks> tracksPerCollisionPresliceP = aod::f1protondaughter::redF1PEventId;
  void processME(aod::RedF1PEvents& collisions, aod::F1Tracks& f1tracks, aod::ProtonTracks& protontracks)
  {
    for (auto& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing, -1, collisions, collisions)) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", collision1.index(), collision2.index());
      if (collision1.index() == collision2.index()) {
        continue;
      }
      if (f1tracks.size() == 0 || protontracks.size() == 0) {
        continue;
      }
      auto groupF1 = f1tracks.sliceBy(tracksPerCollisionPresliceF1, collision1.globalIndex());
      auto groupProton = protontracks.sliceBy(tracksPerCollisionPresliceP, collision2.globalIndex());
      // auto groupF1 = f1tracks.sliceByCached(aod::f1protondaughter::redF1PEventId, collision1.globalIndex(), cache);
      // auto groupProton = protontracks.sliceByCached(aod::f1protondaughter::redF1PEventId, collision2.globalIndex(), cache);
      // for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(f1tracks, protontracks))) {
      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupF1, groupProton))) {
        // LOGF(info, "Mixed event collision1 track1: (%d, %d)", collision1.index(), t1.index());
        if (t1.f1MassKaonKshort() > maxKKS0Mass) {
          continue;
        }
        F1.SetXYZM(t1.f1Px(), t1.f1Py(), t1.f1Pz(), t1.f1Mass());
        Pion.SetXYZM(t1.f1d1Px(), t1.f1d1Py(), t1.f1d1Pz(), 0.139);
        Kaon.SetXYZM(t1.f1d2Px(), t1.f1d2Py(), t1.f1d2Pz(), 0.493);
        Kshort.SetXYZM(t1.f1d3Px(), t1.f1d3Py(), t1.f1d3Pz(), 0.497);
        KaonKshortPair = Kaon + Kshort;
        if (Pion.P() > maxMomentumPion || Kaon.P() > maxMomentumKaon) {
          continue;
        }
        if (strategyPIDPion == 1 && Pion.P() > momentumTOFPion && t1.f1d1TOFHit() != 1) {
          continue;
        }
        if (strategyPIDKaon == 1 && Kaon.P() > momentumTOFKaon && t1.f1d2TOFHit() != 1) {
          continue;
        }
        Proton.SetXYZM(t2.protonPx(), t2.protonPy(), t2.protonPz(), 0.938);
        if (Proton.P() < momentumTOFProton && TMath::Abs(t2.protonNsigmaTPC()) > 3) {
          continue;
        }
        if (Proton.P() >= momentumTOFProton && t2.protonTOFHit() != 1 && TMath::Abs(t2.protonNsigmaTOF()) > 3) {
          continue;
        }
        auto relative_momentum = getkstar(F1, Proton);
        if (t1.f1SignalStat() == 1) {
          histos.fill(HIST("h2MixEventInvariantMassUnlike_mass104"), relative_momentum, F1.Pt(), F1.M()); // F1 sign = 1 unlike, F1 sign = -1 like
          histos.fill(HIST("hNsigmaProtonTPCME"), t2.protonNsigmaTPC(), relative_momentum);
        }
        if (t1.f1SignalStat() == -1) {
          histos.fill(HIST("h2MixEventInvariantMassLike_mass104"), relative_momentum, F1.Pt(), F1.M());
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < 9; nrotbkg++) {
            auto anglestart = 5.0 * TMath::Pi() / 6.0;
            auto angleend = 7.0 * TMath::Pi() / 6.0;
            auto anglestep = (angleend - anglestart) / (1.0 * (9.0 - 1.0));
            auto rotangle = anglestart + nrotbkg * anglestep;
            auto rotPionPx = t1.f1d1Px() * std::cos(rotangle) - t1.f1d1Py() * std::sin(rotangle);
            auto rotPionPy = t1.f1d1Px() * std::sin(rotangle) + t1.f1d1Py() * std::cos(rotangle);
            PionRot.SetXYZM(rotPionPx, rotPionPy, t1.f1d1Pz(), 0.139);
            F1Rot = PionRot + KaonKshortPair;
            auto relative_momentum_rot = getkstar(F1Rot, Proton);
            if (t1.f1SignalStat() == 1) {
              histos.fill(HIST("h2MixEventInvariantMassRot_mass"), relative_momentum_rot, F1Rot.Pt(), F1Rot.M());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(f1protoncorrelation, processME, "Process EventMixing for combinatorial background", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<f1protoncorrelation>(cfgc)}; }
