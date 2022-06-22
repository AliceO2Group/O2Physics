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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <TH1D.h>
#include <TH2D.h>
#include "TLorentzVector.h"

#include "PWGUD/DataModel/UDTables.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define mmuon 0.1056583745
#define muonPDG 13

struct FwdDimuonsAnalysis {
  HistogramRegistry registry{
    "registry",
    {{"BCs/MuonTrackPairs", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     {"BCs/FT0Signals", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     //
     {"TracksWithFT0/Eta", ";#eta;", {HistType::kTH1D, {{80, -4., 4.}}}},
     {"TracksWithFT0/TimeFT0A", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     {"TracksWithFT0/TimeFT0C", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     //
     {"PairsMC/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/MassRaw", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtRaw", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassMC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/Pt", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/MassNoFT0", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0MC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtNoFT0", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}}}};

  void process(o2::aod::EventCandidates const& eventCandidates,
               soa::Join<o2::aod::SkimmedMuons, o2::aod::SkimmedMuonsExtra, o2::aod::SkimmedMuonTrackLabels> const& muonTracks,
               o2::aod::SkimmedMCEvents const& mcEvents,
               o2::aod::SkimmedMCParticles const& mcParticles)
  {
    int32_t nMCEvents = mcEvents.size();
    for (int32_t i = 0; i < nMCEvents; i++) {
      auto mcPartsGroup = mcParticles.sliceBy(o2::aod::skimmcpart::skimmedMCEventId, i);
      std::vector<aod::SkimmedMCParticle> mcPartsFiltered;
      for (const auto& mcPart : mcPartsGroup) {
        if (!mcPart.isPhysicalPrimary() || std::abs(mcPart.pdgCode()) != muonPDG) {
          continue;
        }
        mcPartsFiltered.emplace_back(mcPart);
      }
      if (mcPartsFiltered.size() != 2) {
        continue;
      }
      const auto& part1 = mcPartsFiltered[0];
      const auto& part2 = mcPartsFiltered[1];
      TLorentzVector p1, p2, p;
      p1.SetXYZM(part1.px(), part1.py(), part1.pz(), mmuon);
      p2.SetXYZM(part2.px(), part2.py(), part2.pz(), mmuon);
      p = p1 + p2;
      registry.fill(HIST("PairsMC/Mass"), p.M());
      mcPartsFiltered.clear();
    }

    for (const auto& cand : eventCandidates) {
      const auto& trackIDs = cand.matchedFwdTracksIds();
      const auto& tr1 = muonTracks.iteratorAt(trackIDs[0]);
      const auto& tr2 = muonTracks.iteratorAt(trackIDs[1]);
      TLorentzVector p1, p2, pPair;
      p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mmuon);
      p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mmuon);
      pPair = p1 + p2;
      // distributions without selection
      registry.fill(HIST("PairSelection/MassRaw"), pPair.M());
      registry.fill(HIST("PairSelection/PtRaw"), pPair.Pt());
      //
      if (tr1.sign() * tr2.sign() >= 0) {
        continue;
      }
      if (pPair.Pt() > 0.25) {
        continue;
      }
      //
      int mcPartId1 = tr1.mcParticleId();
      int mcPartId2 = tr2.mcParticleId();
      const auto& mcPart1 = mcParticles.iteratorAt(mcPartId1);
      const auto& mcPart2 = mcParticles.iteratorAt(mcPartId2);
      TLorentzVector mcP1, mcP2, pPairMC;
      mcP1.SetXYZM(mcPart1.px(), mcPart1.py(), mcPart1.pz(), mmuon);
      mcP2.SetXYZM(mcPart2.px(), mcPart2.py(), mcPart2.pz(), mmuon);
      pPairMC = mcP1 + mcP2;
      //
      registry.fill(HIST("PairSelection/Mass"), pPair.M());
      registry.fill(HIST("PairSelection/MassMC"), pPairMC.M());
      registry.fill(HIST("PairSelection/Pt"), pPair.Pt());
      //
      bool passFT0;
      bool hasFT0 = cand.hasFT0();
      if (hasFT0) {
        registry.fill(HIST("TracksWithFT0/Eta"), p1.Eta());
        registry.fill(HIST("TracksWithFT0/Eta"), p2.Eta());
        registry.fill(HIST("TracksWithFT0/TimeFT0A"), cand.timeAFT0());
        registry.fill(HIST("TracksWithFT0/TimeFT0C"), cand.timeCFT0());
        passFT0 = cand.timeCFT0() > -1. && cand.timeCFT0() < 1.;
      } else {
        passFT0 = true;
      }
      if (passFT0) {
        registry.fill(HIST("PairSelection/MassNoFT0"), pPair.M());
        registry.fill(HIST("PairSelection/MassNoFT0MC"), pPairMC.M());
        registry.fill(HIST("PairSelection/PtNoFT0"), pPair.Pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FwdDimuonsAnalysis>(cfgc)};
}
