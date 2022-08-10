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
#include "PWGUD/DataModel/UDTables.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define mmuon 0.1056583745

struct UpcTrackSkimmerQa {
  HistogramRegistry registry{
    "registry",
    {
      {"MC/Kine/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"MC/Kine/Eta", ";#eta;", {HistType::kTH1D, {{100, -4., 4.}}}},
      {"MC/Kine/PtTOFCut", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      //
      {"Muons/Timings/TrackTimeRes", ";Track time resolution, ns;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Muons/Kine/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Muons/Kine/Eta", ";#eta;", {HistType::kTH1D, {{100, -4., 4.}}}},
      //
      {"Barrels/Timings/TrackTimeRes", ";Track time resolution, ns;", {HistType::kTH1D, {{1000, 0., 1000.}}}},
      //
      {"Barrels/Kine/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Barrels/Kine/PtWTOF", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Barrels/Kine/PtMC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Barrels/Kine/PtWTOFMC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Barrels/Kine/Eta", ";#eta;", {HistType::kTH1D, {{100, -4., 4.}}}},
      {"Barrels/ITS/chi2", ";chi2;", {HistType::kTH1D, {{1100, -100., 1000.}}}},
      {"Barrels/ITS/nClusters", ";n clusters;", {HistType::kTH1I, {{10, 0, 10}}}},
      {"Barrels/TPC/chi2", ";chi2;", {HistType::kTH1D, {{110, -10., 100.}}}},
      {"Barrels/TPC/nClusters", ";n clusters;", {HistType::kTH1I, {{200, 0, 200}}}},
      {"Barrels/TPC/nClustersWTOF", ";n clusters;", {HistType::kTH1I, {{200, 0, 200}}}},
      {"Barrels/TOF/chi2", ";chi2;", {HistType::kTH1D, {{1050, -1000., 50.}}}},
      {"Barrels/TOF/ExpMom", ";TOF exp. mom., GeV;", {HistType::kTH1D, {{1000, 0., 100.}}}},
      {"Barrels/MC/Labels", ";MC labels;", {HistType::kTH1I, {{1001, -1, 1000}}}}
    }};

  void process(soa::Join<o2::aod::UDFwdTracks, o2::aod::UDMcFwdTrackLabels> const& muonTracks,
               soa::Join<o2::aod::UDTracks, o2::aod::UDTracksExtra, o2::aod::UDMcTrackLabels> const& barTracks,
               o2::aod::UDMcParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      TLorentzVector p;
      p.SetXYZM(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mmuon);
      registry.fill(HIST("MC/Kine/Pt"), p.Pt());
      registry.fill(HIST("MC/Kine/Eta"), p.Eta());
      if (std::abs(p.Eta()) < 0.8 && std::abs(p.Pt()) > 0.3) {
        registry.fill(HIST("MC/Kine/PtTOFCut"), p.Pt());
      }
    }

    for (const auto& track : muonTracks) {
      int32_t mcPartId = track.udMcParticleId();
      if (mcPartId < 0) {
        continue;
      }
      const auto& mcPart = mcParticles.iteratorAt(mcPartId);
      if (!(mcPart.isPhysicalPrimary() && std::abs(mcPart.pdgCode()) == 13)) {
        continue;
      }
      TLorentzVector p;
      p.SetXYZM(track.px(), track.py(), track.pz(), mmuon);
      registry.fill(HIST("Muons/Timings/TrackTimeRes"), track.trackTimeRes());
      registry.fill(HIST("Muons/Kine/Pt"), p.Pt());
      registry.fill(HIST("Muons/Kine/Eta"), p.Eta());
    }

    for (const auto& track : barTracks) {
      int32_t mcPartId = track.udMcParticleId();
      if (mcPartId < 0) {
        continue;
      }
      const auto& mcPart = mcParticles.iteratorAt(mcPartId);
      if (!(mcPart.isPhysicalPrimary() && std::abs(mcPart.pdgCode()) == 13)) {
        continue;
      }
      registry.fill(HIST("Barrels/MC/Labels"), mcPartId);
      TLorentzVector p;
      p.SetXYZM(track.px(), track.py(), track.pz(), mmuon);
      TLorentzVector pMC;
      pMC.SetXYZM(mcPart.px(), mcPart.py(), mcPart.pz(), mmuon);
      if (std::abs(pMC.Eta()) < 0.8 && std::abs(pMC.Pt()) > 0.3 && track.hasTOF()) {
        registry.fill(HIST("Barrels/Kine/PtWTOF"), p.Pt());
        registry.fill(HIST("Barrels/Kine/PtWTOFMC"), pMC.Pt());
        registry.fill(HIST("Barrels/TPC/nClustersWTOF"), track.tpcNClsCrossedRows());
      }
      registry.fill(HIST("Barrels/Timings/TrackTimeRes"), track.trackTimeRes());
      registry.fill(HIST("Barrels/Kine/Pt"), p.Pt());
      registry.fill(HIST("Barrels/Kine/PtMC"), pMC.Pt());
      registry.fill(HIST("Barrels/Kine/Eta"), p.Eta());
      registry.fill(HIST("Barrels/ITS/chi2"), track.itsChi2NCl());
      registry.fill(HIST("Barrels/ITS/nClusters"), track.itsNCls());
      registry.fill(HIST("Barrels/TPC/chi2"), track.tpcChi2NCl());
      registry.fill(HIST("Barrels/TPC/nClusters"), track.tpcNClsCrossedRows());
      registry.fill(HIST("Barrels/TOF/chi2"), track.tofChi2());
      registry.fill(HIST("Barrels/TOF/ExpMom"), track.tofExpMom());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTrackSkimmerQa>(cfgc)};
}
