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
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;

struct systematicsStudy {
  ConfigurableAxis ptBins{"ptBins", {100, 0.f, 10.f}, "Binning for pT (GeV/c)"};
  Configurable<float> cutEta{"cutEta", 0.8f, "Max |eta| for tracks and V0s"};
  Configurable<int> cutTPCClusters{"cutTPCClusters", 70, "Min TPC clusters for tracks"};
  Configurable<float> cutKaonNSigma{"cutKaonNSigma", 3.f, "Max |nSigma| for kaon PID"};
  Configurable<float> cutK0sMassWindow{"cutK0sMassWindow", 0.01f, "K0s mass window (GeV/c^2)"};
  ConfigurableAxis etaBins{"etaBins", {40, -1.0f, 1.0f}, "Binning for #eta"};
  ConfigurableAxis phiBins{"phiBins", {36, 0.f, 2 * M_PI}, "Binning for #phi (rad)"};
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/c)"};
  const AxisSpec etaAxis{etaBins, "#eta"};
  const AxisSpec phiAxis{phiBins, "#phi (rad)"};
  registry.add("hKaonYieldData", "", HistType::kTH1F, {ptAxis});
  registry.add("hKaonYieldMC", "", HistType::kTH1F, {ptAxis});
  registry.add("hK0sYieldData", "", HistType::kTH1F, {ptAxis});
  registry.add("hK0sYieldMC", "", HistType::kTH1F, {ptAxis});
  registry.add("hKaonYieldMapData", "", HistType::kTH3F, {ptAxis, etaAxis, phiAxis});
  registry.add("hKaonYieldMapMC", "", HistType::kTH3F, {ptAxis, etaAxis, phiAxis});
  registry.add("hK0sYieldMapData", "", HistType::kTH3F, {ptAxis, etaAxis, phiAxis});
  registry.add("hK0sYieldMapMC", "", HistType::kTH3F, {ptAxis, etaAxis, phiAxis});
  }

  void processData(aod::Collisions const& collisions,
                   aod::Tracks const& tracks,
                   aod::V0s const& v0s)
  {
    for (auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.posZ()) > 10)
        continue; // MB selection
      if (collision.isMC()) continue;

      // Kaon loop
      for (auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex())
          continue;
        if (std::abs(track.eta()) > cutEta)
          continue;
        if (track.tpcNClsFound() < cutTPCClusters)
          continue;
        // PID selection for kaons
        if (std::abs(track.tpcNSigmaKa()) < cutKaonNSigma) {
          registry.fill(HIST("hKaonYieldData"), track.pt());
          registry.fill(HIST("hKaonYieldMapData"), track.pt(), track.eta(), track.phi());
        }
      }

      // K0s loop
      for (auto& v0 : v0s) {
        if (v0.collisionId() != collision.globalIndex())
          continue;
        // Basic selection for K0s
        if (std::abs(v0.eta()) > cutEta)
          continue;
        if (std::abs(v0.mass() - constants::physics::MassK0Short) > cutK0sMassWindow)
          continue;
        registry.fill(HIST("hK0sYieldData"), v0.pt());
        registry.fill(HIST("hK0sYieldMapData"), v0.pt(), v0.eta(), v0.phi());
      }
    }
  }

  void processMC(aod::Collisions const& collisions,
                 aod::Tracks const& tracks,
                 aod::V0s const& v0s)
  {
    for (auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.posZ()) > 10)
        continue; // MB selection
      if (!collision.isMC()) continue;

      // Kaon loop
      for (auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex())
          continue;
        if (std::abs(track.eta()) > cutEta)
          continue;
        if (track.tpcNClsFound() < cutTPCClusters)
          continue;
        // PID selection for kaons
        if (std::abs(track.tpcNSigmaKa()) < cutKaonNSigma) {
          registry.fill(HIST("hKaonYieldMC"), track.pt());
          registry.fill(HIST("hKaonYieldMapMC"), track.pt(), track.eta(), track.phi());
        }
      }

      // K0s loop
      for (auto& v0 : v0s) {
        if (v0.collisionId() != collision.globalIndex())
          continue;
        // Basic selection for K0s
        if (std::abs(v0.eta()) > cutEta)
          continue;
        if (std::abs(v0.mass() - constants::physics::MassK0Short) > cutK0sMassWindow)
          continue;
        registry.fill(HIST("hK0sYieldMC"), v0.pt());
        registry.fill(HIST("hK0sYieldMapMC"), v0.pt(), v0.eta(), v0.phi());
      }
    }
  }

  void process(aod::Collisions const& collisions,
               aod::Tracks const& tracks,
               aod::V0s const& v0s)
  {
    processData(collisions, tracks, v0s);
    processMC(collisions, tracks, v0s);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<systematicsStudy>(cfgc)};
}
