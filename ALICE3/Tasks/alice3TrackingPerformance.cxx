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
/// \file alice3TrackingPerformance.cxx
///
/// \brief This task produces the tracking performance
///
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \since  May 27, 2025
///

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
std::shared_ptr<TH1> particlePdgCodes;
std::map<int, std::shared_ptr<TH1>> particlePtDistribution;
std::map<int, std::shared_ptr<TH1>> particleEtaDistribution;
std::map<int, std::shared_ptr<TH1>> ptDistribution;
std::map<int, std::shared_ptr<TH2>> ptResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> invPtResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> dcaXyResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> dcaZResolutionVsPt;

struct Alice3TrackingPerformance {
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {0, 211}, "List of PDG codes to consider for efficiency calculation. (0 means all)"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::pair<float, float>> etaRange{"etaRange", {-5.f, 5.f}, "Eta range for efficiency calculation"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisPt{100, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec axisEta{100, etaRange.value.first, etaRange.value.second, "#eta"};
    const AxisSpec axisPtDelta{100, -1, 1, "p_{T}^{gen} - p_{T}^{reco} (GeV/c)"};
    const AxisSpec axisInvPtDelta{100, -1, 1, "1./p_{T}^{gen} - 1./p_{T}^{reco} (GeV/c)^{-1}"};
    const AxisSpec axisDcaXy{100, -1, 1, "DCA_{xy} (cm)"};
    const AxisSpec axisDcaZ{100, -1, 1, "DCA_{z} (cm)"};
    particlePdgCodes = histos.add<TH1>("particlePdgCodes", "", kTH1D, {AxisSpec{100, -0.5, 99.5, "PDG Code"}});
    for (const int& pdg : pdgCodes.value) {
      std::string tag = Form("_%d", pdg);
      if (pdg < 0) {
        tag = Form("_m%d", -pdg);
      }
      particlePtDistribution[pdg] = histos.add<TH1>("particlePtDistribution" + tag, "", kTH1D, {axisPt});
      particleEtaDistribution[pdg] = histos.add<TH1>("particleEtaDistribution" + tag, "", kTH1D, {axisEta});

      ptDistribution[pdg] = histos.add<TH1>("ptDistribution" + tag, "", kTH1D, {axisPt});
      ptResolutionVsPt[pdg] = histos.add<TH2>("ptResolutionVsPt" + tag, "", kTH2D, {axisPt, axisPtDelta});
      invPtResolutionVsPt[pdg] = histos.add<TH2>("invPtResolutionVsPt" + tag, "", kTH2D, {axisPt, axisInvPtDelta});
      dcaXyResolutionVsPt[pdg] = histos.add<TH2>("dcaXyResolutionVsPt" + tag, "", kTH2D, {axisPt, axisDcaXy});
      dcaZResolutionVsPt[pdg] = histos.add<TH2>("dcaZResolutionVsPt" + tag, "", kTH2D, {axisPt, axisDcaZ});
    }
  }

  void process(soa::Join<aod::Tracks, o2::aod::McTrackLabels, o2::aod::TracksDCA> const& tracks,
               aod::McParticles const& mcParticles)
  {
    auto isParticleSelected = [&](const o2::aod::McParticle& p) {
      if (!p.isPhysicalPrimary()) {
        return false;
      }
      if (p.eta() < etaRange.value.first) {
        return false;
      }
      if (p.eta() > etaRange.value.second) {
        return false;
      }
      return true;
    };

    for (const auto& mcParticle : mcParticles) {
      particlePdgCodes->Fill(Form("%d", mcParticle.pdgCode()), 1);
      if (!isParticleSelected(mcParticle)) {
        continue;
      }
      particlePtDistribution[0]->Fill(mcParticle.pt());
      particleEtaDistribution[0]->Fill(mcParticle.eta());
      if (particlePtDistribution.find(mcParticle.pdgCode()) == particlePtDistribution.end()) {
        continue;
      }
      particlePtDistribution[mcParticle.pdgCode()]->Fill(mcParticle.pt());
      particleEtaDistribution[mcParticle.pdgCode()]->Fill(mcParticle.eta());
    }
    for (const auto& track : tracks) {
      ptDistribution[0]->Fill(track.pt());
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      ptResolutionVsPt[0]->Fill(mcParticle.pt(), mcParticle.pt() - track.pt());
      invPtResolutionVsPt[0]->Fill(mcParticle.pt(), 1.f / mcParticle.pt() - 1.f / track.pt());
      dcaXyResolutionVsPt[0]->Fill(mcParticle.pt(), track.dcaXY());
      dcaZResolutionVsPt[0]->Fill(mcParticle.pt(), track.dcaZ());
      if (!isParticleSelected(mcParticle)) {
        continue;
      }
      if (ptResolutionVsPt.find(mcParticle.pdgCode()) == ptResolutionVsPt.end()) {
        continue;
      }
      ptDistribution[mcParticle.pdgCode()]->Fill(mcParticle.pt());
      ptResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), mcParticle.pt() - track.pt());
      invPtResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), 1.f / mcParticle.pt() - 1.f / track.pt());
      dcaXyResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), track.dcaXY());
      dcaZResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), track.dcaZ());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3TrackingPerformance>(ctx)};
}
