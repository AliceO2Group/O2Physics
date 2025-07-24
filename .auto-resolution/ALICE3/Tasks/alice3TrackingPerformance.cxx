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
#include <vector>

using namespace o2;
using namespace o2::framework;
std::map<int, std::shared_ptr<TH2>> ptResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> invPtResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> dcaXyResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> dcaZResolutionVsPt;

struct alice3TrackingPerformance {
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {211}, "List of PDG codes to consider for efficiency calculation"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::pair<float, float>> etaRange{"etaRange", {-5.f, 5.f}, "Eta range for efficiency calculation"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisPt{100, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec axisPtDelta{100, -1, 1, "p_{T}^{gen} - p_{T}^{reco} (GeV/c)"};
    const AxisSpec axisInvPtDelta{100, -1, 1, "1./p_{T}^{gen} - 1./p_{T}^{reco} (GeV/c)^{-1}"};
    const AxisSpec axisDcaXy{100, -1, 1, "DCA_{xy} (cm)"};
    const AxisSpec axisDcaZ{100, -1, 1, "DCA_{z} (cm)"};
    for (auto pdg : pdgCodes.value) {
      ptResolutionVsPt[pdg] = histos.add<TH2>(Form("ptResolutionVsPt_%d", pdg), "", kTH2D, {axisPt, axisPtDelta});
      invPtResolutionVsPt[pdg] = histos.add<TH2>(Form("invPtResolutionVsPt_%d", pdg), "", kTH2D, {axisPt, axisInvPtDelta});
      dcaXyResolutionVsPt[pdg] = histos.add<TH2>(Form("dcaXyResolutionVsPt_%d", pdg), "", kTH2D, {axisPt, axisDcaXy});
      dcaZResolutionVsPt[pdg] = histos.add<TH2>(Form("dcaZResolutionVsPt_%d", pdg), "", kTH2D, {axisPt, axisDcaZ});
    }
  }

  void process(soa::Join<aod::Tracks, o2::aod::McTrackLabels, o2::aod::TracksDCA> const& tracks,
               aod::McParticles const&)
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
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      if (!isParticleSelected(mcParticle)) {
        continue;
      }
      if (ptResolutionVsPt.find(mcParticle.pdgCode()) == ptResolutionVsPt.end()) {
        continue;
      }
      ptResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), mcParticle.pt() - track.pt());
      invPtResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), 1.f / mcParticle.pt() - 1.f / track.pt());
      dcaXyResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), track.dcaXY());
      dcaZResolutionVsPt[mcParticle.pdgCode()]->Fill(mcParticle.pt(), track.dcaZ());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<alice3TrackingPerformance>(ctx)};
}
