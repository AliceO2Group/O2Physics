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
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  May 27, 2025
///

#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TString.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
std::shared_ptr<TH1> particlePdgCodes;
std::map<int, std::shared_ptr<TH1>> particlePtDistribution;
std::map<int, std::shared_ptr<TH1>> particleEtaDistribution;
std::map<int, std::shared_ptr<TH1>> ptDistribution;
std::map<int, std::shared_ptr<TH2>> ptResolutionVsPt;
std::map<int, std::shared_ptr<TProfile2D>> ptResolutionVsEta;
std::map<int, std::shared_ptr<TH2>> invPtResolutionVsPt;
std::map<int, std::shared_ptr<TProfile2D>> invPtResolutionVsEta;
std::map<int, std::shared_ptr<TH2>> dcaXyResolutionVsPt;
std::map<int, std::shared_ptr<TH2>> dcaZResolutionVsPt;

struct Alice3TrackingPerformance {
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {0, 211}, "List of PDG codes to consider for efficiency calculation. (0 means all)"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::pair<float, float>> etaRange{"etaRange", {-5.f, 5.f}, "Eta range for efficiency calculation"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisPt{500, 0, 100, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisEta{100, etaRange.value.first, etaRange.value.second, "#eta"};
    const AxisSpec axisPtDelta{100, -1, 1, "(#it{p}_{T}^{reco} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}"};
    const AxisSpec axisInvPtDelta{100, -1, 1, "1./#it{p}_{T}^{gen} - 1./#it{p}_{T}^{reco} (GeV/#it{c})^{-1}"};
    const AxisSpec axisDcaXy{100, -1, 1, "DCA_{xy} (cm)"};
    const AxisSpec axisDcaZ{100, -1, 1, "DCA_{z} (cm)"};
    particlePdgCodes = histos.add<TH1>("particlePdgCodes", "", kTH1D, {AxisSpec{100, -0.5, 99.5, "PDG Code"}});
    for (const int& pdg : pdgCodes.value) {
      std::string prefix = Form("%i", pdg);
      if (pdg < 0) {
        prefix = Form("m%i", -pdg);
      }
      const std::string tag = "_" + prefix;
      prefix += "/";
      particlePtDistribution[pdg] = histos.add<TH1>(prefix + "particlePtDistribution" + tag, "", kTH1D, {axisPt});
      particleEtaDistribution[pdg] = histos.add<TH1>(prefix + "particleEtaDistribution" + tag, "", kTH1D, {axisEta});

      ptDistribution[pdg] = histos.add<TH1>(prefix + "ptDistribution" + tag, "", kTH1D, {axisPt});
      ptResolutionVsPt[pdg] = histos.add<TH2>(prefix + "ptResolutionVsPt" + tag, "", kTH2D, {axisPt, axisPtDelta});
      ptResolutionVsEta[pdg] = histos.add<TProfile2D>(prefix + "ptResolutionVsEta" + tag, "", kTProfile2D, {axisPt, axisEta});
      invPtResolutionVsPt[pdg] = histos.add<TH2>(prefix + "invPtResolutionVsPt" + tag, "", kTH2D, {axisPt, axisInvPtDelta});
      invPtResolutionVsEta[pdg] = histos.add<TProfile2D>(prefix + "invPtResolutionVsEta" + tag, "", kTProfile2D, {axisPt, axisEta});
      dcaXyResolutionVsPt[pdg] = histos.add<TH2>(prefix + "dcaXyResolutionVsPt" + tag, "", kTH2D, {axisPt, axisDcaXy});
      dcaZResolutionVsPt[pdg] = histos.add<TH2>(prefix + "dcaZResolutionVsPt" + tag, "", kTH2D, {axisPt, axisDcaZ});
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
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      const float ptResolution = (track.pt() - mcParticle.pt()) / mcParticle.pt();
      const float invptResolution = 1.f / track.pt() - 1.f / mcParticle.pt();

      auto fillResolutionHistograms = [&](const int p) {
        ptDistribution[p]->Fill(track.pt());
        ptResolutionVsPt[p]->Fill(mcParticle.pt(), ptResolution);
        ptResolutionVsEta[p]->Fill(mcParticle.pt(), mcParticle.eta(), ptResolution);
        invPtResolutionVsPt[p]->Fill(mcParticle.pt(), invptResolution);
        invPtResolutionVsEta[p]->Fill(mcParticle.pt(), mcParticle.eta(), invptResolution);
        dcaXyResolutionVsPt[p]->Fill(mcParticle.pt(), track.dcaXY());
        dcaZResolutionVsPt[p]->Fill(mcParticle.pt(), track.dcaZ());
      };

      fillResolutionHistograms(0);

      if (!isParticleSelected(mcParticle)) {
        continue;
      }
      if (ptResolutionVsPt.find(mcParticle.pdgCode()) == ptResolutionVsPt.end()) {
        continue;
      }
      fillResolutionHistograms(mcParticle.pdgCode());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3TrackingPerformance>(ctx)};
}
