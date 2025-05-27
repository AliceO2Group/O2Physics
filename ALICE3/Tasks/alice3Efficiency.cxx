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
/// \file alice3Efficiency.cxx
///
/// \brief This task produces the efficiency
///
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \since  May 27, 2025
///

#include <map>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ConfigParamRegistry.h"
#include "TEfficiency.h"
#include "THashList.h"

using namespace o2;
using namespace o2::framework;
std::map<int, TEfficiency*> effVsPt;
std::map<int, TEfficiency*> effVsEta;

struct alice3Efficiency {
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {211}, "List of PDG codes to consider for efficiency calculation"};
  OutputObj<THashList> outList{"output"};
  Configurable<std::pair<float, float>> etaRange{"etaRange", {-5.f, 5.f}, "Eta range for efficiency calculation"};
  void init(o2::framework::InitContext&)
  {
    outList.setObject(new THashList);
    auto createEff = [&](const char* baseName, const char* axisTitle, int pdg, int nBins, double min, double max) {
      auto eff = new TEfficiency(Form("%s_pdg%d", baseName, pdg),
                                 Form("Efficiency for PDG %d; %s; Efficiency", pdg, axisTitle),
                                 nBins, min, max);
      outList->Add(eff);
      return eff;
    };
    for (auto pdg : pdgCodes.value) {
      effVsPt[pdg] = createEff("efficiency", "p_{T} (GeV/c)", pdg, 100, 0, 10);
      effVsEta[pdg] = createEff("efficiency_eta", "#eta", pdg, 100, -5, 5);
    }
  }

  void process(soa::Join<aod::Tracks, o2::aod::McTrackLabels> const& tracks,
               aod::McParticles const& mcParticles)
  {
    std::map<int, std::vector<int64_t>> pdgIndices;

    // Lambda function to select particles after all cuts
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
      pdgIndices[mcParticle.pdgCode()].push_back(mcParticle.globalIndex());
    }

    for (auto& mc : mcParticles) {
      if (effVsPt.find(mc.pdgCode()) == effVsPt.end()) {
        continue;
      }
      if (!isParticleSelected(mc)) {
        continue;
      }
      std::vector<int64_t>& indices = pdgIndices[mc.pdgCode()];
      // Fill efficiency histogram
      const bool found = std::find(indices.begin(), indices.end(), mc.globalIndex()) != indices.end();
      effVsPt[mc.pdgCode()]->Fill(found, mc.pt());
      effVsEta[mc.pdgCode()]->Fill(found, mc.eta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<alice3Efficiency>(ctx)};
}
