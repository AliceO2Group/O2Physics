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

struct alice3Efficiency {
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {211}, "List of PDG codes to consider for efficiency calculation"};
  OutputObj<THashList> outList{"output"};
  void init(o2::framework::InitContext&)
  {
    outList.setObject(new THashList);
    for (auto pdg : pdgCodes.value) {
      effVsPt[pdg] = new TEfficiency(Form("efficiency_pdg%d", pdg),
                                     Form("Efficiency for PDG %d; p_{T} (GeV/c); Efficiency", pdg),
                                     100, 0, 10);
      outList->Add(effVsPt[pdg]);
    }
  }

  void process(soa::Join<aod::Tracks, o2::aod::McTrackLabels> const& tracks,
               aod::McParticles const& mcParticles)
  {
    std::map<int, std::vector<int64_t>> pdgIndices;
    for (auto& mc : mcParticles) {
      if (effVsPt.find(mc.pdgCode()) == effVsPt.end()) {
        continue;
      }
      pdgIndices[mc.pdgCode()].push_back(mc.globalIndex());
    }
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      // Check that the PDG is in the list of PDGs
      if (pdgIndices.find(mcParticle.pdgCode()) == pdgIndices.end()) {
        continue;
      }

      std::vector<int64_t>& indices = pdgIndices[mcParticle.pdgCode()];
      // Fill efficiency histogram
      auto& eff = effVsPt[mcParticle.pdgCode()];
      const bool found = std::find(indices.begin(), indices.end(), track.globalIndex()) != indices.end();
      eff->Fill(track.pt(), found);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<alice3Efficiency>(ctx)};
}
