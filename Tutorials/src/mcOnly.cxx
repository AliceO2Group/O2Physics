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
/// \brief Accessing MC data and the related MC truth.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace constants::math;

// Simple access to collision
struct VertexDistribution {
  OutputObj<TH1F> vertex{TH1F("vertex", "vertex", 100, -10, 10)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // loop over MC truth McCollisions
  void process(aod::McCollision const& mcCollision)
  {
    if (reduceOutput < 2) {
      LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
    }
    vertex->Fill(mcCollision.posZ());
  }
};

// Grouping between MC particles and collisions
struct AccessMcData {
  OutputObj<TH1F> phiH{TH1F("phi", "phi", 100, 0., TwoPI)};
  OutputObj<TH1F> etaH{TH1F("eta", "eta", 102, -2.01, 2.01)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // group according to McCollisions
  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    if (reduceOutput < 2) {
      LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
      LOGF(info, "First: %d | Length: %d", mcParticles.begin().index(), mcParticles.size());
    }
    int count = 0;
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        phiH->Fill(mcParticle.phi());
        etaH->Fill(mcParticle.eta());
        count++;
        // Loop over mothers and daughters
        if (mcParticle.has_mothers()) {
          // Check first mother
          auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
          if (reduceOutput == 0) {
            LOGF(info, "First mother: %d has pdg code %d", mother.globalIndex(), mother.pdgCode());
          }
          // Loop over all mothers (needed for some MCs with junctions etc.)
          for (auto& m : mcParticle.mothers_as<aod::McParticles>()) {
            LOGF(debug, "M2 %d %d", mcParticle.globalIndex(), m.globalIndex());
          }
        }
        if (mcParticle.has_daughters()) {
          for (auto& d : mcParticle.daughters_as<aod::McParticles>()) {
            LOGF(debug, "D2 %d %d", mcParticle.globalIndex(), d.globalIndex());
          }
        }
      }
    }
    if (reduceOutput < 2) {
      LOGF(info, "Primaries for this collision: %d", count);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<VertexDistribution>(cfgc),
    adaptAnalysisTask<AccessMcData>(cfgc)};
}
