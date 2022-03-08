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
/// \brief Load all tables in the data model to check if they can be read correctly
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Check the complete list of indices of final-state daughters of an MC particle.
/// \param particlesMC  table with MC particles
/// \param particle  MC particle
template <typename T>
void checkDaughters(const T& particlesMC,
                    const typename T::iterator& particle)
{
  auto firstDauIdx = particle.daughtersIds().front();
  auto lastDauIdx = particle.daughtersIds().back();
  if ((firstDauIdx < 0 && lastDauIdx >= 0) || (lastDauIdx < 0 && firstDauIdx >= 0)) {
    LOG(fatal) << "MC particle " << particle.globalIndex() << " with PDG " << particle.pdgCode() << " has first and last daughter indices " << firstDauIdx << ", " << lastDauIdx;
  }
  for (auto& idxDau : particle.daughtersIds()) {
    if (idxDau > particlesMC.size()) {
      LOG(fatal) << "MC particle " << particle.globalIndex() << " with PDG " << particle.pdgCode() << " has daughter with index " << idxDau << " > MC particle table size (" << particlesMC.size() << ")";
    }
  }
}

template <typename Table>
struct LoadTable {
  OutputObj<TH1F> counter{TH1F("counter", "counter", 2, 0., 2)};
  void process(Table const& table)
  {
    LOGF(info, "Table has %d entries", table.size());
    counter->Fill(0.5);
    counter->Fill(1.5, table.size());
  }
};

struct CheckMcParticlesIndices {
  void process(aod::McParticles const& particlesMC)
  {
    for (auto& particle : particlesMC) {
      checkDaughters(particlesMC, particle);
    }
  }
};

struct CheckMcParticlesIndicesGrouped {
  void process(aod::McCollision const& collision,
               aod::McParticles const& particlesMC)
  {
    for (auto& particle : particlesMC) {
      checkDaughters(particlesMC, particle);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LoadTable<aod::McParticles>>(cfgc, TaskName("McParticles")),
    adaptAnalysisTask<LoadTable<aod::McCollisions>>(cfgc, TaskName("McCollisions")),
    adaptAnalysisTask<LoadTable<aod::McTrackLabels>>(cfgc, TaskName("McTrackLabels")),
    // adaptAnalysisTask<LoadTable<aod::McFwdTrackLabels>>(cfgc, TaskName("McFwdTrackLabels")),
    // adaptAnalysisTask<LoadTable<aod::McMFTTrackLabels>>(cfgc, TaskName("McMFTTrackLabels")),
    adaptAnalysisTask<CheckMcParticlesIndices>(cfgc),
    adaptAnalysisTask<CheckMcParticlesIndicesGrouped>(cfgc)};
}
