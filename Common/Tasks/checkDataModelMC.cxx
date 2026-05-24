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

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <limits>
#include <memory>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Check the complete list of indices of final-state daughters of an MC particle.
/// \param particlesMC  table with MC particles
/// \param offset       offset in case of grouped tables
/// \param particle     MC particle
/// \param debugMode    flag for debug mode
/// \param debugHisto   histo to fill for debug
template <soa::is_table T>
void checkDaughters(const T& particlesMC,
                    const typename T::iterator& particle,
                    uint64_t offset,
                    bool debugMode,
                    std::shared_ptr<TH1>& debugHisto)
{
  auto firstDauIdx = particle.daughtersIds().front();
  auto lastDauIdx = particle.daughtersIds().back();
  if ((firstDauIdx < 0 && lastDauIdx >= 0) || (lastDauIdx < 0 && firstDauIdx >= 0)) {
    if (debugMode) {
      debugHisto->Fill(0);
    } else {
      LOG(fatal) << "MC particle " << particle.globalIndex() << " with PDG " << particle.pdgCode() << " has first and last daughter indices " << firstDauIdx << ", " << lastDauIdx;
    }
  }
  for (const auto& idxDau : particle.daughtersIds()) {
    if (idxDau >= 0 && ((uint64_t)idxDau > offset + particlesMC.size() || (uint64_t)idxDau < offset)) {
      if (debugMode) {
        debugHisto->Fill(1);
      } else {
        LOG(fatal) << "MC particle " << particle.globalIndex() << " with PDG " << particle.pdgCode() << " has daughter with index " << idxDau << " > MC particle table size (+ offset) (" << particlesMC.size() + offset << ")";
      }
    } else if (idxDau < -1 * std::numeric_limits<int>::max() / 2) {
      if (debugMode) {
        debugHisto->Fill(2);
      } else {
        LOG(fatal) << "MC particle " << particle.globalIndex() << " with PDG " << particle.pdgCode() << " has daughter with too negative index " << idxDau;
      }
    }
  }
}

/// Check for cycles in the decay chain
/// \param particlesMC  table with MC particles
template <soa::is_table T>
bool checkCycles(T const& particlesMC)
{
  auto slow = particlesMC.begin();
  auto fast = particlesMC.begin();
  auto probe = particlesMC.begin();
  bool hasCycle = false;
  for (auto start = 0U; start < particlesMC.size(); ++start) {
    probe.setCursor(start);
    if (!probe.has_mothers() || probe.mothersIds()[0] == -1) {
      continue;
    }
    slow.setCursor(start);
    fast.setCursor(start);

    while (true) {
      if (!slow.has_mothers() || slow.mothersIds()[0] < 0) {
        break;
      }
      slow.setCursor(slow.mothersIds()[0]);

      if (!fast.has_mothers() || fast.mothersIds()[0] < 0) {
        break;
      }
      fast.setCursor(fast.mothersIds()[0]);
      if (!fast.has_mothers() || fast.mothersIds()[0] < 0) {
        break;
      }
      fast.setCursor(fast.mothersIds()[0]);

      if (slow.globalIndex() == fast.globalIndex()) {
        LOGP(error, "Cycle found: particle: {} | cycle starts at: {}", start, slow.globalIndex());
        hasCycle = true;
        break;
      }
    }
  }
  return hasCycle;
}

template <soa::is_table Table>
struct LoadTable {
  OutputObj<TH1F> counter{TH1F("counter", "counter", 2, 0., 2)};
  void process(Table const& table)
  {
    LOGF(info, "Table has %d entries", table.size());
    counter->Fill(0.5);
    counter->Fill(1.5, table.size());
  }
};

std::shared_ptr<TH1> defineHistogram(HistogramRegistry& registry)
{
  std::shared_ptr<TH1> hDebug;
  hDebug = registry.add<TH1>("hDebug", "debug histo;;counts", HistType::kTH1F, {{3, -0.5, 2.5}});
  hDebug->GetXaxis()->SetBinLabel(1, "(-X, Y) or (X, -Y)");
  hDebug->GetXaxis()->SetBinLabel(2, "out of range");
  hDebug->GetXaxis()->SetBinLabel(3, "#minus max integer");
  return hDebug;
}

struct CheckMcParticlesIndices {

  Configurable<bool> debugMode{"debugMode", false, "flag to enable debug mode"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  std::shared_ptr<TH1> hDebug;

  void init(o2::framework::InitContext&)
  {
    hDebug = defineHistogram(registry);
  }

  void process(aod::McParticles const& particlesMC)
  {
    if (checkCycles(particlesMC)) {
      LOG(fatal) << "Cycles found, aborting.";
    }
    uint64_t offset = 0;
    for (const auto& particle : particlesMC) {
      checkDaughters(particlesMC, particle, offset, debugMode.value, hDebug);
    }
  }
};

struct CheckMcParticlesIndicesGrouped {

  Configurable<bool> debugMode{"debugMode", false, "flag to enable debug mode"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  std::shared_ptr<TH1> hDebug;

  void init(o2::framework::InitContext&)
  {
    hDebug = defineHistogram(registry);
  }

  void process(aod::McCollision const&,
               aod::McParticles const& particlesMC)
  {
    for (const auto& particle : particlesMC) {
      checkDaughters(particlesMC, particle, particlesMC.offset(), debugMode.value, hDebug);
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
