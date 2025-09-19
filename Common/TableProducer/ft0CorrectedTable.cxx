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

#include "Common/Core/CollisionTypeHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>

#include <bitset>
#include <chrono>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

struct ft0CorrectedTable {
  // Configurables
  Configurable<bool> addHistograms{"addHistograms", false, "Add QA histograms"};

  // Producer
  Produces<o2::aod::FT0sCorrected> table;
  using BCsWithMatchings = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  static constexpr float invLightSpeedCm2NS = 1.f / o2::constants::physics::LightSpeedCm2NS;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    if (!addHistograms) {
      return;
    }
    histos.add("t0A", "t0A", kTH1D, {{1000, -1, 1, "t0A (ns)"}});
    histos.add("t0C", "t0C", kTH1D, {{1000, -1, 1, "t0C (ns)"}});
    histos.add("t0AC", "t0AC", kTH1D, {{1000, -1000, 1000, "t0AC (ns)"}});
    histos.add("deltat0AC", "deltat0AC", kTH1D, {{1000, -1, 1, "#Deltat0AC (ns)"}});
    histos.add("deltat0ACps", "deltat0ACps", kTH1D, {{1000, -1000, 1000, "#Deltat0AC (ps)"}});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
               BCsWithMatchings const&,
               aod::FT0s const&)
  {
    table.reserve(collisions.size());
    float t0A = 1e10f;
    float t0C = 1e10f;
    for (const auto& collision : collisions) {
      t0A = 1e10f;
      t0C = 1e10f;
      const float vertexPV = collision.posZ();
      const float vertex_corr = vertexPV * invLightSpeedCm2NS;
      constexpr float dummyTime = 30.; // Due to HW limitations time can be only within range (-25,25) ns, dummy time is around 32 ns
      if (collision.has_foundFT0()) {
        const auto& ft0 = collision.foundFT0();
        const std::bitset<8>& triggers = ft0.triggerMask();
        const bool ora = triggers[o2::ft0::Triggers::bitA];
        const bool orc = triggers[o2::ft0::Triggers::bitC];
        LOGF(debug, "triggers OrA %i OrC %i ", ora, orc);
        LOGF(debug, " T0A = %f, T0C %f, vertex_corr %f", ft0.timeA(), ft0.timeC(), vertex_corr);
        if (ora && ft0.timeA() < dummyTime) {
          t0A = ft0.timeA() + vertex_corr;
        }
        if (orc && ft0.timeC() < dummyTime) {
          t0C = ft0.timeC() - vertex_corr;
        }
      }
      LOGF(debug, " T0 collision time T0A = %f, T0C = %f", t0A, t0C);
      if (addHistograms) {
        histos.fill(HIST("t0A"), t0A);
        histos.fill(HIST("t0C"), t0C);
        if (t0A < 1e10f && t0C < 1e10f) {
          histos.fill(HIST("t0AC"), (t0A + t0C) * 0.5f);
          histos.fill(HIST("deltat0AC"), (t0A - t0C) * 0.5f);
          histos.fill(HIST("deltat0ACps"), (t0A - t0C) * 500.f);
        }
      }
      table(t0A, t0C);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ft0CorrectedTable>(cfgc)}; }
