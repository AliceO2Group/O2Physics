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

#include <bitset>
#include <string>
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsFT0/Digit.h"
#include "CCDB/BasicCCDBManager.h"
#include "CollisionTypeHelper.h"
#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

struct ft0CorrectedTable {
  // Configurables
  Configurable<float> resoFT0A{"resoFT0A", 20.f, "FT0A resolution in ps for the MC override"};
  Configurable<float> resoFT0C{"resoFT0C", 20.f, "FT0C resolution in ps for the MC override"};
  Configurable<bool> addHistograms{"addHistograms", false, "Add QA histograms"};
  Configurable<int> cfgCollisionSystem{"collisionSystem", -2, "Collision system: -2 (use cfg values), -1 (autoset), 0 (pp), 1 (PbPb), 2 (XeXe), 3 (pPb)"};
  Configurable<std::string> cfgUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgPathGrpLhcIf{"ccdb-path-grplhcif", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
  Configurable<int64_t> cfgTimestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Producer
  Produces<o2::aod::FT0sCorrected> table;
  using BCsWithMatchings = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  static constexpr float invLightSpeedCm2NS = 1.f / o2::constants::physics::LightSpeedCm2NS;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    if (doprocessStandard && doprocessWithBypassFT0timeInMC) {
      LOG(fatal) << "Both processStandard and processWithBypassFT0timeInMC are enabled. Pick one of the two";
    }
    if (!doprocessStandard && !doprocessWithBypassFT0timeInMC) {
      LOG(fatal) << "No process is enabled. Pick one";
    }
    ccdb->setURL(cfgUrl);
    ccdb->setTimestamp(cfgTimestamp);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    if (doprocessWithBypassFT0timeInMC) {
      // From ps to ns
      resoFT0A.value = resoFT0A.value / 1000.f;
      resoFT0C.value = resoFT0C.value / 1000.f;
    }
    if (!addHistograms) {
      return;
    }
    histos.add("t0A", "t0A", kTH1D, {{1000, -1, 1, "t0A (ns)"}});
    histos.add("t0C", "t0C", kTH1D, {{1000, -1, 1, "t0C (ns)"}});
    histos.add("t0AC", "t0AC", kTH1D, {{1000, -1000, 1000, "t0AC (ns)"}});
    histos.add("deltat0AC", "deltat0AC", kTH1D, {{1000, -10, 10, "#Deltat0AC (ns)"}});
    if (doprocessWithBypassFT0timeInMC) {
      histos.add("MC/deltat0A", "t0A", kTH1D, {{1000, -50, 50, "t0A (ps)"}});
      histos.add("MC/deltat0C", "t0C", kTH1D, {{1000, -50, 50, "t0C (ps)"}});
      histos.add("MC/deltat0AC", "t0AC", kTH1D, {{1000, -50, 50, "t0AC (ps)"}});
    }
  }

  void processStandard(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                       BCsWithMatchings const&,
                       aod::FT0s const&)
  {
    table.reserve(collisions.size());
    float t0A = 1e10f;
    float t0C = 1e10f;
    for (auto& collision : collisions) {
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
        histos.fill(HIST("t0AC"), (t0A + t0C) * 0.5f);
        histos.fill(HIST("deltat0AC"), t0A - t0C);
      }
      table(t0A, t0C);
    }
  }
  PROCESS_SWITCH(ft0CorrectedTable, processStandard, "Process standard table (default)", true);

  void processWithBypassFT0timeInMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
                                    soa::Join<BCsWithMatchings, aod::Timestamps> const& bcs,
                                    aod::FT0s const&,
                                    aod::McCollisions const&)
  {
    if (cfgCollisionSystem.value == -1) {
      o2::parameters::GRPLHCIFData* grpo = ccdb->template getForTimeStamp<o2::parameters::GRPLHCIFData>(cfgPathGrpLhcIf,
                                                                                                        bcs.iteratorAt(0).timestamp());
      cfgCollisionSystem.value = CollisionSystemType::getCollisionTypeFromGrp(grpo);
      switch (cfgCollisionSystem.value) {
        case CollisionSystemType::kCollSyspp:
          resoFT0A.value = 24.f;
          resoFT0C.value = 24.f;
          break;
        case CollisionSystemType::kCollSysPbPb:
          resoFT0A.value = 5.65f;
          resoFT0C.value = 5.65f;
          break;
        default:
          break;
      }
      // Resolution is given in ps
      resoFT0A.value = resoFT0A.value / 1000.f;
      resoFT0C.value = resoFT0C.value / 1000.f;
    }
    table.reserve(collisions.size());
    float t0A = 1e10f;
    float t0C = 1e10f;
    float eventtimeMC = 1e10f;
    float posZMC = 0;
    bool hasMCcoll = false;

    for (auto& collision : collisions) {
      hasMCcoll = false;
      eventtimeMC = 1e10f;
      t0A = 1e10f;
      t0C = 1e10f;
      posZMC = 0;
      const float vertexPV = collision.posZ();
      const float vertex_corr = vertexPV * invLightSpeedCm2NS;
      constexpr float dummyTime = 30.; // Due to HW limitations time can be only within range (-25,25) ns, dummy time is around 32 ns
      if (collision.has_mcCollision()) {
        hasMCcoll = true;
        const auto& collisionMC = collision.mcCollision();
        eventtimeMC = collisionMC.t();
        posZMC = collisionMC.posZ();
      }
      if (collision.has_foundFT0()) {
        const auto& ft0 = collision.foundFT0();
        const std::bitset<8>& triggers = ft0.triggerMask();
        const bool ora = triggers[o2::ft0::Triggers::bitA];
        const bool orc = triggers[o2::ft0::Triggers::bitC];

        if (ora && ft0.timeA() < dummyTime) {
          t0A = ft0.timeA();
          if (hasMCcoll) {
            const float diff = eventtimeMC - posZMC * invLightSpeedCm2NS + gRandom->Gaus(0.f, resoFT0A);
            t0A = diff;
          }
          t0A += vertex_corr;
        }
        if (orc && ft0.timeC() < dummyTime) {
          t0C = ft0.timeC();
          if (hasMCcoll) {
            const float diff = eventtimeMC + posZMC * invLightSpeedCm2NS + gRandom->Gaus(0.f, resoFT0C);
            t0C = diff;
          }
          t0C -= vertex_corr;
        }
      }
      LOGF(debug, " T0 collision time T0A = %f, T0C = %f", t0A, t0C);
      if (addHistograms) {
        histos.fill(HIST("t0A"), t0A);
        histos.fill(HIST("t0C"), t0C);
        histos.fill(HIST("t0AC"), (t0A + t0C) * 0.5f);
        histos.fill(HIST("deltat0AC"), t0A - t0C);
        if (hasMCcoll) {
          histos.fill(HIST("MC/deltat0A"), (t0A - eventtimeMC) * 1000.f);
          histos.fill(HIST("MC/deltat0C"), (t0C - eventtimeMC) * 1000.f);
          histos.fill(HIST("MC/deltat0AC"), ((t0A + t0C) * 0.5f - eventtimeMC) * 1000.f);
        }
      }
      table(t0A, t0C);
    }
  }
  PROCESS_SWITCH(ft0CorrectedTable, processWithBypassFT0timeInMC, "Process MC with bypass of the AO2D information. Use with care!", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ft0CorrectedTable>(cfgc)}; }
