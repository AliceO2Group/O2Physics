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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/runDataProcessing.h"

#include "Gencentralities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

constexpr float FT0Alo = -3.3;
constexpr float FT0Ahi = -2.1;
constexpr float FT0Clo = 3.5;
constexpr float FT0Chi = 4.9;

struct Binner {
  Service<o2::framework::O2DatabasePDG> pdg;
  Produces<aod::GenCents> gencents;
  HistogramRegistry h{"histograms", {}};

  using Particles = soa::Filtered<aod::McParticles>;
  Preslice<Particles> perMcCol = aod::mcparticle::mcCollisionId;

  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};

  // The objects are uploaded with https://alimonitor.cern.ch/ccdb/upload.jsp
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> path{"ccdb-path", "Users/a/aalkin/gencentralities", "base path to the ccdb object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  bool isChargedParticle(int code)
  {
    static auto p = pdg->GetParticle(code);
    if (p == nullptr) {
      return 0;
    }
    return std::abs(p->Charge()) >= 3.;
  }

  TH1F* multFT0C = nullptr;
  TH1F* multFT0M = nullptr;

  void init(InitContext const&)
  {
    AxisSpec MultAxis = {multBinning};
    h.add({"hFT0M", "; N_{part} in FT0M acc", {HistType::kTH1F, {MultAxis}}});
    h.add({"hFT0C", "; N_{part} in FT0C acc", {HistType::kTH1F, {MultAxis}}});
    if (dobin) {
      ccdb->setURL(url.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
    } else if (docalibrate) {
    } else {
      LOGP(fatal, "Need to have either calibration or binning enabled");
    }
  }

  Filter primaries = ncheckbit(aod::mcparticle::flags, (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);
  Partition<Particles> pFT0M = ((aod::mcparticle::eta > FT0Alo) && (aod::mcparticle::eta < FT0Ahi)) || ((aod::mcparticle::eta > FT0Clo) && (aod::mcparticle::eta < FT0Chi));
  Partition<Particles> pFT0C = (aod::mcparticle::eta > FT0Clo) && (aod::mcparticle::eta < FT0Chi);

  void calibrate(aod::McCollisions const& mccollisions, Particles const&)
  {
    for (auto& mcc : mccollisions) {
      auto pcFT0M = pFT0M.sliceBy(perMcCol, mcc.globalIndex());
      auto pcFT0C = pFT0C.sliceBy(perMcCol, mcc.globalIndex());
      int nFT0M = 0;
      int nFT0C = 0;
      for (auto& p : pcFT0M) {
        if (isChargedParticle(p.pdgCode())) {
          ++nFT0M;
        }
      }
      h.fill(HIST("hFT0M"), nFT0M);
      for (auto& p : pcFT0C) {
        if (isChargedParticle(p.pdgCode())) {
          ++nFT0C;
        }
      }
      h.fill(HIST("hFT0C"), nFT0C);
    }
  }

  PROCESS_SWITCH(Binner, calibrate, "Create binnings", false);

  void bin(aod::BCsWithTimestamps const& bcs, aod::McCollisions const& mccollisions, Particles const&)
  {
    auto bc = bcs.begin();
    multFT0M = ccdb->getForTimeStamp<TH1F>(path.value + "/hFT0M", bc.timestamp());
    multFT0C = ccdb->getForTimeStamp<TH1F>(path.value + "/hFT0C", bc.timestamp());
    if (multFT0C == nullptr && multFT0M == nullptr) {
      LOGP(fatal, "Unable to get the distributions from CCDB");
    }
    for (auto& mcc : mccollisions) {
      auto pcFT0M = pFT0M.sliceBy(perMcCol, mcc.globalIndex());
      auto pcFT0C = pFT0C.sliceBy(perMcCol, mcc.globalIndex());
      int nFT0M = 0;
      int nFT0C = 0;
      for (auto& p : pcFT0M) {
        if (isChargedParticle(p.pdgCode())) {
          ++nFT0M;
        }
      }
      h.fill(HIST("hFT0M"), nFT0M);
      for (auto& p : pcFT0C) {
        if (isChargedParticle(p.pdgCode())) {
          ++nFT0C;
        }
      }
      h.fill(HIST("hFT0C"), nFT0C);
      float percentileFT0M = 100.f * multFT0M->Integral(multFT0M->FindFixBin(nFT0C), multFT0M->FindLastBinAbove());
      float percentileFT0C = 100.f * multFT0C->Integral(multFT0C->FindFixBin(nFT0C), multFT0C->FindLastBinAbove());
      gencents(percentileFT0C, percentileFT0M);
    }
  }

  PROCESS_SWITCH(Binner, bin, "Bin collisions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<Binner>(cfgc)};
}
