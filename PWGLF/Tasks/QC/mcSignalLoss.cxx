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
//
// Compute signal loss for pions, protons, kaons, k0s, lambdas

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
//
// Build hypertriton candidates from V0s and tracks

#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

namespace
{
static const std::vector<int> pdgCodes{211, 321, 2212, 1000010020, 1000020030, 310, 3122, 1010010030};
} // namespace

struct mcsignalloss {
  ConfigurableAxis ptAxis{"ptAxis", {80., 0.f, 10.f}, "pt axis binning"};
  ConfigurableAxis decLAxis{"decLAxis", {100., 0.f, 50.f}, "decay length binning"};
  ConfigurableAxis zAxis{"zAxis", {100, -20, 20}, "z axis binning"};
  ConfigurableAxis pdgAxis{"pdgAxis", {static_cast<double>(pdgCodes.size()), 0, static_cast<double>(pdgCodes.size())}, "pdg axis binning"};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<bool> isRecoCollision;
  Service<o2::framework::O2DatabasePDG> pdgDB;
  SliceCache cache;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;

  float calcDecL(const aod::McParticle& mcPart, aod::McParticles const&)
  {
    // look for the pion daughter and compute the decay length
    float decL = -1;
    if (mcPart.has_daughters()) {
      for (auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) != 211) {
          continue;
        }
        if (daughter.getProcess() == 4) { // TMCProcess::kPDecay
          auto posPrimVtx = array{mcPart.vx(), mcPart.vy(), mcPart.vz()};
          auto secVtx = array{daughter.vx(), daughter.vy(), daughter.vz()};
          decL = std::hypot(secVtx[0] - posPrimVtx[0], secVtx[1] - posPrimVtx[1], secVtx[2] - posPrimVtx[2]);
          return decL;
        }
      }
    }
    return decL;
  }

  void init(o2::framework::InitContext&)
  {
    histos.add<TH1>("mcGenCollisionVtx", ";#it{z} (cm)", HistType::kTH1F, {zAxis});
    histos.add<TH1>("recCollisionVtx", ";#it{z} (cm)", HistType::kTH1F, {zAxis});
    histos.add<TH2>("mcGenAll", ";#it{p_{T}(GeV/#it{c}); pdg", HistType::kTH2F, {ptAxis, pdgAxis});
    histos.add<TH2>("mcGenAllRecoColl", ";#it{p_{T}(GeV/#it{c}); pdg", HistType::kTH2F, {ptAxis, pdgAxis});
    histos.add<TH3>("mcGenV0s", ";#it{p_{T}(GeV/#it{c}); Decay Length (cm); pdg", HistType::kTH3F, {ptAxis, decLAxis, pdgAxis});
    histos.add<TH3>("mcGenV0sRecoColl", ";#it{p_{T}(GeV/#it{c}); Decay Length (cm); pdg", HistType::kTH3F, {ptAxis, decLAxis, pdgAxis});
  };

  void process(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC)
  {
    isRecoCollision.clear();
    isRecoCollision.resize(mcCollisions.size(), false);

    for (const auto& collision : collisions) {

      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > 10)
        continue;

      histos.fill(HIST("recCollisionVtx"), collision.posZ());

      if (collision.has_mcCollision()) {
        isRecoCollision[collision.mcCollisionId()] = true;
      }
    }

    for (const auto& mcCollision : mcCollisions) {
      auto mcParticles_per_coll = particlesMC.sliceBy(perMCCollision, mcCollision.globalIndex());
      if (!pwglf::isINELgtNmc(mcParticles_per_coll, 0, pdgDB)) {
        continue;
      }

      histos.fill(HIST("mcGenCollisionVtx"), mcCollision.posZ());
      for (const auto& mcPart : mcParticles_per_coll) {
        if (!mcPart.isPhysicalPrimary() || mcPart.mcCollisionId() < 0) {
          continue;
        }
        for (unsigned int i = 0; i < pdgCodes.size(); i++) {
          if (std::abs(mcPart.pdgCode()) == pdgCodes[i]) {
            bool fillV0s = i > 4;
            histos.fill(HIST("mcGenAll"), mcPart.pt(), i);
            if (fillV0s) {
              histos.fill(HIST("mcGenV0s"), mcPart.pt(), calcDecL(mcPart, particlesMC), i);
            }
            if (isRecoCollision[mcPart.mcCollisionId()]) {
              histos.fill(HIST("mcGenAllRecoColl"), mcPart.pt(), i);
              if (fillV0s) {
                histos.fill(HIST("mcGenV0sRecoColl"), mcPart.pt(), calcDecL(mcPart, particlesMC), i);
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcsignalloss>(cfgc)};
}
