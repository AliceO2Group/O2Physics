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

/// \file   flowMcUpc.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch), Yongxi Du (yongxi.du@cern.ch)
/// \since  Apr/2/2026
/// \brief  flow efficiency analysis on UPC MC

#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowMcUpc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> minB{"minB", 0.0f, "min impact parameter"};
  Configurable<float> maxB{"maxB", 20.0f, "max impact parameter"};
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.1f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 1000.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgDcaxy, bool, true, "choose dcaxy")

  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, ""};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};
  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  double epsilon = 1e-6;

  // using McParts = soa::Join<aod::UDMcParticles, aod::UDMcTrackLabels>;

  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    const AxisSpec axisVertex{20, -10, 10, "Vtxz (cm)"};
    const AxisSpec axisEta{20, -1., 1., "#eta"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    // QA histograms
    histos.add<TH1>("mcEventCounter", "Monte Carlo Truth EventCounter", HistType::kTH1F, {{5, 0, 5}});
    histos.add<TH1>("RecoProcessEventCounter", "Reconstruction EventCounter", HistType::kTH1F, {{5, 0, 5}});
    histos.add<TH1>("hImpactParameter", "hImpactParameter", HistType::kTH1D, {axisB});
    histos.add<TH1>("RecoProcessTrackCounter", "Reconstruction TrackCounter", HistType::kTH1F, {{5, 0, 5}});
    histos.add<TH1>("numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{100, -0.5f, 99.5f}});

    histos.add<TH1>("hPtMCGen", "Monte Carlo Truth; pT (GeV/c);", {HistType::kTH1D, {axisPt}});
    histos.add<TH3>("hEtaPtVtxzMCGen", "Monte Carlo Truth; #eta; p_{T} (GeV/c); V_{z} (cm);", {HistType::kTH3D, {axisEta, axisPt, axisVertex}});
    histos.add<TH1>("hPtReco", "Monte Carlo Reco Global; pT (GeV/c);", {HistType::kTH1D, {axisPt}});
    histos.add<TH3>("hEtaPtVtxzMCReco", "Monte Carlo Global; #eta; p_{T} (GeV/c); V_{z} (cm);", {HistType::kTH3D, {axisEta, axisPt, axisVertex}});
  }

  // template <typename TCollision>
  // bool eventSelected(TCollision collision)
  // {
  //   return true;
  // }

  template <typename TTrack>
  bool trackSelected(TTrack const& track)
  {
    if (!track.hasTPC())
      return false;
    if (!track.isPVContributor())
      return false;
    // auto momentum = std::array<double, 3>{track.px(), track.py(), track.pz()};
    auto pt = track.pt();
    if (pt < cfgPtCutMin || pt > cfgPtCutMax) {
      return false;
    }
    double dcaLimit = 0.0105 + 0.035 / std::pow(pt, 1.1);
    if (cfgDcaxy && !(std::fabs(track.dcaXY()) < dcaLimit)) {
      return false;
    }
    return true;
  }

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;

  void processMCTrue(aod::UDMcCollisions const& mcCollisions, aod::UDMcParticles const& mcParts, aod::BCs const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    for (const auto& mcCollision : mcCollisions) {
      histos.fill(HIST("mcEventCounter"), 0.5);
      float imp = mcCollision.impactParameter();
      float vtxz = mcCollision.posZ();

      if (imp >= minB && imp <= maxB) {
        // event within range
        histos.fill(HIST("hImpactParameter"), imp);

        auto const& tempParts = mcParts.sliceBy(partPerMcCollision, static_cast<int64_t>(mcCollision.globalIndex()));

        for (auto const& mcParticle : tempParts) {
          auto momentum = std::array<double, 3>{mcParticle.px(), mcParticle.py(), mcParticle.pz()};
          int pdgCode = std::abs(mcParticle.pdgCode());

          double pt = RecoDecay::pt(momentum);
          double eta = RecoDecay::eta(momentum);

          if (pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != PDG_t::kKPlus && pdgCode != PDG_t::kProton)
            continue;

          if (!mcParticle.isPhysicalPrimary())
            continue;
          if (std::fabs(eta) > cfgCutEta) // main acceptance
            continue;

          histos.fill(HIST("hPtMCGen"), pt);
          histos.fill(HIST("hEtaPtVtxzMCGen"), eta, pt, vtxz);
        }
      }
    }
  }
  PROCESS_SWITCH(FlowMcUpc, processMCTrue, "process pure simulation information", true);

  using MCRecoTracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;
  using MCRecoCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>;

  // PresliceUnsorted<MCRecoTracks> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;
  Preslice<MCRecoTracks> trackPerCollision = aod::udtrack::udCollisionId; // sorted preslice used because the pair track-collision is already sorted in processDataSG function

  void processReco(MCRecoCollisions const& collisions, MCRecoTracks const& tracks, aod::UDMcParticles const& mcParticles)
  {
    histos.fill(HIST("numberOfRecoCollisions"), mcParticles.size()); // number of times coll was reco-ed
    // std::cout << "process reco" << std::endl;
    for (const auto& collision : collisions) {
      Partition<MCRecoTracks> pvContributors = aod::udtrack::isPVContributor == true;
      pvContributors.bindTable(tracks);
      // std::cout << "collision loop" << std::endl;
      histos.fill(HIST("RecoProcessEventCounter"), 0.5);
      // if (!eventSelected(collision))
      //   return;
      histos.fill(HIST("RecoProcessEventCounter"), 1.5);
      // if (!collision.has_udMcCollision())
      //   return;
      histos.fill(HIST("RecoProcessEventCounter"), 2.5);
      histos.fill(HIST("RecoProcessEventCounter"), 3.5);

      float vtxz = collision.posZ();

      // auto const& tempTracks = tracks.sliceBy(trackPerCollision, static_cast<int64_t>(collision.globalIndex()));
      // std::cout << "sliced" << std::endl;

      for (const auto& track : tracks) {
        histos.fill(HIST("RecoProcessTrackCounter"), 0.5);
        // std::cout << "track loop" << std::endl;
        // focus on bulk: e, mu, pi, k, p
        auto momentum = std::array<double, 3>{track.px(), track.py(), track.pz()};
        double pt = RecoDecay::pt(momentum);
        double eta = RecoDecay::eta(momentum);
        // double phi = RecoDecay::phi(momentum);
        if (!trackSelected(track) || (!track.has_udMcParticle()))
          continue;
        histos.fill(HIST("RecoProcessTrackCounter"), 1.5);
        // std::cout << "track selected" << std::endl;
        auto mcParticle = track.udMcParticle();
        // std::cout << "mc particle" << std::endl;
        int pdgCode = std::abs(mcParticle.pdgCode());
        // std::cout <<  "pdg code" << std::endl;

        // double pt = recoMC.Pt();
        // double eta = recoMC.Eta();
        if (pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != PDG_t::kKPlus && pdgCode != PDG_t::kProton) {
          // std::cout << "pdg code not in list" << std::endl;
          continue;
        }
        histos.fill(HIST("RecoProcessTrackCounter"), 2.5);
        if (std::fabs(eta) > cfgCutEta) {
          // std::cout << "cfgcuteta" << std::endl;
          continue;
        } // main acceptance
        histos.fill(HIST("RecoProcessTrackCounter"), 3.5);
        if (!mcParticle.isPhysicalPrimary()) {
          // std::cout << "not physical primary" << std::endl;
          continue;
        }
        histos.fill(HIST("RecoProcessTrackCounter"), 4.5);

        histos.fill(HIST("hPtReco"), pt);
        histos.fill(HIST("hEtaPtVtxzMCReco"), eta, pt, vtxz);
        // std::cout << "first loop end" << std::endl;
      }
    }
  }
  PROCESS_SWITCH(FlowMcUpc, processReco, "process reconstructed information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowMcUpc>(cfgc)};
}
