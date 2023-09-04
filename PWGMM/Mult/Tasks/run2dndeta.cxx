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
#include <cmath>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"

#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PseudorapidityDensity {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  ConfigurableAxis percentileBinning{"pBins",
                                     {VARIABLE_WIDTH, 0., 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100},
                                     "Centrality/multiplicity percentile binning"};

  Configurable<std::vector<float>> parRatio{"parRatio", {0.965475, -0.0010506, 0.00242833, 1.42873e-05, -1.49124e-05, 6.45095e-09, -1.0309e-08}, "parameters of pol6 of mes /reco tracks versus zvtx"};

  TF1* fRatio = 0;

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},        //
      {"TracksEtaZvtx_gt0", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},    //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},         //
      {"EventSelection", ";status;events", {HistType::kTH1F, {{4, 0.5, 4.5}}}}                                       //
    }                                                                                                                //
  };

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Selected INEL>0");
    x->SetBinLabel(4, "Rejected");

    if (doprocessGen) {
      registry.add({"EventsNtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}});
      registry.add({"EventsNtrkZvtxGen_t", "; N_{part}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_t", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_gt0", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_gt0t", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});

      registry.add({"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}});
      registry.add({"EventEfficiency", "; status; events", {HistType::kTH1F, {{5, 0.5, 5.5}}}});
      registry.add({"NotFoundEventZvtx", " ; Z_{vtx}", {HistType::kTH1F, {{201, -20.1, 20.1}}}});

      auto heff = registry.get<TH1>(HIST("EventEfficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");

      fRatio = new TF1("fRatio", "pol6", -15, 15);
      auto param = (std::vector<float>)parRatio;
      fRatio->SetParameters(param[0], param[1], param[2], param[3], param[4], param[5], param[6]);
    }

    if (doprocessBinned) {
      registry.add({"EventsNtrkZvtxBin", "; N_{trk}; Z_{vtx}; Percentile", {HistType::kTH3F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksEtaZvtxBin", "; #eta; Z_{vtx}; Percentile", {HistType::kTH3F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksPhiEtaBin", "; #varphi; #eta; Percentile", {HistType::kTH3F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}, percentileBinning}}});
      registry.add({"EventSelectionBin", "; status; Percentile; events", {HistType::kTH2F, {{3, 0.5, 3.5}, percentileBinning}}});

      hstat = registry.get<TH2>(HIST("EventSelectionBin"));
      x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");
    }
  }

  expressions::Filter trackTypeFilter = (aod::track::trackType == (uint8_t)aod::track::TrackTypeEnum::Run2Tracklet);

  using Trks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA>>;
  Partition<Trks> sample = nabs(aod::track::eta) < estimatorEta;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, Trks const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (doprocessGen) // we are sudying the MC, that needs correction
    {
      if (!useEvSel || (useEvSel && collision.sel7())) {

        auto z = collision.posZ();
        auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        double w = 1;
        if ((z > -12) && (z < 12)) {
          w = fRatio->Eval(z);
        }

        registry.fill(HIST("EventSelection"), 2., w);
        if (perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), 3., w);
        }
        registry.fill(HIST("EventsNtrkZvtx"), perCollisionSample.size(), z, w);
        for (auto& track : tracks) {
          registry.fill(HIST("TracksEtaZvtx"), track.eta(), z, w);
          registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta(), w);
          if (perCollisionSample.size() > 0) {
            registry.fill(HIST("TracksEtaZvtx_gt0"), track.eta(), z, w);
          }
        }
      } else {
        registry.fill(HIST("EventSelection"), 4.);
      }
    } else {
      if (!useEvSel || (useEvSel && collision.sel7())) {
        registry.fill(HIST("EventSelection"), 2.);
        auto z = collision.posZ();
        auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), 3.);
        }
        registry.fill(HIST("EventsNtrkZvtx"), perCollisionSample.size(), z);
        for (auto& track : tracks) {
          registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
          registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
          if (perCollisionSample.size() > 0) {
            registry.fill(HIST("TracksEtaZvtx_gt0"), track.eta(), z);
          }
        }
      } else {
        registry.fill(HIST("EventSelection"), 4.);
      }
    }
  }

  void processBinned(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA>> const& tracks)
  {
    auto p = collision.centRun2V0M();
    registry.fill(HIST("EventSelectionBin"), 1., p);
    if (!useEvSel || (useEvSel && collision.sel7())) {
      registry.fill(HIST("EventSelectionBin"), 2., p);
      auto z = collision.posZ();
      auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      registry.fill(HIST("EventsNtrkZvtxBin"), perCollisionSample.size(), z, p);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtxBin"), track.eta(), z, p);
        registry.fill(HIST("TracksPhiEtaBin"), track.phi(), track.eta(), p);
      }
    } else {
      registry.fill(HIST("EventSelectionBin"), 4., p);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity, processBinned, "Process centrality/mult. percentile binned", false);

  using Particles = aod::McParticles;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<soa::Filtered<Particles>> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;

  void processGen(aod::McCollisions::iterator const& mcCollision, soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, soa::Filtered<Particles> const& particles, soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA>> const& /*tracks*/)
  {
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;
    for (auto& particle : perCollisionMCSample) {
      auto charge = 0;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (charge != 0) {
        nCharged++;
      }
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());
    registry.fill(HIST("EventEfficiency"), 1.);

    if (nCharged > 0) {
      registry.fill(HIST("EventEfficiency"), 2.);
    }
    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 3.);
      if (!useEvSel || (useEvSel && collision.sel7())) {
        atLeastOne = true;
        auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        registry.fill(HIST("EventEfficiency"), 4.);
        if (perCollisionSample.size() > 0) {
          atLeastOne_gt0 = true;
          registry.fill(HIST("EventEfficiency"), 5.);
        }
        registry.fill(HIST("EventsNtrkZvtxGen"), perCollisionSample.size(), collision.posZ());
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (charge != 0) {
        registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
        if (perCollisionMCSample.size() > 0) {
          registry.fill(HIST("TracksEtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ());
        }
        if (atLeastOne) {
          registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), mcCollision.posZ());
          if (atLeastOne_gt0) {
            registry.fill(HIST("TracksEtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ());
          }
        }
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity, processGen, "Process generator-level info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensity>(cfgc)};
}
