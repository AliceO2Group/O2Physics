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

/// \file   flowPtEfficiency.cxx
/// \author Mingrui Zhao (mingrui.zhao@cern.ch), Zhiyong Lu (zhiyong.lu@cern.ch), Tao Jiang (tao.jiang@cern.ch)
/// \since  Jun/08/2023
/// \brief  a task to calculate the pt efficiency

#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowPtEfficiency {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 1000.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgTriggerkTVXinTRD, bool, false, "TRD triggered")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, false, "cut time intervals with dead ITS staves")
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, false, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, false, "V0A T0A 5 sigma cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, false, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 3000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelRun3ITSMatch, bool, false, "GlobalTrackRun3ITSMatching::Run3ITSall7Layers selection")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 70.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCcrossedrows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxyppPass3Enabled, bool, false, "switch of ppPass3 DCAxy pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAzPtDepEnabled, bool, false, "switch of DCAz pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgSelRunNumberEnabled, bool, false, "switch of run number selection")
  Configurable<std::vector<int>> cfgRunNumberList{"cfgRunNumberList", std::vector<int>{-1}, "runnumber list in consideration for analysis"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};

  // Filter the tracks
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  // Filter for collisions
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;

  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using MyMcParticles = soa::Filtered<aod::McParticles>;

  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using MyMcCollisions = soa::Filtered<aod::McCollisions>;

  // Additional filters for tracks
  TrackSelection myTrackSel;
  // Additional track selection cuts - Copy from flowGenericFramework.cxx
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  // Define the output
  HistogramRegistry registry{"registry"};

  bool isStable(int pdg)
  {
    if (std::abs(pdg) == 211)
      return true;
    if (std::abs(pdg) == 321)
      return true;
    if (std::abs(pdg) == 2212)
      return true;
    if (std::abs(pdg) == 11)
      return true;
    if (std::abs(pdg) == 13)
      return true;
    return false;
  }

  void init(InitContext const&)
  {
    const AxisSpec axisCounter{1, 0, +1, ""};
    // create histograms
    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});

    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});

    if (cfgTrkSelRun3ITSMatch) {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    } else {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    }
    if (cfgCutDCAxyppPass3Enabled) {
      myTrackSel.SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; });
    } else {
      myTrackSel.SetMaxDcaXY(cfgCutDCAxy);
    }
    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNCrossedRowsTPC(cfgCutTPCcrossedrows);
    if (!cfgCutDCAzPtDepEnabled)
      myTrackSel.SetMaxDcaZ(cfgCutDCAz);

    if (cfgEvSelMultCorrelation) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }
    if (cfgEvSelV0AT0ACut) {
      fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float centrality)
  {
    if (cfgTriggerkTVXinTRD && collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return 0;
    }
    if (cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }

    if (cfgEvSelOccupancy) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
        return 0;
    }

    if (cfgEvSelMultCorrelation) {
      auto multNTracksPV = collision.multNTracksPV();
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
    }

    // V0A T0A 5 sigma cut
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;

    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgCutDCAzPtDepEnabled && (track.dcaZ() > (0.004f + 0.013f / track.pt())))
      return false;
    return myTrackSel.IsSelected(track);
  }

  void processReco(MyCollisions::iterator const& collision, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::McParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    if (cfgSelRunNumberEnabled) {
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), runNumber))
        return;
    }
    const auto cent = collision.centFT0C();
    if (!eventSelected(collision, tracks.size(), cent))
      return;
    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCRec"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(FlowPtEfficiency, processReco, "process reconstructed information", true);

  void processSim(MyMcCollisions::iterator const& collision, aod::BCsWithTimestamps const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, MyMcParticles const& mcParticles)
  {
    if (cfgSelRunNumberEnabled) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      int runNumber = bc.runNumber();
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), runNumber))
        return;
    }
    if (collisions.size() > -1) {
      registry.fill(HIST("mcEventCounter"), 0.5);
      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.isPhysicalPrimary() && isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCGen"), mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(FlowPtEfficiency, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowPtEfficiency>(cfgc)};
}
