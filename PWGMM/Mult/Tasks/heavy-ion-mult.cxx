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
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author Abhi Modak (contact: abhi.modak@cern.ch)
/// \help: To develop this code, I took help from the following codes and O2 analysis tutorial
// 1. https://github.com/AliceO2Group/O2Physics/blob/master/PWGMM/Mult/Tasks/dndeta.cxx
// 2. https://github.com/AliceO2Group/O2Physics/blob/master/PWGMM/Mult/Tasks/dndeta-hi.cxx
// 3. https://github.com/AliceO2Group/O2Physics/blob/master/PWGMM/Mult/Tasks/puremc-dndeta.cxx
// 4. O2 analysis tutorial: https://indico.cern.ch/event/1267433/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "bestCollisionTable.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype trackSelectionDCAXYonly =
  TrackSelectionFlags::kDCAxy;

AxisSpec axisEvent{4, -0.5, 3.5, "#Event"};
AxisSpec axisVtxZ{800, -20, 20, "Vertex Z"};
AxisSpec axisDCA = {601, -3.01, 3.01};
AxisSpec axisPT = {1000, -0.05, 49.95};
AxisSpec axisEta{200, -5, 5, "#eta"};
AxisSpec axisPhi{629, 0, 2 * M_PI, "#phi"};
AxisSpec axisMCEvent_ambiguity{6, -0.5, 5.5, "reco collisions per true collision"};

struct HeavyIonMultiplicity {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};
  Configurable<float> VtxRange{"vertex-range", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  ConfigurableAxis multHistBin{"MultDistBinning", {501, -0.5, 500.5}, ""};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin};
    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    if (doprocessData) {
      histos.add("MultHist", "MultHist", kTH1D, {axisMult}, true);
      histos.add("MultHist_Inelgt0", "MultHist_Inelgt0", kTH1D, {axisMult}, true);
      histos.add("EtaHist", "EtaHist", kTH1D, {axisEta}, true);
      histos.add("PhiHist", "PhiHist", kTH1D, {axisPhi}, true);
      histos.add("EtaVsVtxZHist", "EtaVsVtxZHist", kTH2D, {axisEta, axisVtxZ}, false);
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2D, {axisPhi, axisEta}, false);
      histos.add("DCAXYHist", "DCAXYHist", kTH1D, {axisDCA}, false);
      histos.add("DCAZHist", "DCAZHist", kTH1D, {axisDCA}, false);
      histos.add("pTHist", "pTHist", kTH1D, {axisPT}, true);
    }

    if (doprocessMC) {
      histos.add("MCEventHist_ambiguity", "MCEventHist_ambiguity", kTH1D, {axisMCEvent_ambiguity}, false);
      histos.add("DCAXYMCRecHist", "DCAXYMCRecHist", kTH1D, {axisDCA}, false);
      histos.add("DCAZMCRecHist", "DCAZMCRecHist", kTH1D, {axisDCA}, false);
      histos.add("pTMCRecHist", "pTMCRecHist", kTH1D, {axisPT}, true);
      histos.add("EtaVsVtxZMCRecHist", "EtaVsVtxZMCRecHist", kTH2D, {axisEta, axisVtxZ}, true);
      histos.add("MCRecEtaHist", "MCRecEtaHist", kTH1D, {axisEta}, true);
      histos.add("MCGenEtaHist", "MCGenEtaHist", kTH1D, {axisEta}, true);
      histos.add("MCRecPhiHist", "MCRecPhiHist", kTH1D, {axisPhi}, true);
      histos.add("MCGenPhiHist", "MCGenPhiHist", kTH1D, {axisPhi}, true);
      histos.add("MCRecPhiVsEtaHist", "MCRecPhiVsEtaHist", kTH2D, {axisPhi, axisEta}, false);
      histos.add("MCGenPhiVsEtaHist", "MCGenPhiVsEtaHist", kTH2D, {axisPhi, axisEta}, false);
      histos.add("MCRecMultHist", "MCRecMultHist", kTH1D, {axisMult}, true);
      histos.add("MCGenMultHist", "MCGenMultHist", kTH1D, {axisMult}, true);
      histos.add("MCGenVsRecMultHist", "MCGenVsRecMultHist", kTH2D, {axisMult, axisMult}, true);
      histos.add("MCRecMultHist_Inelgt0", "MCRecMultHist_Inelgt0", kTH1D, {axisMult}, true);
      histos.add("MCGenMultHist_Inelgt0", "MCGenMultHist_Inelgt0", kTH1D, {axisMult}, true);
      histos.add("MCGenVsRecMultHist_Inelgt0", "MCGenVsRecMultHist_Inelgt0", kTH2D, {axisMult, axisMult}, true);
    }
  }

  expressions::Filter trackSelectionProperMixed = ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionITS) &&
                                                  ifnode(ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true) &&
                                                  ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, trackSelectionDCAXYonly),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionDCA));

  void processData(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    auto NchTracks = 0;
    bool Inelgt0 = false;
    histos.fill(HIST("EventHist"), 0);
    if (collision.sel8()) {
      histos.fill(HIST("EventHist"), 1);
      if (std::abs(collision.posZ()) < VtxRange) {
        histos.fill(HIST("EventHist"), 2);
        histos.fill(HIST("VtxZHist"), collision.posZ());
        for (auto& track : tracks) {
          if (std::abs(track.eta()) < etaRange) {
            NchTracks++;
            histos.fill(HIST("EtaHist"), track.eta());
            histos.fill(HIST("PhiHist"), track.phi());
            histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
            histos.fill(HIST("EtaVsVtxZHist"), track.eta(), collision.posZ());
            histos.fill(HIST("DCAXYHist"), track.dcaXY());
            histos.fill(HIST("DCAZHist"), track.dcaZ());
            histos.fill(HIST("pTHist"), track.pt());
          }
        }
        histos.fill(HIST("MultHist"), NchTracks);
        if (NchTracks > 0) {
          Inelgt0 = true;
        }
        if (Inelgt0) {
          histos.fill(HIST("EventHist"), 3);
          histos.fill(HIST("MultHist_Inelgt0"), NchTracks);
        }
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processData, "process data", false);

  void processMC(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("MCEventHist_ambiguity"), RecCollisions.size());
    if (RecCollisions.size() == 0 || RecCollisions.size() > 1) {
      return;
    }
    auto NchGenTracks = 0;
    bool Inelgt0Gen = false;
    auto NchRecTracks = 0;
    bool Inelgt0Rec = false;
    for (auto& particle : GenParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (!particle.producedByGenerator()) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle == nullptr) {
        continue;
      }
      if (std::abs(pdgParticle->Charge()) >= 3) {
        if (std::abs(particle.eta()) < etaRange) {
          NchGenTracks++;
          histos.fill(HIST("MCGenEtaHist"), particle.eta());
          histos.fill(HIST("MCGenPhiHist"), particle.phi());
          histos.fill(HIST("MCGenPhiVsEtaHist"), particle.phi(), particle.eta());
        }
      }
    }
    histos.fill(HIST("MCGenMultHist"), NchGenTracks);
    if (NchGenTracks > 0) {
      Inelgt0Gen = true;
    }
    if (Inelgt0Gen) {
      histos.fill(HIST("MCGenMultHist_Inelgt0"), NchGenTracks);
    }
    for (auto& RecCollision : RecCollisions) {
      histos.fill(HIST("EventHist"), 0);
      if (RecCollision.sel8()) {
        histos.fill(HIST("EventHist"), 1);
        if (std::abs(RecCollision.posZ()) < VtxRange) {
          histos.fill(HIST("EventHist"), 2);
          histos.fill(HIST("VtxZHist"), RecCollision.posZ());

          auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
          for (auto& Rectrack : Rectrackspart) {
            if (std::abs(Rectrack.eta()) < etaRange) {
              NchRecTracks++;
              histos.fill(HIST("MCRecEtaHist"), Rectrack.eta());
              histos.fill(HIST("MCRecPhiHist"), Rectrack.phi());
              histos.fill(HIST("MCRecPhiVsEtaHist"), Rectrack.phi(), Rectrack.eta());
              histos.fill(HIST("EtaVsVtxZMCRecHist"), Rectrack.eta(), RecCollision.posZ());
              histos.fill(HIST("DCAXYMCRecHist"), Rectrack.dcaXY());
              histos.fill(HIST("DCAZMCRecHist"), Rectrack.dcaZ());
              histos.fill(HIST("pTMCRecHist"), Rectrack.pt());
            }
          }
          histos.fill(HIST("MCRecMultHist"), NchRecTracks);
          histos.fill(HIST("MCGenVsRecMultHist"), NchRecTracks, NchGenTracks);
          if (NchRecTracks > 0) {
            Inelgt0Rec = true;
            histos.fill(HIST("EventHist"), 3);
          }
          if (Inelgt0Rec) {
            histos.fill(HIST("MCRecMultHist_Inelgt0"), NchRecTracks);
          }
          if (Inelgt0Gen && Inelgt0Rec) {
            histos.fill(HIST("MCGenVsRecMultHist_Inelgt0"), NchRecTracks, NchGenTracks);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HeavyIonMultiplicity, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyIonMultiplicity>(cfgc)};
}
