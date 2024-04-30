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
#include "Common/DataModel/Multiplicity.h"
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
using namespace o2::aod::evsel;

using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionDataTableCorrelation = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using CollisionDataTableCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
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

AxisSpec axisEvent{5, -0.5, 4.5, "#Event"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z"};
AxisSpec axisDCA = {601, -3.01, 3.01};
AxisSpec axisEta{40, -2, 2, "#eta"};
AxisSpec axisPhi{629, 0, 2 * M_PI, "#phi"};
AxisSpec axisPhi2{{0, M_PI / 2, M_PI, M_PI * 3. / 2, 2 * M_PI}, "#phi"};
AxisSpec axisMCEvent_ambiguity{10, -0.5, 9.5, "reco collisions per true collision"};
AxisSpec axisCent{100, 0, 100, "#Cent"};

struct HeavyIonMultiplicity {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};
  Configurable<float> VtxRange{"vertex-range", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  ConfigurableAxis multHistBin{"MultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis PVHistBin{"PVDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FV0multHistBin{"FV0MultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FT0multHistBin{"FT0MultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FT0AmultHistBin{"FT0AMultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FT0CmultHistBin{"FT0CMultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis pTHistBin{"pTHistBin", {VARIABLE_WIDTH, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100}, ""};
  ConfigurableAxis CentralityBinning{"CentralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin};
    AxisSpec axisPV = {PVHistBin};
    AxisSpec axisFV0Mult = {FV0multHistBin};
    AxisSpec axisFT0Mult = {FT0multHistBin};
    AxisSpec axisFT0AMult = {FT0AmultHistBin};
    AxisSpec axisFT0CMult = {FT0CmultHistBin};
    AxisSpec CentAxis = {CentralityBinning, "Centrality"};
    AxisSpec axisPT = {pTHistBin};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoITSROFrameBorder");
    x->SetBinLabel(4, "kNoTimeFrameBorder");
    x->SetBinLabel(5, "|VtxZ|<10");

    if (doprocessData) {
      histos.add("MultHist", "MultHist", kTH1D, {axisMult}, true);
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
      histos.add("MCRecEtaHist", "MCRecEtaHist", kTH1D, {axisEta}, true);
      histos.add("MCRecPhiHist", "MCRecPhiHist", kTH1D, {axisPhi}, true);
      histos.add("MCRecPhiVsEtaHist", "MCRecPhiVsEtaHist", kTH2D, {axisPhi, axisEta}, false);
      histos.add("EtaVsVtxZMCRecHist", "EtaVsVtxZMCRecHist", kTH2D, {axisEta, axisVtxZ}, true);
      histos.add("DCAXYMCRecHist", "DCAXYMCRecHist", kTH1D, {axisDCA}, false);
      histos.add("DCAZMCRecHist", "DCAZMCRecHist", kTH1D, {axisDCA}, false);
      histos.add("pTMCRecHist", "pTMCRecHist", kTH1D, {axisPT}, true);
      histos.add("MCRecMultHist", "MCRecMultHist", kTH1D, {axisMult}, true);

      histos.add("MCGenEtaHist", "MCGenEtaHist", kTH1D, {axisEta}, true);
      histos.add("MCGenPhiHist", "MCGenPhiHist", kTH1D, {axisPhi}, true);
      histos.add("MCGenPhiVsEtaHist", "MCGenPhiVsEtaHist", kTH2D, {axisPhi, axisEta}, false);
      histos.add("EtaVsVtxZMCGenHist", "EtaVsVtxZMCGenHist", kTH2D, {axisEta, axisVtxZ}, true);
      histos.add("MCGenMultHist", "MCGenMultHist", kTH1D, {axisMult}, true);
      histos.add("MCGenVsRecMultHist", "MCGenVsRecMultHist", kTH2D, {axisMult, axisMult}, true);
      histos.add("pTMCGenHist", "pTMCGenHist", kTH1D, {axisPT}, true);
      histos.add("VtxZGenHist", "VtxZGenHist", kTH1D, {axisVtxZ}, false);

      histos.add("MCGenEtaHistAll", "MCGenEtaHistAll", kTH1D, {axisEta}, true);
      histos.add("MCGenPhiHistAll", "MCGenPhiHistAll", kTH1D, {axisPhi}, true);
      histos.add("MCGenPhiVsEtaHistAll", "MCGenPhiVsEtaHistAll", kTH2D, {axisPhi, axisEta}, false);
      histos.add("EtaVsVtxZMCGenHistAll", "EtaVsVtxZMCGenHistAll", kTH2D, {axisEta, axisVtxZ}, true);
      histos.add("MCGenMultHistAll", "MCGenMultHistAll", kTH1D, {axisMult}, true);
    }

    if (doprocessDataCentFT0C) {
      histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
      histos.add("CentHistInsideTrackloop", "CentHistInsideTrackloop", kTH1D, {axisCent}, false);
      histos.add("hdatamult", "hdatamult", kTHnSparseF, {axisVtxZ, axisMult, CentAxis}, false);
      histos.add("hdatadndeta", "hdatadndeta", kTHnSparseF, {axisVtxZ, CentAxis, axisEta, axisPhi2}, false);
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, CentAxis}, false);
    }

    if (doprocessMCCentFT0C) {
      histos.add("CentPercentileMCRecHist", "CentPercentileMCRecHist", kTH1D, {axisCent}, false);
      histos.add("hmczvtxcent", "hmczvtxcent", kTH2D, {axisVtxZ, CentAxis}, false);
      histos.add("hrecdndeta_pion", "hrecdndeta_pion", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hrecdndeta_kaon", "hrecdndeta_kaon", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hrecdndeta_proton", "hrecdndeta_proton", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hrecdndeta_other", "hrecdndeta_other", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hrecdndeta_bkg", "hrecdndeta_bkg", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hrecdndeta_NotPrimary", "hrecdndeta_NotPrimary", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hgendndeta_pion", "hgendndeta_pion", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hgendndeta_kaon", "hgendndeta_kaon", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hgendndeta_proton", "hgendndeta_proton", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
      histos.add("hgendndeta_other", "hgendndeta_other", kTHnSparseF, {axisVtxZ, axisEta, axisPhi2, CentAxis, axisPT}, false);
    }

    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisFT0AMult, axisMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisFT0CMult, axisMult}, true);
      histos.add("GlobalMult_vs_FT0", "GlobalMult_vs_FT0", kTH2F, {axisFT0Mult, axisMult}, true);
      histos.add("GlobalMult_vs_FV0", "GlobalMult_vs_FV0", kTH2F, {axisFV0Mult, axisMult}, true);
      histos.add("GlobalMult_vs_NumPVContributor", "GlobalMult_vs_NumPVContributor", kTH2F, {axisPV, axisMult}, true);
    }
  }

  expressions::Filter trackSelectionProperMixed = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionITS) &&
                                                  ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true) &&
                                                  ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, trackSelectionDCAXYonly),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionDCA));

  void processData(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    auto NchTracks = 0;
    histos.fill(HIST("EventHist"), 0);
    if (collision.sel8()) {
      histos.fill(HIST("EventHist"), 1);
      if (collision.selection_bit(kNoITSROFrameBorder)) {
        histos.fill(HIST("EventHist"), 2);
        if (collision.selection_bit(kNoTimeFrameBorder)) {
          histos.fill(HIST("EventHist"), 3);
          if (std::abs(collision.posZ()) < VtxRange) {
            histos.fill(HIST("EventHist"), 4);
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
          }
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

    auto NchRecTracks = 0;
    auto NchGenTracks = 0;
    for (auto& RecCollision : RecCollisions) {
      histos.fill(HIST("EventHist"), 0);
      if (RecCollision.sel8()) {
        histos.fill(HIST("EventHist"), 1);
        if (RecCollision.selection_bit(kNoITSROFrameBorder)) {
          histos.fill(HIST("EventHist"), 2);
          if (RecCollision.selection_bit(kNoTimeFrameBorder)) {
            histos.fill(HIST("EventHist"), 3);
            if (std::abs(RecCollision.posZ()) < VtxRange) {
              histos.fill(HIST("EventHist"), 4);
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
                    histos.fill(HIST("EtaVsVtxZMCGenHist"), particle.eta(), RecCollision.posZ());
                    histos.fill(HIST("pTMCGenHist"), particle.pt());
                  }
                }
              }

              histos.fill(HIST("MCRecMultHist"), NchRecTracks);
              histos.fill(HIST("MCGenMultHist"), NchGenTracks);
              histos.fill(HIST("MCGenVsRecMultHist"), NchRecTracks, NchGenTracks);
            }
          }
        }
      }
    }

    auto NchGenTracksAll = 0;
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
        if (std::abs(TrueCollision.posZ()) < VtxRange) {
          histos.fill(HIST("VtxZGenHist"), TrueCollision.posZ());
          if (std::abs(particle.eta()) < etaRange) {
            NchGenTracksAll++;
            histos.fill(HIST("MCGenEtaHistAll"), particle.eta());
            histos.fill(HIST("MCGenPhiHistAll"), particle.phi());
            histos.fill(HIST("MCGenPhiVsEtaHistAll"), particle.phi(), particle.eta());
            histos.fill(HIST("EtaVsVtxZMCGenHistAll"), particle.eta(), TrueCollision.posZ());
          }
        }
      }
    }
    histos.fill(HIST("MCGenMultHistAll"), NchGenTracksAll);
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processMC, "process MC", false);

  void processDataCentFT0C(CollisionDataTableCentFT0C::iterator const& collision, FilTrackDataTable const& tracks)
  {
    float cent = -1;
    auto NchTracks = 0;
    constexpr auto hasCentrality = CollisionDataTableCentFT0C::template contains<aod::CentFT0Cs>();
    histos.fill(HIST("EventHist"), 0);
    if (collision.sel8()) {
      histos.fill(HIST("EventHist"), 1);
      if (collision.selection_bit(kNoITSROFrameBorder)) {
        histos.fill(HIST("EventHist"), 2);
        if (collision.selection_bit(kNoTimeFrameBorder)) {
          histos.fill(HIST("EventHist"), 3);
          if constexpr (hasCentrality) {
            cent = collision.centFT0C();
            if (cent < 0.0 || cent > 100.0) {
              return;
            }
            histos.fill(HIST("EventHist"), 4);
            histos.fill(HIST("VtxZHist"), collision.posZ());
            histos.fill(HIST("CentPercentileHist"), cent);
            histos.fill(HIST("hdatazvtxcent"), collision.posZ(), cent);
            for (auto& track : tracks) {
              if (std::abs(track.eta()) < etaRange) {
                NchTracks++;
                histos.fill(HIST("CentHistInsideTrackloop"), cent);
                histos.fill(HIST("hdatadndeta"), collision.posZ(), cent, track.eta(), track.phi());
              }
            }
            histos.fill(HIST("hdatamult"), collision.posZ(), NchTracks, cent);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processDataCentFT0C, "process data CentFT0C", false);

  void processMCCentFT0C(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    if (RecCollisions.size() == 0 || RecCollisions.size() > 1) {
      return;
    }

    float cent = -1;
    constexpr auto hasCentrality = CollisionMCRecTableCentFT0C::template contains<aod::CentFT0Cs>();

    for (auto& RecCollision : RecCollisions) {
      histos.fill(HIST("EventHist"), 0);
      if (RecCollision.sel8()) {
        histos.fill(HIST("EventHist"), 1);
        if (RecCollision.selection_bit(kNoITSROFrameBorder)) {
          histos.fill(HIST("EventHist"), 2);
          if (RecCollision.selection_bit(kNoTimeFrameBorder)) {
            histos.fill(HIST("EventHist"), 3);
            if constexpr (hasCentrality) {
              cent = RecCollision.centFT0C();
              if (cent < 0.0 || cent > 100.0) {
                continue;
              }
              histos.fill(HIST("EventHist"), 4);
              histos.fill(HIST("VtxZHist"), RecCollision.posZ());
              histos.fill(HIST("CentPercentileMCRecHist"), cent);
              histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), cent);

              auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
              for (auto& Rectrack : Rectrackspart) {
                if (std::abs(Rectrack.eta()) > etaRange) {
                  continue;
                }
                if (Rectrack.has_mcParticle()) {
                  auto mcpart = Rectrack.mcParticle();
                  if (std::abs(mcpart.eta()) > etaRange) {
                    continue;
                  }
                  if (mcpart.isPhysicalPrimary()) {
                    if (std::abs(mcpart.pdgCode()) == 211) {
                      histos.fill(HIST("hrecdndeta_pion"), RecCollision.posZ(), mcpart.eta(), mcpart.phi(), cent, mcpart.pt());
                    } else if (std::abs(mcpart.pdgCode()) == 321) {
                      histos.fill(HIST("hrecdndeta_kaon"), RecCollision.posZ(), mcpart.eta(), mcpart.phi(), cent, mcpart.pt());
                    } else if (std::abs(mcpart.pdgCode()) == 2212) {
                      histos.fill(HIST("hrecdndeta_proton"), RecCollision.posZ(), mcpart.eta(), mcpart.phi(), cent, mcpart.pt());
                    } else {
                      histos.fill(HIST("hrecdndeta_other"), RecCollision.posZ(), mcpart.eta(), mcpart.phi(), cent, mcpart.pt());
                    }
                  } else {
                    histos.fill(HIST("hrecdndeta_NotPrimary"), RecCollision.posZ(), mcpart.eta(), mcpart.phi(), cent, mcpart.pt());
                  }
                } else {
                  histos.fill(HIST("hrecdndeta_bkg"), RecCollision.posZ(), Rectrack.eta(), Rectrack.phi(), cent, Rectrack.pt());
                }
              }

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
                  if (std::abs(particle.eta()) > etaRange) {
                    continue;
                  }
                  if (std::abs(particle.pdgCode()) == 211) {
                    histos.fill(HIST("hgendndeta_pion"), RecCollision.posZ(), particle.eta(), particle.phi(), cent, particle.pt());
                  } else if (std::abs(particle.pdgCode()) == 321) {
                    histos.fill(HIST("hgendndeta_kaon"), RecCollision.posZ(), particle.eta(), particle.phi(), cent, particle.pt());
                  } else if (std::abs(particle.pdgCode()) == 2212) {
                    histos.fill(HIST("hgendndeta_proton"), RecCollision.posZ(), particle.eta(), particle.phi(), cent, particle.pt());
                  } else {
                    histos.fill(HIST("hgendndeta_other"), RecCollision.posZ(), particle.eta(), particle.phi(), cent, particle.pt());
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processMCCentFT0C, "process MC CentFT0C", false);

  void processCorrelation(CollisionDataTableCorrelation::iterator const& collision, FilTrackDataTable const& tracks)
  {
    auto NchTracks = 0;
    histos.fill(HIST("EventHist"), 0);
    if (collision.sel8()) {
      histos.fill(HIST("EventHist"), 1);
      if (collision.selection_bit(kNoITSROFrameBorder)) {
        histos.fill(HIST("EventHist"), 2);
        if (collision.selection_bit(kNoTimeFrameBorder)) {
          histos.fill(HIST("EventHist"), 3);
          if (std::abs(collision.posZ()) < VtxRange) {
            histos.fill(HIST("EventHist"), 4);
            histos.fill(HIST("VtxZHist"), collision.posZ());
            for (auto& track : tracks) {
              if (std::abs(track.eta()) < etaRange) {
                NchTracks++;
              }
            }
            histos.fill(HIST("GlobalMult_vs_FT0A"), collision.multFT0A(), NchTracks);
            histos.fill(HIST("GlobalMult_vs_FT0C"), collision.multFT0C(), NchTracks);
            histos.fill(HIST("GlobalMult_vs_FT0"), collision.multFT0A() + collision.multFT0C(), NchTracks);
            histos.fill(HIST("GlobalMult_vs_FV0"), collision.multFV0A(), NchTracks);
            histos.fill(HIST("GlobalMult_vs_NumPVContributor"), collision.numContrib(), NchTracks);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processCorrelation, "do correlation between FT0/FV0 Mult vs GlobalMult", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyIonMultiplicity>(cfgc)};
}
