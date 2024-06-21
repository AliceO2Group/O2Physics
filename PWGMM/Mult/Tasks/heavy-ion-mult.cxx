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

using CollisionDataTable = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::Mults>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::Mults>>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

enum {
  kTrackTypebegin = 0,
  kGlobalplusITS = 1,
  kGlobalonly,
  kITSonly,
  kTrackTypeend
};

enum {
  kGenpTbegin = 0,
  kNoGenpTVar = 1,
  kGenpTup,
  kGenpTdown,
  kGenpTend
};

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

AxisSpec axisEvent{9, 0.5, 9.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisPhi{{0, M_PI / 4, M_PI / 2, M_PI * 3. / 4, M_PI, M_PI * 5. / 4, M_PI * 3. / 2, M_PI * 7. / 4, 2 * M_PI}, "#phi", "PhiAxis"};
AxisSpec axisPhi2{629, 0, 2 * M_PI, "#phi"};
AxisSpec axisCent{100, 0, 100, "#Cent"};
AxisSpec AxisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
AxisSpec AxisGenpTVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
AxisSpec axisTracks{9, 0.5, 9.5, "#tracks", "TrackAxis"};

struct HeavyIonMultiplicity {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};
  Configurable<float> VtxRange{"vertex-range", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  ConfigurableAxis multHistBin{"MultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis PVHistBin{"PVDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FV0AmultHistBin{"FV0AMultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FT0AmultHistBin{"FT0AMultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis FT0CmultHistBin{"FT0CMultDistBinning", {501, -0.5, 500.5}, ""};
  ConfigurableAxis pTHistBin{"pTHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis CentralityBinning{"CentralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};

  Configurable<bool> IsApplySameBunchPileup{"IsApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> IsApplyGoodZvtxFT0vsPV{"IsApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> IsApplyVertexITSTPC{"IsApplyVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<bool> IsApplyVertexTOFmatched{"IsApplyVertexTOFmatched", true, "Enable VertexTOFmatched cut"};
  Configurable<bool> IsApplyVertexTRDmatched{"IsApplyVertexTRDmatched", true, "Enable VertexTRDmatched cut"};
  Configurable<bool> IsApplyExtraCorrCut{"IsApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> IsApplyExtraPhiCut{"IsApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<float> NPVtracksCut{"NPVtracksCut", 1.0f, "Apply extra NPVtracks cut"};
  Configurable<float> FT0CCut{"FT0CCut", 1.0f, "Apply extra FT0C cut"};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin};
    AxisSpec axisPV = {PVHistBin};
    AxisSpec axisFV0AMult = {FV0AmultHistBin};
    AxisSpec axisFT0AMult = {FT0AmultHistBin};
    AxisSpec axisFT0CMult = {FT0CmultHistBin};
    AxisSpec CentAxis = {CentralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPT = {pTHistBin};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoSameBunchPileup");  // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV");  // small difference between z-vertex from PV and from FT0
    x->SetBinLabel(5, "kIsVertexITSTPC");     // at least one ITS-TPC track (reject vertices built from ITS-only tracks)
    x->SetBinLabel(6, "kIsVertexTOFmatched"); // at least one of vertex contributors is matched to TOF
    x->SetBinLabel(7, "kIsVertexTRDmatched"); // at least one of vertex contributors is matched to TRD
    x->SetBinLabel(8, "Centrality");
    x->SetBinLabel(9, "ApplyExtraCorrCut");

    if (doprocessData) {
      histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
      histos.add("CentHistInsideTrackloop", "CentHistInsideTrackloop", kTH1D, {axisCent}, false);
      histos.add("hdatamult", "hdatamult", kTHnSparseD, {axisVtxZ, axisMult, CentAxis}, false);
      histos.add("hdatadndeta", "hdatadndeta", kTHnSparseD, {axisVtxZ, CentAxis, axisEta, axisPhi, axisPT, AxisTrackType}, false);
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, CentAxis}, false);
      histos.add("PhiVsEtaHistBeforePhiCut", "PhiVsEtaHistBeforePhiCut", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("PhiVsEtaHistAfterPhiCut", "PhiVsEtaHistAfterPhiCut", kTH2D, {axisPhi2, axisEta}, false);
    }

    if (doprocessMonteCarlo || doprocessMCpTefficiency || doprocessMCcheckFakeTracks) {
      histos.add("CentPercentileMCRecHist", "CentPercentileMCRecHist", kTH1D, {axisCent}, false);
      histos.add("hmczvtxcent", "hmczvtxcent", kTH2D, {axisVtxZ, CentAxis}, false);
    }

    if (doprocessMonteCarlo) {
      histos.add("MCCentHistInsideTrackloop", "MCCentHistInsideTrackloop", kTH1D, {axisCent}, false);
      histos.add("MCrecPhiVsEtaHistBeforePhiCut", "MCrecPhiVsEtaHistBeforePhiCut", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("MCrecPhiVsEtaHistAfterPhiCut", "MCrecPhiVsEtaHistAfterPhiCut", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hmcrecdndeta", "hmcrecdndeta", kTHnSparseD, {axisVtxZ, CentAxis, axisEta, axisPhi}, false);
      histos.add("hmcgendndeta", "hmcgendndeta", kTHnSparseD, {axisVtxZ, CentAxis, axisEta, axisPhi, AxisGenpTVary}, false);
      histos.add("hdndeta0_5", "hdndeta0_5", kTH1D, {axisEta}, false);
    }

    if (doprocessMCpTefficiency) {
      histos.add("hmcrecdndpt", "hmcrecdndpt", kTHnSparseD, {CentAxis, axisPT}, false);
      histos.add("hmcgendndpt", "hmcgendndpt", kTHnSparseD, {CentAxis, axisPT, AxisGenpTVary}, false);
    }

    if (doprocessMCcheckFakeTracks) {
      histos.add("hTracksCount", "hTracksCount", kTHnSparseD, {CentAxis, axisTracks}, false);
      auto htrack = histos.get<THnSparse>(HIST("hTracksCount"));
      auto* x2 = htrack->GetAxis(1);
      x2->SetBinLabel(1, "All tracks");
      x2->SetBinLabel(2, "Non-fake tracks");
      for (int i = 0; i < 7; i++) {
        x2->SetBinLabel(i + 3, Form("layer %d", i));
      }
    }

    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisMult, axisFT0AMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFT0CMult}, true);
      histos.add("NPVtracks_vs_FT0C", "NPVtracks_vs_FT0C", kTH2F, {axisPV, axisFT0CMult}, true);
      histos.add("GlobalMult_vs_FV0A", "GlobalMult_vs_FV0A", kTH2F, {axisMult, axisFV0AMult}, true);
      histos.add("NPVtracks_vs_GlobalMult", "NPVtracks_vs_GlobalMult", kTH2F, {axisPV, axisMult}, true);
    }
  }

  template <typename CheckCol>
  bool IsEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);

    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);

    if (IsApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);

    if (IsApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);

    if (IsApplyVertexITSTPC && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);

    if (IsApplyVertexTOFmatched && !col.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);

    if (IsApplyVertexTRDmatched && !col.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);

    if (col.centFT0C() < 0. || col.centFT0C() > 100.) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);

    if (IsApplyExtraCorrCut && col.multNTracksPV() > NPVtracksCut && col.multFT0C() < (10 * col.multNTracksPV() - FT0CCut)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);

    return true;
  }

  expressions::Filter trackSelectionProperMixed = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionITS) &&
                                                  ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true) &&
                                                  ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, trackSelectionDCAXYonly),
                                                         ncheckbit(aod::track::trackCutFlag, trackSelectionDCA));

  void processData(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {

    if (!IsEventSelected(collision)) {
      return;
    }

    histos.fill(HIST("VtxZHist"), collision.posZ());
    histos.fill(HIST("CentPercentileHist"), collision.centFT0C());
    histos.fill(HIST("hdatazvtxcent"), collision.posZ(), collision.centFT0C());

    auto NchTracks = 0;
    for (auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      histos.fill(HIST("PhiVsEtaHistBeforePhiCut"), track.phi(), track.eta());
      if (IsApplyExtraPhiCut) {
        if ((track.phi() > 3.07666 && track.phi() < 3.12661) || track.phi() <= 0.03 || track.phi() >= 6.253) {
          continue;
        }
      }
      histos.fill(HIST("PhiVsEtaHistAfterPhiCut"), track.phi(), track.eta());
      NchTracks++;
      histos.fill(HIST("CentHistInsideTrackloop"), collision.centFT0C());
      histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0C(), track.eta(), track.phi(), track.pt(), kGlobalplusITS);

      if (track.hasTPC()) {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0C(), track.eta(), track.phi(), track.pt(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0C(), track.eta(), track.phi(), track.pt(), kITSonly);
      }
    }
    histos.fill(HIST("hdatamult"), collision.posZ(), NchTracks, collision.centFT0C());
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processData, "process data CentFT0C", false);

  void processCorrelation(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {

    if (!IsEventSelected(collision)) {
      return;
    }

    if (std::abs(collision.posZ()) > VtxRange) {
      return;
    }

    histos.fill(HIST("VtxZHist"), collision.posZ());
    auto NchTracks = 0;
    for (auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      NchTracks++;
    }

    histos.fill(HIST("GlobalMult_vs_FT0A"), NchTracks, collision.multFT0A());
    histos.fill(HIST("GlobalMult_vs_FT0C"), NchTracks, collision.multFT0C());
    histos.fill(HIST("NPVtracks_vs_FT0C"), collision.multNTracksPV(), collision.multFT0C());
    histos.fill(HIST("GlobalMult_vs_FV0A"), NchTracks, collision.multFV0A());
    histos.fill(HIST("NPVtracks_vs_GlobalMult"), collision.multNTracksPV(), NchTracks);
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processCorrelation, "do correlation study in data", false);

  void processMonteCarlo(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {

    for (auto& RecCollision : RecCollisions) {

      if (!IsEventSelected(RecCollision)) {
        continue;
      }

      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), RecCollision.centFT0C());
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), RecCollision.centFT0C());

      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (auto& Rectrack : Rectrackspart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        histos.fill(HIST("MCrecPhiVsEtaHistBeforePhiCut"), Rectrack.phi(), Rectrack.eta());
        if (IsApplyExtraPhiCut) {
          if ((Rectrack.phi() > 3.07666 && Rectrack.phi() < 3.12661) || Rectrack.phi() <= 0.03 || Rectrack.phi() >= 6.253) {
            continue;
          }
        }
        histos.fill(HIST("MCrecPhiVsEtaHistAfterPhiCut"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("MCCentHistInsideTrackloop"), RecCollision.centFT0C());
        histos.fill(HIST("hmcrecdndeta"), RecCollision.posZ(), RecCollision.centFT0C(), Rectrack.eta(), Rectrack.phi());
        if (RecCollision.centFT0C() >= 0 && RecCollision.centFT0C() < 5.0) {
          histos.fill(HIST("hdndeta0_5"), Rectrack.eta());
        }
      } // track (mcrec) loop

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

        if (std::abs(pdgParticle->Charge()) < 3) {
          continue;
        }

        if (std::abs(particle.eta()) >= etaRange) {
          continue;
        }
        if (IsApplyExtraPhiCut) {
          if ((particle.phi() > 3.07666 && particle.phi() < 3.12661) || particle.phi() <= 0.03 || particle.phi() >= 6.253) {
            continue;
          }
        }
        histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), RecCollision.centFT0C(), particle.eta(), particle.phi(), kNoGenpTVar);

        if (particle.pt() < 0.1) {
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), RecCollision.centFT0C(), particle.eta(), particle.phi(), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), RecCollision.centFT0C(), particle.eta(), particle.phi(), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), RecCollision.centFT0C(), particle.eta(), particle.phi(), kGenpTup);
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), RecCollision.centFT0C(), particle.eta(), particle.phi(), kGenpTdown);
        }
      } // track (mcgen) loop
    }   // collision loop
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processMonteCarlo, "process MC CentFT0C", false);

  void processMCpTefficiency(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {

    for (auto& RecCollision : RecCollisions) {

      if (!IsEventSelected(RecCollision)) {
        continue;
      }

      if (std::abs(RecCollision.posZ()) > VtxRange) {
        continue;
      }

      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), RecCollision.centFT0C());
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), RecCollision.centFT0C());

      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (auto& Rectrack : Rectrackspart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (Rectrack.has_mcParticle()) {
          auto mcpart = Rectrack.mcParticle();
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hmcrecdndpt"), RecCollision.centFT0C(), mcpart.pt());
          }
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

        if (std::abs(pdgParticle->Charge()) < 3) {
          continue;
        }

        if (std::abs(particle.eta()) >= etaRange) {
          continue;
        }
        histos.fill(HIST("hmcgendndpt"), RecCollision.centFT0C(), particle.pt(), kNoGenpTVar);
        if (particle.pt() < 0.1) {
          histos.fill(HIST("hmcgendndpt"), RecCollision.centFT0C(), particle.pt(), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndpt"), RecCollision.centFT0C(), particle.pt(), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndpt"), RecCollision.centFT0C(), particle.pt(), kGenpTup);
          histos.fill(HIST("hmcgendndpt"), RecCollision.centFT0C(), particle.pt(), kGenpTdown);
        }
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processMCpTefficiency, "process MC pTefficiency", false);

  void processMCcheckFakeTracks(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, FilTrackMCRecTable const& RecTracks)
  {

    for (auto& RecCollision : RecCollisions) {

      if (!IsEventSelected(RecCollision)) {
        continue;
      }

      if (std::abs(RecCollision.posZ()) > VtxRange) {
        continue;
      }

      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), RecCollision.centFT0C());
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), RecCollision.centFT0C());

      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (auto& Rectrack : Rectrackspart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (!Rectrack.hasTPC()) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), RecCollision.centFT0C(), 1);
        bool IsFakeITStracks = false;
        for (int i = 0; i < 7; i++) {
          if (Rectrack.mcMask() & 1 << i) {
            IsFakeITStracks = true;
            histos.fill(HIST("hTracksCount"), RecCollision.centFT0C(), i + 3);
            break;
          }
        }
        if (IsFakeITStracks) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), RecCollision.centFT0C(), 2);
      }
    }
  }

  PROCESS_SWITCH(HeavyIonMultiplicity, processMCcheckFakeTracks, "Check Fake tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyIonMultiplicity>(cfgc)};
}
