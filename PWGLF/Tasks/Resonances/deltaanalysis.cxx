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
// Analysis task for delta analysis

#include "PWGLF/DataModel/LFLithium4Tables.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

namespace
{
constexpr float pionMass = o2::constants::physics::MassPionCharged;
constexpr float protonMass = o2::constants::physics::MassProton;
constexpr int deltaPlusPlusPDG = 2224;
constexpr int deltaZeroPDG = 2114;
constexpr int protonPDG = 2212;
constexpr int pionPDG = 211;
} // namespace

struct deltaAnalysis {

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5, "Y cut on reconstructed delta"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<float> cfgCutPtProtonTPC{"cfgCutPtProtonTPC", 0.8, "Pt cut of proton in TPC"};
  Configurable<float> cfgCutPtPionTPC{"cfgCutPtPionTPC", 0.8, "Pt cut of pion in TPC"};

  // Histogram axes
  AxisSpec nSigmaTPCaxis = {100, -5., 5., "n#sigma_{TPC}"};
  AxisSpec nSigmaTOFaxis = {100, -5., 5., "n#sigma_{TOF}"};
  ConfigurableAxis cfgPtAxis{"cfgPtAxis", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis cfgMassAxis{"cfgDeltaPlusPlusAxis", {75, 1.05, 1.8}, ""};

  void init(o2::framework::InitContext&)
  {
    // Define axes
    AxisSpec deltaPlusPlusAxis{cfgMassAxis, "m(p + #pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec antiDeltaPlusPlusAxis{cfgMassAxis, "m(#bar{p} + #pi^{-}) (GeV/#it{c}^{2})"};
    AxisSpec deltaZeroAxis{cfgMassAxis, "m(p + #pi^{-}) (GeV/#it{c}^{2})"};
    AxisSpec antiDeltaZeroAxis{cfgMassAxis, "m(#bar{p} + #pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec ptAxis{cfgPtAxis, "#it{p}_{T} (GeV/#it{c})"};

    // Collision
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{2001, -0.5, 2000.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributors; n", kTH1F, {{2001, -0.5f, 2000.5f}});

    // Single track
    histos.add("hPiPlusDCAxy", "DCA_{xy} distribution for #pi^{+}; DCA_{xy} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPiPlusDCAz", "DCA_{z} distribution for #pi^{+}; DCA_{z} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPiPlusNsigmaTPCvsPt", "n#sigma_{TPC} distribution vs #it{p}_{T} for #pi^{+}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPiPlusNsigmaTPCvsPt_TPC_only", "n#sigma_{TPC} distribution vs #it{p}_{T} for #pi^{+}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPiPlusNsigmaTOFvsPt", "n#sigma_{TOF} distribution vs #it{p}_{T} for #pi^{+}", kTH2F, {ptAxis, nSigmaTOFaxis});

    histos.add("hPiMinusDCAxy", "DCA_{xy} distribution for #pi^{-}; DCA_{xy} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPiMinusDCAz", "DCA_{z} distribution for #pi^{-}; DCA_{z} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPiMinusNsigmaTPCvsPt", "n#sigma_{TPC} distribution vs #it{p}_{T} for #pi^{-}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPiMinusNsigmaTPCvsPt_TPC_only", "n#sigma_{TPC} distribution vs #it{p}_{T} for #pi^{-}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPiMinusNsigmaTOFvsPt", "n#sigma_{TOF} distribution vs #it{p}_{T} for #pi^{-}", kTH2F, {ptAxis, nSigmaTOFaxis});

    histos.add("hPrPlusDCAxy", "DCA_{xy} distribution for p; DCA_{xy} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPrPlusDCAz", "DCA_{z} distribution for p; DCA_{z} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPrPlusNsigmaTPCvsPt", "n#sigma_{TPC} distribution vs #it{p}_{T} for p", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPrPlusNsigmaTPCvsPt_TPC_only", "n#sigma_{TPC} distribution vs #it{p}_{T} for p", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPrPlusNsigmaTOFvsPt", "n#sigma_{TOF} distribution vs #it{p}_{T} for p", kTH2F, {ptAxis, nSigmaTOFaxis});

    histos.add("hPrMinusDCAxy", "DCA_{xy} distribution for #bar{p}; DCA_{xy} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPrMinusDCAz", "DCA_{z} distribution for #bar{p};  DCA_{z} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPrMinusNsigmaTPCvsPt", "n#sigma_{TPC} distribution vs #it{p}_{T} for #bar{p}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPrMinusNsigmaTPCvsPt_TPC_only", "n#sigma_{TPC} distribution vs #it{p}_{T} for #bar{p}", kTH2F, {ptAxis, nSigmaTPCaxis});
    histos.add("hPrMinusNsigmaTOFvsPt", "n#sigma_{TOF} distribution vs #it{p}_{T} for #bar{p}", kTH2F, {ptAxis, nSigmaTOFaxis});

    // Deltas
    histos.add("hDeltaPlusPlusInvMass", "Invariant mass distribution for #Delta^{++}", kTH2F, {ptAxis, deltaPlusPlusAxis});
    histos.add("hAntiDeltaPlusPlusInvMass", "Invariant mass distribution for #bar{#Delta^{++}}", kTH2F, {ptAxis, antiDeltaPlusPlusAxis});

    histos.add("hDeltaZeroInvMass", "Invariant mass distribution for #Delta^{0}", kTH2F, {ptAxis, deltaZeroAxis});
    histos.add("hAntiDeltaZeroInvMass", "Invariant mass distribution for #bar{#Delta^{0}}", kTH2F, {ptAxis, antiDeltaZeroAxis});

    if (doprocessMixedEvent) {
      // Deltas - Event mixing
      histos.add("hDeltaPlusPlusInvMassEM", "Invariant mass distribution for #Delta^{++} - event mixing", kTH2F, {ptAxis, deltaPlusPlusAxis});
      histos.add("hAntiDeltaPlusPlusInvMassEM", "Invariant mass distribution for #bar{#Delta^{++}} - event mixing", kTH2F, {ptAxis, antiDeltaPlusPlusAxis});

      histos.add("hDeltaZeroInvMassEM", "Invariant mass distribution for #Delta^{0} - event mixing", kTH2F, {ptAxis, deltaZeroAxis});
      histos.add("hAntiDeltaZeroInvMassEM", "Invariant mass distribution for #bar{#Delta^{0}} - event mixing", kTH2F, {ptAxis, antiDeltaZeroAxis});
    }

    if (doprocessMC) {
      // generated quantities
      histos.add("hDeltaPlusPlusInvMassGen", "Invariant mass distribution for #Delta^{++} - generated", kTH2F, {ptAxis, deltaPlusPlusAxis});
      histos.add("hAntiDeltaPlusPlusInvMassGen", "Invariant mass distribution for #bar{#Delta^{++}} - generated", kTH2F, {ptAxis, antiDeltaPlusPlusAxis});

      histos.add("hDeltaZeroInvMassGen", "Invariant mass distribution for #Delta^{0} - generated", kTH2F, {ptAxis, deltaZeroAxis});
      histos.add("hAntiDeltaZeroInvMassGen", "Invariant mass distribution for #bar{#Delta^{0}} - generated", kTH2F, {ptAxis, antiDeltaZeroAxis});

      histos.add("hRecProtonAntiDeltaPlusPlus", "Proton from #bar{#Delta^{++}}", kTH1F, {ptAxis});
      histos.add("hRecProtonDeltaPlusPlus", "Proton from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("hRecProtonAntiDeltaZero", "Proton from #bar{#Delta^{0}}", kTH1F, {ptAxis});
      histos.add("hRecProtonDeltaZero", "Proton from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("hRecPionAntiDeltaPlusPlus", "Pion from #bar{#Delta^{++}}", kTH1F, {ptAxis});
      histos.add("hRecPionDeltaPlusPlus", "Pion from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("hRecPionAntiDeltaZero", "Pion from #bar{#Delta^{0}}", kTH1F, {ptAxis});
      histos.add("hRecPionDeltaZero", "Pion from #Delta^{0}", kTH1F, {ptAxis});

      histos.add("hGenProtonAntiDeltaPlusPlus", "Proton from #bar{#Delta^{++}}", kTH1F, {ptAxis});
      histos.add("hGenProtonDeltaPlusPlus", "Proton from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("hGenProtonAntiDeltaZero", "Proton from #bar{#Delta^{0}}", kTH1F, {ptAxis});
      histos.add("hGenProtonDeltaZero", "Proton from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("hGenPionAntiDeltaPlusPlus", "Pion from #bar{#Delta^{++}}", kTH1F, {ptAxis});
      histos.add("hGenPionDeltaPlusPlus", "Pion from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("hGenPionAntiDeltaZero", "Pion from #bar{#Delta^{0}}", kTH1F, {ptAxis});
      histos.add("hGenPionDeltaZero", "Pion from #Delta^{0}", kTH1F, {ptAxis});
    }
  }

  template <typename T>
  bool selectionTrack(const T& track)
  {
    // if (!track.isGlobalTrack()) {
    //   return false;
    // }
    if (track.itsNCls() < 5) {
      return false;
    }
    if (track.tpcNClsShared() > 0) {
      return false;
    }
    if (track.tpcNClsFound() < 70) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable()) {
      return false;
    }
    if (track.tpcChi2NCl() > 4.f) {
      return false;
    }
    if (track.itsChi2NCl() > 36.f) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPIDPion(const T& track)
  {
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaPi()) < nsigmaCutTOF && std::abs(track.tpcNSigmaPi()) < nsigmaCutTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("hPiPlusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPi());
          histos.fill(HIST("hPiPlusNsigmaTOFvsPt"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("hPiPlusDCAxy"), track.dcaXY());
          histos.fill(HIST("hPiPlusDCAz"), track.dcaZ());
        } else {
          histos.fill(HIST("hPiMinusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPi());
          histos.fill(HIST("hPiMinusNsigmaTOFvsPt"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("hPiMinusDCAxy"), track.dcaXY());
          histos.fill(HIST("hPiMinusDCAz"), track.dcaZ());
        }
        return true;
      }
    } else if (std::abs(track.tpcNSigmaPi()) < nsigmaCutTPC) {
      if (track.sign() > 0) {
        histos.fill(HIST("hPiPlusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("hPiPlusNsigmaTPCvsPt_TPC_only"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("hPiPlusDCAxy"), track.dcaXY());
        histos.fill(HIST("hPiPlusDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("hPiMinusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("hPiMinusNsigmaTPCvsPt_TPC_only"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("hPiMinusDCAxy"), track.dcaXY());
        histos.fill(HIST("hPiMinusDCAz"), track.dcaZ());
      }
      if (track.pt() < cfgCutPtPionTPC) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDProton(const T& track)
  {
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaPr()) < nsigmaCutTOF && std::abs(track.tpcNSigmaPr()) < nsigmaCutTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("hPrPlusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPr());
          histos.fill(HIST("hPrPlusNsigmaTOFvsPt"), track.pt(), track.tofNSigmaPr());
          histos.fill(HIST("hPrPlusDCAxy"), track.dcaXY());
          histos.fill(HIST("hPrPlusDCAz"), track.dcaZ());
        } else {
          histos.fill(HIST("hPrMinusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPr());
          histos.fill(HIST("hPrMinusNsigmaTOFvsPt"), track.pt(), track.tofNSigmaPr());
          histos.fill(HIST("hPrMinusDCAxy"), track.dcaXY());
          histos.fill(HIST("hPrMinusDCAz"), track.dcaZ());
        }
        return true;
      }
    } else if (std::abs(track.tpcNSigmaPr()) < nsigmaCutTPC) {
      if (track.sign() > 0) {
        histos.fill(HIST("hPrPlusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPr());
        histos.fill(HIST("hPrPlusNsigmaTPCvsPt_TPC_only"), track.pt(), track.tpcNSigmaPr());
        histos.fill(HIST("hPrPlusDCAxy"), track.dcaXY());
        histos.fill(HIST("hPrPlusDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("hPrMinusNsigmaTPCvsPt"), track.pt(), track.tpcNSigmaPr());
        histos.fill(HIST("hPrMinusNsigmaTPCvsPt_TPC_only"), track.pt(), track.tpcNSigmaPr());
        histos.fill(HIST("hPrMinusDCAxy"), track.dcaXY());
        histos.fill(HIST("hPrMinusDCAz"), track.dcaZ());
      }
      if (track.pt() < cfgCutPtProtonTPC) {
        return true;
      }
    }
    return false;
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>>;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;

  // Binning for event mixing
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1.0, 5.0, 10., 15., 20., 30., 40., 50., 70., 100.0, 105.}, "Binning of the centrality axis"};
  ConfigurableAxis cfgVtxAxis{"cfgVtxAxis", {VARIABLE_WIDTH, -20, -15, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7, 10, 15, 20}, "Mixing bins - z-vertex"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;

  BinningType binningOnPositions{{cfgVtxAxis}, true};
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  void processSameEvent(EventCandidates const& collisions, TrackCandidates const& tracks, aod::BCs const&)
  {
    for (auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto perColTracks = tracks.sliceBy(perCol, collIdx);
      perColTracks.bindExternalIndices(&tracks);

      for (auto& [t0, t1] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(perColTracks, perColTracks))) {

        if (t0.globalIndex() == t1.globalIndex()) {
          continue;
        }
        if (!selectionTrack(t0) || !selectionTrack(t1)) {
          continue;
        }
        if (selectionPIDProton(t0)) {

          TLorentzVector proton;
          proton.SetPtEtaPhiM(t0.pt(), t0.eta(), t0.phi(), protonMass);

          if (selectionPIDPion(t1)) {

            TLorentzVector pion;
            pion.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), pionMass);

            TLorentzVector delta = proton + pion;

            if (abs(delta.Rapidity()) > cfgCutY) {
              continue;
            }

            if (t0.sign() < 0) {
              if (t1.sign() < 0) {
                histos.fill(HIST("hAntiDeltaPlusPlusInvMass"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hAntiDeltaZeroInvMass"), delta.Pt(), delta.M());
              }
            } else {
              if (t1.sign() < 0) {
                histos.fill(HIST("hDeltaZeroInvMass"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hDeltaPlusPlusInvMass"), delta.Pt(), delta.M());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(deltaAnalysis, processSameEvent, "Process same event", false);

  void processMixedEvent(EventCandidates const& /*collisions*/, TrackCandidates const& /*tracks*/)
  {
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      for (auto& [t0, t1] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t0) || !selectionTrack(t1)) {
          continue;
        }

        if (selectionPIDProton(t0)) {

          TLorentzVector proton;
          proton.SetPtEtaPhiM(t0.pt(), t0.eta(), t0.phi(), protonMass);

          if (selectionPIDPion(t1)) {

            TLorentzVector pion;
            pion.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), pionMass);

            TLorentzVector delta = proton + pion;

            if (abs(delta.Rapidity()) > cfgCutY) {
              continue;
            }

            if (t0.sign() < 0) {
              if (t1.sign() < 0) {
                histos.fill(HIST("hAntiDeltaPlusPlusInvMassEM"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hAntiDeltaZeroInvMassEM"), delta.Pt(), delta.M());
              }
            } else {
              if (t1.sign() < 0) {
                histos.fill(HIST("hDeltaZeroInvMassEM"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hDeltaPlusPlusInvMassEM"), delta.Pt(), delta.M());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(deltaAnalysis, processMixedEvent, "Process Mixed event", false);

  void processMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const&, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles)
  {

    for (auto& collision : collisions) {

      if (!collision.sel8() || std::abs(collision.posZ()) > cfgCutVertex) {
        continue;
      }

      histos.fill(HIST("hCentrality"), 1);
      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto perColTracks = tracks.sliceBy(perColMC, collIdx);
      perColTracks.bindExternalIndices(&tracks);

      for (auto& t0 : perColTracks) {

        if (!selectionTrack(t0)) {
          continue;
        }
        if (!t0.has_mcParticle()) {
          continue;
        }

        const auto mcTrack = t0.mcParticle();

        if (std::abs(mcTrack.pdgCode()) == protonPDG) {
          if (selectionPIDProton(t0)) {
            for (auto& motherTrackProton : mcTrack.mothers_as<aod::McParticles>()) {
              if (std::abs(motherTrackProton.pdgCode()) == deltaPlusPlusPDG) {
                if (t0.sign() < 0) {
                  histos.fill(HIST("hRecProtonAntiDeltaPlusPlus"), t0.pt());
                } else {
                  histos.fill(HIST("hRecProtonDeltaPlusPlus"), t0.pt());
                }
              } else if (std::abs(motherTrackProton.pdgCode()) == deltaZeroPDG) {
                if (t0.sign() < 0) {
                  histos.fill(HIST("hRecProtonAntiDeltaZero"), t0.pt());
                } else {
                  histos.fill(HIST("hRecProtonDeltaZero"), t0.pt());
                }
              }
            }
          }
        } else if (std::abs(mcTrack.pdgCode()) == pionPDG) {
          if (selectionPIDPion(t0)) {
            for (auto& motherTrackPion : mcTrack.mothers_as<aod::McParticles>()) {
              if (std::abs(motherTrackPion.pdgCode()) == deltaPlusPlusPDG) {
                if (t0.sign() < 0) {
                  histos.fill(HIST("hRecPionAntiDeltaPlusPlus"), t0.pt());
                } else {
                  histos.fill(HIST("hRecPionDeltaPlusPlus"), t0.pt());
                }
              } else if (std::abs(motherTrackPion.pdgCode()) == deltaZeroPDG) {
                if (t0.sign() < 0) {
                  histos.fill(HIST("hRecPionDeltaZero"), t0.pt());
                } else {
                  histos.fill(HIST("hRecPionAntiDeltaZero"), t0.pt());
                }
              }
            }
          }
        }
      }

      for (auto& [t0, t1] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(perColTracks, perColTracks))) {

        if (t0.globalIndex() == t1.globalIndex()) {
          continue;
        }
        if (!selectionTrack(t0) || !selectionTrack(t1)) {
          continue;
        }
        if (selectionPIDProton(t0)) {

          if (!t0.has_mcParticle()) {
            continue;
          }
          TLorentzVector proton;
          proton.SetPtEtaPhiM(t0.pt(), t0.eta(), t0.phi(), protonMass);

          if (selectionPIDPion(t1)) {

            if (!t1.has_mcParticle()) {
              continue;
            }
            TLorentzVector pion;
            pion.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), pionMass);

            // select real protons and pion
            const auto mcTrackProton = t0.mcParticle();
            const auto mcTrackPion = t1.mcParticle();

            if (std::abs(mcTrackProton.pdgCode()) != protonPDG) {
              continue;
            }

            if (std::abs(mcTrackPion.pdgCode()) != pionPDG) {
              continue;
            }

            bool found_mother = false;

            for (auto& motherTrackProton : mcTrackProton.mothers_as<aod::McParticles>()) {
              for (auto& motherTrackPion : mcTrackPion.mothers_as<aod::McParticles>()) {
                if (motherTrackProton != motherTrackPion) {
                  continue;
                }
                if (std::abs(motherTrackProton.pdgCode()) != deltaPlusPlusPDG && std::abs(motherTrackProton.pdgCode()) != deltaZeroPDG) {
                  continue;
                }
                found_mother = true;
                break;
              }
              if (found_mother) {
                break;
              }
            }

            if (!found_mother) {
              continue;
            }

            TLorentzVector delta = proton + pion;

            if (abs(delta.Rapidity()) > cfgCutY) {
              continue;
            }

            if (t0.sign() < 0) {
              if (t1.sign() < 0) {
                histos.fill(HIST("hAntiDeltaPlusPlusInvMass"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hAntiDeltaZeroInvMass"), delta.Pt(), delta.M());
              }
            } else {
              if (t1.sign() < 0) {
                histos.fill(HIST("hDeltaZeroInvMass"), delta.Pt(), delta.M());
              } else {
                histos.fill(HIST("hDeltaPlusPlusInvMass"), delta.Pt(), delta.M());
              }
            }
          }
        }
      }
    }

    for (auto& mcParticle : mcParticles) {

      int pdg = mcParticle.pdgCode();
      if (std::abs(pdg) != deltaPlusPlusPDG && std::abs(pdg) != deltaZeroPDG) {
        continue;
      }

      if (std::abs(mcParticle.y()) > cfgCutY) {
        continue;
      }

      // if (!mcParticle.isPhysicalPrimary()) {
      //   continue;
      // }

      auto daughters = mcParticle.daughters_as<aod::McParticles>();
      bool daughtPr = false;
      bool daughtPi = false;
      float ptPr = -999.;
      float ptPi = -999.;
      double eDelta = 0.;
      for (auto daught : daughters) {
        if (std::abs(daught.pdgCode()) == protonPDG) {
          daughtPr = true;
          ptPr = daught.pt();
          eDelta += daught.e();
        } else if (std::abs(daught.pdgCode()) == pionPDG) {
          daughtPi = true;
          ptPi = daught.pt();
          eDelta += daught.e();
        }
      }

      if (!daughtPr || !daughtPi) {
        continue;
      }

      float gen_mass = TMath::Sqrt(eDelta * eDelta - mcParticle.p() * mcParticle.p());

      if (pdg > 0) {
        if (std::abs(pdg) == deltaPlusPlusPDG) {
          histos.fill(HIST("hDeltaPlusPlusInvMassGen"), mcParticle.pt(), gen_mass);
          histos.fill(HIST("hGenProtonDeltaPlusPlus"), ptPr);
          histos.fill(HIST("hGenPionDeltaPlusPlus"), ptPi);
        } else if (std::abs(pdg) == deltaZeroPDG) {
          histos.fill(HIST("hDeltaZeroInvMassGen"), mcParticle.pt(), gen_mass);
          histos.fill(HIST("hGenProtonDeltaZero"), ptPr);
          histos.fill(HIST("hGenPionDeltaZero"), ptPi);
        }
      } else {
        if (std::abs(pdg) == deltaPlusPlusPDG) {
          histos.fill(HIST("hAntiDeltaPlusPlusInvMassGen"), mcParticle.pt(), gen_mass);
          histos.fill(HIST("hGenProtonAntiDeltaPlusPlus"), ptPr);
          histos.fill(HIST("hGenPionAntiDeltaPlusPlus"), ptPi);
        } else if (std::abs(pdg) == deltaZeroPDG) {
          histos.fill(HIST("hAntiDeltaZeroInvMassGen"), mcParticle.pt(), gen_mass);
          histos.fill(HIST("hGenProtonAntiDeltaZero"), ptPr);
          histos.fill(HIST("hGenPionAntiDeltaZero"), ptPi);
        }
      }
    }
  }
  PROCESS_SWITCH(deltaAnalysis, processMC, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<deltaAnalysis>(cfgc, TaskName{"deltaAnalysis"})};
}
