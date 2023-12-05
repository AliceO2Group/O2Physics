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
/// \file lfpropStudy.cxx
/// \since 27-11-2023
/// \author Carolina Reetz <c.reetz@cern.ch>
/// \brief QA task to study properties of propagated tracks

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksType = soa::Join<aod::StoredTracks, aod::TracksDCA, aod::TracksExtra>;
using TracksCovMatType = soa::Join<TracksType, aod::StoredTracksCov, aod::TracksDCACov>;
using TracksLabeledType = soa::Join<TracksType, aod::McTrackLabels>;

struct lfpropStudy {

  ConfigurableAxis axisDCAxy{"axisDCAxy", {1000, -10.f, 10.f}, "DCA_{xy} (cm)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {1000, -10.f, 10.f}, "DCA_{z} (cm)"};
  ConfigurableAxis axisMom{"axisMom", {1000, 0.0f, 10.f}, "p (GeV/c)"};

  Configurable<bool> d_noITS{"d_noITS", true, "flag to select tracks without ITS"};
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum track momentum"};
  Configurable<float> d_TPCrowsMin{"d_TPCrowsMin", 70, "minimum number of TPC crossed rows"};
  Configurable<float> d_TPCrowsOverFindMin{"d_TPCrowsOverFindMin", 0.8, "minimum for ratio of TPC crossed rows over findable"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = (aod::evsel::sel8 == true);

  void init(InitContext const&)
  {
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0.0f, 1.0f}});

    const AxisSpec AxisPx{axisMom, "#it{p}_{x} (GeV/#it{c})"};
    const AxisSpec AxisPy{axisMom, "#it{p}_{y} (GeV/#it{c})"};
    const AxisSpec AxisPz{axisMom, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec AxisPt{axisMom, "#it{p}_{T} (GeV/#it{c})"};

    // momentum
    histos.add("hPx", "hPx", kTH1F, {AxisPx});
    histos.add("hPy", "hPy", kTH1F, {AxisPy});
    histos.add("hPz", "hPz", kTH1F, {AxisPz});
    histos.add("hPt", "hPt", kTH1F, {AxisPt});
    histos.add("hPxEl", "hPxEl", kTH1F, {AxisPx});
    histos.add("hPyEl", "hPyEl", kTH1F, {AxisPy});
    histos.add("hPzEl", "hPzEl", kTH1F, {AxisPz});
    histos.add("hPtEl", "hPtEl", kTH1F, {AxisPt});
    histos.add("hPxMu", "hPxMu", kTH1F, {AxisPx});
    histos.add("hPyMu", "hPyMu", kTH1F, {AxisPy});
    histos.add("hPzMu", "hPzMu", kTH1F, {AxisPz});
    histos.add("hPtMu", "hPtMu", kTH1F, {AxisPt});
    histos.add("hPxPi", "hPxPi", kTH1F, {AxisPx});
    histos.add("hPyPi", "hPyPi", kTH1F, {AxisPy});
    histos.add("hPzPi", "hPzPi", kTH1F, {AxisPz});
    histos.add("hPtPi", "hPtPi", kTH1F, {AxisPt});
    histos.add("hPxKa", "hPxKa", kTH1F, {AxisPx});
    histos.add("hPyKa", "hPyKa", kTH1F, {AxisPy});
    histos.add("hPzKa", "hPzKa", kTH1F, {AxisPz});
    histos.add("hPtKa", "hPtKa", kTH1F, {AxisPt});
    histos.add("hPxPr", "hPxPr", kTH1F, {AxisPx});
    histos.add("hPyPr", "hPyPr", kTH1F, {AxisPy});
    histos.add("hPzPr", "hPzPr", kTH1F, {AxisPz});
    histos.add("hPtPr", "hPtPr", kTH1F, {AxisPt});
    histos.add("hPxDe", "hPxDe", kTH1F, {AxisPx});
    histos.add("hPyDe", "hPyDe", kTH1F, {AxisPy});
    histos.add("hPzDe", "hPzDe", kTH1F, {AxisPz});
    histos.add("hPtDe", "hPtDe", kTH1F, {AxisPt});
    histos.add("hPxTr", "hPxTr", kTH1F, {AxisPx});
    histos.add("hPyTr", "hPyTr", kTH1F, {AxisPy});
    histos.add("hPzTr", "hPzTr", kTH1F, {AxisPz});
    histos.add("hPtTr", "hPtTr", kTH1F, {AxisPt});
    histos.add("hPxHe", "hPxHe", kTH1F, {AxisPx});
    histos.add("hPyHe", "hPyHe", kTH1F, {AxisPy});
    histos.add("hPzHe", "hPzHe", kTH1F, {AxisPz});
    histos.add("hPtHe", "hPtHe", kTH1F, {AxisPt});
    histos.add("hPxAl", "hPxAl", kTH1F, {AxisPx});
    histos.add("hPyAl", "hPyAl", kTH1F, {AxisPy});
    histos.add("hPzAl", "hPzAl", kTH1F, {AxisPz});
    histos.add("hPtAl", "hPtAl", kTH1F, {AxisPt});

    const AxisSpec AxisDCAxy{axisDCAxy, "DCA_{xy} (cm)"};
    const AxisSpec AxisDCAz{axisDCAz, "DCA_{z} (cm)"};

    // DCA
    histos.add("hDCAxy", "hDCAxy", kTH1F, {AxisDCAxy});
    histos.add("hDCAz", "hDCAz", kTH1F, {AxisDCAz});
    histos.add("hDCAxyEl", "hDCAxyEl", kTH1F, {AxisDCAxy});
    histos.add("hDCAzEl", "hDCAzEl", kTH1F, {AxisDCAz});
    histos.add("hDCAxyMu", "hDCAxyMu", kTH1F, {AxisDCAxy});
    histos.add("hDCAzMu", "hDCAzMu", kTH1F, {AxisDCAz});
    histos.add("hDCAxyPi", "hDCAxyPi", kTH1F, {AxisDCAxy});
    histos.add("hDCAzPi", "hDCAzPi", kTH1F, {AxisDCAz});
    histos.add("hDCAxyKa", "hDCAxyKa", kTH1F, {AxisDCAxy});
    histos.add("hDCAzKa", "hDCAzKa", kTH1F, {AxisDCAz});
    histos.add("hDCAxyPr", "hDCAxyPr", kTH1F, {AxisDCAxy});
    histos.add("hDCAzPr", "hDCAzPr", kTH1F, {AxisDCAz});
    histos.add("hDCAxyDe", "hDCAxyDe", kTH1F, {AxisDCAxy});
    histos.add("hDCAzDe", "hDCAzDe", kTH1F, {AxisDCAz});
    histos.add("hDCAxyTr", "hDCAxyTr", kTH1F, {AxisDCAxy});
    histos.add("hDCAzTr", "hDCAzTr", kTH1F, {AxisDCAz});
    histos.add("hDCAxyHe", "hDCAxyHe", kTH1F, {AxisDCAxy});
    histos.add("hDCAzHe", "hDCAzHe", kTH1F, {AxisDCAz});
    histos.add("hDCAxyAl", "hDCAxyAl", kTH1F, {AxisDCAxy});
    histos.add("hDCAzAl", "hDCAzAl", kTH1F, {AxisDCAz});
  }

  template <typename T>
  void fillPerTrack(const T& track)
  {
    histos.fill(HIST("hPx"), track.px());
    histos.fill(HIST("hPy"), track.py());
    histos.fill(HIST("hPz"), track.pz());
    histos.fill(HIST("hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
    histos.fill(HIST("hDCAxy"), track.dcaXY());
    histos.fill(HIST("hDCAz"), track.dcaZ());

    switch (track.pidForTracking()) {
      case o2::track::PID::Electron:
        histos.fill(HIST("hPxEl"), track.px());
        histos.fill(HIST("hPyEl"), track.py());
        histos.fill(HIST("hPzEl"), track.pz());
        histos.fill(HIST("hPtEl"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyEl"), track.dcaXY());
        histos.fill(HIST("hDCAzEl"), track.dcaZ());
        break;
      case o2::track::PID::Muon:
        histos.fill(HIST("hPxMu"), track.px());
        histos.fill(HIST("hPyMu"), track.py());
        histos.fill(HIST("hPzMu"), track.pz());
        histos.fill(HIST("hPtMu"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyMu"), track.dcaXY());
        histos.fill(HIST("hDCAzMu"), track.dcaZ());
        break;
      case o2::track::PID::Pion:
        histos.fill(HIST("hPxPi"), track.px());
        histos.fill(HIST("hPyPi"), track.py());
        histos.fill(HIST("hPzPi"), track.pz());
        histos.fill(HIST("hPtPi"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyPi"), track.dcaXY());
        histos.fill(HIST("hDCAzPi"), track.dcaZ());
        break;
      case o2::track::PID::Kaon:
        histos.fill(HIST("hPxKa"), track.px());
        histos.fill(HIST("hPyKa"), track.py());
        histos.fill(HIST("hPzKa"), track.pz());
        histos.fill(HIST("hPtKa"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyKa"), track.dcaXY());
        histos.fill(HIST("hDCAzKa"), track.dcaZ());
        break;
      case o2::track::PID::Proton:
        histos.fill(HIST("hPxPr"), track.px());
        histos.fill(HIST("hPyPr"), track.py());
        histos.fill(HIST("hPzPr"), track.pz());
        histos.fill(HIST("hPtPr"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyPr"), track.dcaXY());
        histos.fill(HIST("hDCAzPr"), track.dcaZ());
        break;
      case o2::track::PID::Deuteron:
        histos.fill(HIST("hPxDe"), track.px());
        histos.fill(HIST("hPyDe"), track.py());
        histos.fill(HIST("hPzDe"), track.pz());
        histos.fill(HIST("hPtDe"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyDe"), track.dcaXY());
        histos.fill(HIST("hDCAzDe"), track.dcaZ());
        break;
      case o2::track::PID::Triton:
        histos.fill(HIST("hPxTr"), track.px());
        histos.fill(HIST("hPyTr"), track.py());
        histos.fill(HIST("hPzTr"), track.pz());
        histos.fill(HIST("hPtTr"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyTr"), track.dcaXY());
        histos.fill(HIST("hDCAzTr"), track.dcaZ());
        break;
      case o2::track::PID::Helium3:
        histos.fill(HIST("hPxHe"), track.px());
        histos.fill(HIST("hPyHe"), track.py());
        histos.fill(HIST("hPzHe"), track.pz());
        histos.fill(HIST("hPtHe"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyHe"), track.dcaXY());
        histos.fill(HIST("hDCAzHe"), track.dcaZ());
        break;
      case o2::track::PID::Alpha:
        histos.fill(HIST("hPxAl"), track.px());
        histos.fill(HIST("hPyAl"), track.py());
        histos.fill(HIST("hPzAl"), track.pz());
        histos.fill(HIST("hPtAl"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("hDCAxyAl"), track.dcaXY());
        histos.fill(HIST("hDCAzAl"), track.dcaZ());
        break;
      default:
        LOG(fatal) << "Unknown PID: " << track.pidForTracking();
    }
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   TracksType const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      // track selection
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.trackType() == aod::track::TrackIU)
        continue; // only look at tracks which were propagated
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      fillPerTrack(track);
    }
  }
  PROCESS_SWITCH(lfpropStudy, processData, "process data", true);

  void processDataCovMat(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                         TracksCovMatType const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      // track selection
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.trackType() == aod::track::TrackIU)
        continue; // only look at tracks which were propagated
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      fillPerTrack(track);
    }
  }
  PROCESS_SWITCH(lfpropStudy, processDataCovMat, "process data with Cov. Mat.", false);

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                 TracksLabeledType const& tracks,
                 aod::McParticles const& particlesMC)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      if (!track.has_mcParticle() || track.mcParticleId() <= -1 || track.mcParticleId() > particlesMC.size())
        continue;
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.trackType() == aod::track::TrackIU)
        continue; // only look at tracks which were propagated

      histos.fill(HIST("hPx"), track.px());
      histos.fill(HIST("hPy"), track.py());
      histos.fill(HIST("hPz"), track.pz());
      histos.fill(HIST("hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
      histos.fill(HIST("hDCAxy"), track.dcaXY());
      histos.fill(HIST("hDCAz"), track.dcaZ());

      switch (std::abs(track.mcParticle().pdgCode())) {
        case kElectron:
          histos.fill(HIST("hPxEl"), track.px());
          histos.fill(HIST("hPyEl"), track.py());
          histos.fill(HIST("hPzEl"), track.pz());
          histos.fill(HIST("hPtEl"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyEl"), track.dcaXY());
          histos.fill(HIST("hDCAzEl"), track.dcaZ());
          break;
        case kMuonPlus:
          histos.fill(HIST("hPxMu"), track.px());
          histos.fill(HIST("hPyMu"), track.py());
          histos.fill(HIST("hPzMu"), track.pz());
          histos.fill(HIST("hPtMu"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyMu"), track.dcaXY());
          histos.fill(HIST("hDCAzMu"), track.dcaZ());
          break;
        case kPiPlus:
          histos.fill(HIST("hPxPi"), track.px());
          histos.fill(HIST("hPyPi"), track.py());
          histos.fill(HIST("hPzPi"), track.pz());
          histos.fill(HIST("hPtPi"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyPi"), track.dcaXY());
          histos.fill(HIST("hDCAzPi"), track.dcaZ());
          break;
        case kKPlus:
          histos.fill(HIST("hPxKa"), track.px());
          histos.fill(HIST("hPyKa"), track.py());
          histos.fill(HIST("hPzKa"), track.pz());
          histos.fill(HIST("hPtKa"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyKa"), track.dcaXY());
          histos.fill(HIST("hDCAzKa"), track.dcaZ());
          break;
        case kProton:
          histos.fill(HIST("hPxPr"), track.px());
          histos.fill(HIST("hPyPr"), track.py());
          histos.fill(HIST("hPzPr"), track.pz());
          histos.fill(HIST("hPtPr"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyPr"), track.dcaXY());
          histos.fill(HIST("hDCAzPr"), track.dcaZ());
          break;
        case 1000010020:
          histos.fill(HIST("hPxDe"), track.px());
          histos.fill(HIST("hPyDe"), track.py());
          histos.fill(HIST("hPzDe"), track.pz());
          histos.fill(HIST("hPtDe"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyDe"), track.dcaXY());
          histos.fill(HIST("hDCAzDe"), track.dcaZ());
          break;
        case 1000010030:
          histos.fill(HIST("hPxTr"), track.px());
          histos.fill(HIST("hPyTr"), track.py());
          histos.fill(HIST("hPzTr"), track.pz());
          histos.fill(HIST("hPtTr"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyTr"), track.dcaXY());
          histos.fill(HIST("hDCAzTr"), track.dcaZ());
          break;
        case 1000020030:
          histos.fill(HIST("hPxHe"), track.px());
          histos.fill(HIST("hPyHe"), track.py());
          histos.fill(HIST("hPzHe"), track.pz());
          histos.fill(HIST("hPtHe"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyHe"), track.dcaXY());
          histos.fill(HIST("hDCAzHe"), track.dcaZ());
          break;
        case 1000020040:
          histos.fill(HIST("hPxAl"), track.px());
          histos.fill(HIST("hPyAl"), track.py());
          histos.fill(HIST("hPzAl"), track.pz());
          histos.fill(HIST("hPtAl"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("hDCAxyAl"), track.dcaXY());
          histos.fill(HIST("hDCAzAl"), track.dcaZ());
          break;
        default:
          continue;
      }
    }
  }
  PROCESS_SWITCH(lfpropStudy, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<lfpropStudy>(cfgc)}; }
