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
/// \file LFTreeCreatorNuclei.cxx
/// \brief Writer of the nuclei candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch> and Francesca Bellini <fbellini@cern.ch>
///

#include "PWGLF/DataModel/LFNucleiTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

// #include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Writes the full information in an output TTree
struct LfTreeCreatorNuclei {
  Produces<o2::aod::LfNuclEvents> tableEvents;
  Produces<o2::aod::LfCandNucleus> tableCandidate;
  Produces<o2::aod::LfCandNucleusExtra> tableCandidateExtra;
  Produces<o2::aod::LfCandNucleusMC> tableCandidateMC;
  HistogramRegistry hEvents{"hEvents", {}, OutputObjHandlingPolicy::AnalysisObject};
  TrackSelection customTrackSelection;
  void init(o2::framework::InitContext&)
  {
    if (doprocessData == true && doprocessMC == true) {
      LOGF(fatal, "Cannot enable processData and processMC at the same time. Please choose one.");
    }

    hEvents.add("eventSelection", "eventSelection", kTH1D, {{9, -0.5, 8.5}});
    auto h = hEvents.get<TH1>(HIST("eventSelection"));
    h->GetXaxis()->SetBinLabel(1, "Custom z-vertex cut");
    h->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
    h->GetXaxis()->SetBinLabel(3, "TF border cut");
    h->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
    h->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
    h->GetXaxis()->SetBinLabel(6, "Sel8 cut");
    h->GetXaxis()->SetBinLabel(7, "Not empty events");
    h->GetXaxis()->SetBinLabel(8, "|z|<10 norm");
    h->GetXaxis()->SetBinLabel(9, "With a good track");
    customTrackSelection = myTrackSelection();
  }

  // track
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -8.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +8.0, "Value of the Nsigma cut"};
  Configurable<float> ptcutLow{"ptcutLow", 0.f, "Value of the lower pt cut for the filtering (option 3)"};
  Configurable<float> ptcutHigh{"ptcutHigh", 10.f, "Value of the upper pt cut for the filtering (option 3)"};
  Configurable<float> filterDeTPC{"filterDeTPC", 15.0, "Value of the Nsigma cut for deuterons for the filtering (option 3)"};
  Configurable<float> filterHeTPC{"filterHeTPC", 15.0, "Value of the Nsigma cut for helium3 for the filtering (option 3)"};
  Configurable<int> trackSelType{"trackSelType", 0, "Option for the track cut: 0 isGlobalTrackWoDCA, 1 isGlobalTrack, 2 is for custom track selection, 3 is for filtered mode"};
  Configurable<int> nITSInnerBarrelHits{"nITSInnerBarrelHits", 0, "Option for ITS inner barrel hits maximum: 3"};

  // events
  Configurable<float> cfgHighCutVertex{"cfgHighCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgLowCutVertex{"cfgLowCutVertex", -10.0f, "Accepted z-vertex range"};
  Configurable<bool> useSel8{"useSel8", true, "Use Sel8 for run3 Event Selection"};
  Configurable<bool> TVXtrigger{"TVXtrigger", false, "Use TVX for Event Selection (default w/ Sel8)"};
  Configurable<bool> removeTFBorder{"removeTFBorder", false, "Remove TimeFrame border (default w/ Sel8)"};
  Configurable<bool> removeITSROFBorder{"removeITSROFBorder", false, "Remove ITS Read-Out Frame border (default w/ Sel8)"};

  Configurable<bool> doSkim{"doSkim", false, "Save events that contains only selected tracks (for filtered mode)"};

  // custom track cut
  Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<int> minITScl{"minITScl", 4, "Additional cut on the ITS cluster"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};

  Filter collisionFilter = (aod::collision::posZ < cfgHighCutVertex && aod::collision::posZ > cfgLowCutVertex);

  Filter etaFilter = (nabs(aod::track::eta) < cfgCutEta);
  Filter trackFilter = (trackSelType.value == 0 && requireGlobalTrackWoDCAInFilter()) ||
                       (trackSelType.value == 1 && requireGlobalTrackInFilter()) ||
                       (trackSelType.value == 2) ||
                       (trackSelType.value == 3);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFV0As>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension,
                                    aod::pidTOFbeta, aod::TOFSignal, aod::pidEvTimeFlags,
                                    aod::pidTPCLfFullPi, aod::pidTOFFullPi,
                                    aod::pidTPCLfFullKa, aod::pidTOFFullKa,
                                    aod::pidTPCLfFullPr, aod::pidTOFFullPr,
                                    aod::pidTPCLfFullDe, aod::pidTOFFullDe,
                                    aod::pidTPCLfFullTr, aod::pidTOFFullTr,
                                    aod::pidTPCLfFullHe, aod::pidTOFFullHe,
                                    aod::pidTPCLfFullAl, aod::pidTOFFullAl>;

  template <bool isMC, typename TrackType, typename CollisionType>
  bool checkQuality(CollisionType const& /*collision*/, TrackType const& tracks)
  {
    bool out = kFALSE;
    for (auto& track : tracks) {

      if (track.itsNClsInnerBarrel() < nITSInnerBarrelHits) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      if (track.tpcNClsCrossedRows() < 90) {
        continue;
      }
      if ((track.pt() < ptcutLow.value) || (track.pt() > ptcutHigh.value)) {
        continue;
      }
      if ((TMath::Abs(track.tpcNSigmaDe()) > filterDeTPC.value) && (TMath::Abs(track.tpcNSigmaHe()) > filterHeTPC.value)) {
        continue;
      }
      out = kTRUE;
    }
    return out;
  }
  template <typename trackType>
  bool isITStrack(trackType const& track)
  {
    return (track.passedITSNCls() &&
            track.passedITSChi2NDF() &&
            track.passedITSRefit() &&
            track.passedITSHits() &&
            track.hasITS());
  }
  template <typename trackType>
  bool isTPCtrack(trackType const& track)
  {
    return (track.passedTPCNCls() &&
            track.passedTPCCrossedRows() &&
            track.passedTPCCrossedRowsOverNCls() &&
            track.passedTPCChi2NDF() &&
            track.passedTPCRefit() &&
            track.hasTPC());
  }

  TrackSelection myTrackSelection()
  {
    TrackSelection selectedTracks;
    selectedTracks = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
    LOG(info) << "Customizing track cuts:";
    selectedTracks.SetRequireITSRefit(requireITS.value);
    selectedTracks.SetRequireTPCRefit(requireTPC.value);
    selectedTracks.SetRequireGoldenChi2(requireGoldenChi2.value);
    selectedTracks.SetRequireHitsInITSLayers(minITScl.value, {0, 1, 2, 3, 4, 5, 6});
    selectedTracks.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
    selectedTracks.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
    selectedTracks.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
    selectedTracks.SetMinNClustersTPC(minTPCNClsFound.value);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
    selectedTracks.print();
    return selectedTracks;
  }

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillForOneEvent(CollisionType const& collision, TrackType const& tracks)
  {
    // Filling event properties
    tableEvents(collision.numContrib(), // collision.bcId(),
                collision.posX(),
                collision.posY(),
                collision.posZ(),
                collision.centFV0A(),
                collision.centFT0M(),
                collision.sel8(),
                collision.bc().runNumber());

    // Filling candidate properties
    tableCandidate.reserve(tracks.size());
    tableCandidateExtra.reserve(tracks.size());
    if constexpr (isMC) {
      tableCandidateMC.reserve(tracks.size());
    }
    for (auto& track : tracks) {
      if (track.itsNClsInnerBarrel() < nITSInnerBarrelHits) {
        continue;
      }
      if (trackSelType.value == 2) { // custom track selection mode
        if (!customTrackSelection.IsSelected(track)) {
          continue;
        }
      }
      if (trackSelType.value == 3) { // Filtering mode
        if (!track.hasTPC()) {
          continue;
        }
        if (track.tpcNClsCrossedRows() < 90) {
          continue;
        }
        if ((track.pt() < ptcutLow.value) || (track.pt() > ptcutHigh.value)) {
          continue;
        }
        if ((TMath::Abs(track.tpcNSigmaDe()) > filterDeTPC.value) && (TMath::Abs(track.tpcNSigmaHe()) > filterHeTPC.value)) {
          continue;
        }
      }

      tableCandidate(
        tableEvents.lastIndex(),
        track.dcaXY(), track.dcaZ(),
        track.tpcNSigmaDe(), track.tpcNSigmaHe(),
        track.tofNSigmaDe(), track.tofNSigmaHe(),
        track.isEvTimeTOF(),
        track.isEvTimeT0AC(),
        track.hasTOF(),
        track.hasTRD(),
        track.tpcInnerParam(),
        track.beta(),
        track.tpcSignal(),
        track.pt(), track.eta(), track.phi(),
        track.sign(),
        track.itsClusterSizes(),
        track.itsNCls(),
        track.tpcNClsFindable(),
        track.tpcNClsFindableMinusFound(),
        track.tpcNClsFindableMinusCrossedRows(),
        track.tpcChi2NCl(),
        track.itsChi2NCl(),
        track.itsClusterMap(),
        track.isPVContributor());

      tableCandidateExtra(
        track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
        track.tpcNSigmaTr(), track.tpcNSigmaAl(),
        track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
        track.tofNSigmaTr(), track.tofNSigmaAl(),
        track.tpcExpSignalDiffPr(), track.tpcExpSignalDiffDe(), track.tpcExpSignalDiffHe(),
        track.tofExpSignalDiffPr(), track.tofExpSignalDiffDe(), track.tofExpSignalDiffHe(),
        track.tofExpMom());

      if constexpr (isMC) { // Filling MC reco information
        if (track.has_mcParticle()) {
          bool itsPassed = isITStrack(track);
          bool tpcPassed = isTPCtrack(track);
          bool hasFakeHit = false;
          for (int i = 0; i < 10; i++) { // From ITS to TPC
            if (track.mcMask() & 1 << i) {
              hasFakeHit = true;
              break;
            }
          }

          const auto& particle = track.mcParticle();
          tableCandidateMC(particle.pdgCode(),
                           particle.isPhysicalPrimary(),
                           particle.producedByGenerator(),
                           particle.getProcess(),
                           itsPassed,
                           tpcPassed,
                           particle.px(),
                           particle.py(),
                           particle.pz(),
                           hasFakeHit);

          continue;
        }
        tableCandidateMC(0, -1, -1, 0, 0, 0, 0, 0, 0, 0);
      }
    }
  }

  Preslice<soa::Filtered<TrackCandidates>> perCollision = aod::track::collisionId;

  void processData(soa::Filtered<EventCandidates> const& collisions,
                   soa::Filtered<TrackCandidates> const& tracks,
                   aod::BCs const&)
  {
    for (const auto& collision : collisions) {
      hEvents.fill(HIST("eventSelection"), 0);

      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        if (TVXtrigger)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 1);
      }

      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        if (removeTFBorder)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 2);
      }

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        if (removeITSROFBorder)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 3);
      }

      if ((collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        hEvents.fill(HIST("eventSelection"), 4);
      }

      if (useSel8 && !collision.sel8())
        continue;
      hEvents.fill(HIST("eventSelection"), 5);

      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      if (doSkim && tracksInCollision.size() == 0)
        continue;
      hEvents.fill(HIST("eventSelection"), 6);

      // Fill the norm. column with good events with |z| < 10 cm before skimming
      if (collision.posZ() < 10 && collision.posZ() > -10) {
        hEvents.fill(HIST("eventSelection"), 7);
      }

      if (doSkim && (trackSelType.value == 3) && !checkQuality<false>(collision, tracksInCollision))
        continue;
      fillForOneEvent<false>(collision, tracksInCollision);
      hEvents.fill(HIST("eventSelection"), 8);
    }
  }

  PROCESS_SWITCH(LfTreeCreatorNuclei, processData, "process Data", true);

  void processMC(soa::Filtered<soa::Join<EventCandidates, aod::McCollisionLabels>> const& collisions,
                 soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>> const& tracks,
                 aod::BCs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {

      hEvents.fill(HIST("eventSelection"), 0);

      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        if (TVXtrigger)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 1);
      }

      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        if (removeTFBorder)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 2);
      }

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        if (removeITSROFBorder)
          continue;
      } else {
        hEvents.fill(HIST("eventSelection"), 3);
      }

      if ((collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
          (collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        hEvents.fill(HIST("eventSelection"), 4);
      }

      if (useSel8 && !collision.sel8())
        continue;
      hEvents.fill(HIST("eventSelection"), 5);

      // Fill the norm. column with good events with |z| < 10 cm before skimming
      if (collision.posZ() < 10 && collision.posZ() > -10) {
        hEvents.fill(HIST("eventSelection"), 7);
      }

      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      fillForOneEvent<true>(collision, tracksInCollision);
      hEvents.fill(HIST("eventSelection"), 8);
    }
  }

  PROCESS_SWITCH(LfTreeCreatorNuclei, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorNuclei>(cfgc)};
}
