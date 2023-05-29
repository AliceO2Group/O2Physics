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

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"

// #include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Writes the full information in an output TTree
struct LfTreeCreatorNuclei {
  Produces<o2::aod::LfCandNucleusFullEvents> tableEvents;
  Produces<o2::aod::LfCandNucleusFull> tableCandidate;
  Produces<o2::aod::LfCandNucleusMC> tableCandidateMC;

  void init(o2::framework::InitContext&)
  {
    if (doprocessData == true && doprocessMC == true) {
      LOGF(fatal, "Cannot enable processData and processMC at the same time. Please choose one.");
    }
  }

  // track
  Configurable<float> yCut{"yMin", 1.f, "Rapidity cut"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -8.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +8.0, "Value of the Nsigma cut"};
  Configurable<int> trackSelType{"trackSelType", 0, "Option for the track cut: 0 isGlobalTrackWoDCA, 1 isGlobalTrack"};
  Configurable<int> nITSInnerBarrelHits{"nITSInnerBarrelHits", 0, "Option for ITS inner barrel hits maximum: 3"};

  // events
  Configurable<float> cfgHighCutVertex{"cfgHighCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgLowCutVertex{"cfgLowCutVertex", -10.0f, "Accepted z-vertex range"};
  Configurable<bool> useEvsel{"useEvsel", true, "Use sel8 for run3 Event Selection"};

  Filter collisionFilter = (aod::collision::posZ < cfgHighCutVertex && aod::collision::posZ > cfgLowCutVertex);
  // Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (requireGlobalTrackInFilter());
  Filter etaFilter = (nabs(aod::track::eta) < cfgCutEta);
  Filter trackFilter = (trackSelType.value == 0 && requireGlobalTrackWoDCAInFilter()) || (trackSelType.value == 1 && requireGlobalTrackInFilter());
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTOFbeta, aod::TOFSignal, aod::pidEvTimeFlags,
                                    aod::pidTPCFullPi, aod::pidTOFFullPi,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa,
                                    aod::pidTPCFullPr, aod::pidTOFFullPr,
                                    aod::pidTPCFullDe, aod::pidTOFFullDe,
                                    aod::pidTPCFullTr, aod::pidTOFFullTr,
                                    aod::pidTPCFullHe, aod::pidTOFFullHe,
                                    aod::pidTPCFullAl, aod::pidTOFFullAl>;

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillForOneEvent(CollisionType const& collision, TrackType const& tracks)
  {
    // Filling event properties
    tableEvents(collision.bcId(),
                collision.numContrib(),
                collision.posX(),
                collision.posY(),
                collision.posZ(),
                collision.multFV0M(),
                collision.sel8(),
                collision.bc().runNumber());

    // Filling candidate properties
    tableCandidate.reserve(tracks.size());
    if constexpr (isMC) {
      tableCandidateMC.reserve(tracks.size());
    }
    for (auto& track : tracks) {
      if (track.itsNClsInnerBarrel() < nITSInnerBarrelHits)
        continue;
      // auto const& mcParticle = track.mcParticle();
      tableCandidate(
        tableEvents.lastIndex(),
        track.dcaXY(),
        track.dcaZ(),
        track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
        track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe(), track.tpcNSigmaAl(),
        track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
        track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe(), track.tofNSigmaAl(),
        track.tpcExpSignalDiffPr(), track.tpcExpSignalDiffDe(), track.tpcExpSignalDiffHe(),
        track.tofExpSignalDiffPr(), track.tofExpSignalDiffDe(), track.tofExpSignalDiffHe(),
        track.isEvTimeTOF(),
        track.isEvTimeT0AC(),
        track.hasTOF(),
        track.hasTRD(),
        track.tpcInnerParam(),
        track.tofExpMom(),
        track.tpcSignal(),
        track.beta(),
        track.px(),
        track.py(),
        track.pz(),
        track.pt(),
        track.p(),
        track.eta(),
        track.phi(),
        track.sign(),
        track.itsNCls(),
        track.tpcNClsCrossedRows(),
        track.tpcCrossedRowsOverFindableCls(),
        track.tpcNClsFound(),
        track.tpcChi2NCl(),
        track.itsChi2NCl(),
        track.itsClusterMap(),
        track.isPVContributor());

      if constexpr (isMC) { // Filling MC reco information
        if (track.has_mcParticle()) {
          const auto& particle = track.mcParticle();
          tableCandidateMC(particle.pdgCode(),
                           particle.isPhysicalPrimary(),
                           particle.producedByGenerator());
          continue;
        }
        tableCandidateMC(0, -1, -1);
      }
    }
  }

  Preslice<soa::Filtered<TrackCandidates>> perCollision = aod::track::collisionId;

  void processData(soa::Filtered<EventCandidates> const& collisions,
                   soa::Filtered<TrackCandidates> const& tracks, aod::BCs const&)
  {
    for (const auto& collision : collisions) {
      if (useEvsel && !collision.sel8()) {
        continue;
      }
      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      fillForOneEvent<false>(collision, tracksInCollision);
    }
  }

  PROCESS_SWITCH(LfTreeCreatorNuclei, processData, "process Data", true);

  void processMC(soa::Filtered<soa::Join<EventCandidates, aod::McCollisionLabels>> const& collisions,
                 soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>> const& tracks,
                 aod::BCs const&, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (useEvsel && !collision.sel8()) {
        continue;
      }
      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      fillForOneEvent<true>(collision, tracksInCollision);
    }
  }

  PROCESS_SWITCH(LfTreeCreatorNuclei, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorNuclei>(cfgc)};
}
