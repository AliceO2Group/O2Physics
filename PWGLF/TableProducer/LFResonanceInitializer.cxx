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

/// \file LFResonanceInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers
struct reso2initializer {

  Produces<aod::ResoCollisions> resoCollisions;
  Produces<aod::ResoDaughters> reso2tracks;

  // Configurables
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"}; // Choose if running on converted data or pilot beam
  Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0s"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC and TOF PID cut (loose, improve performance)"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> isRun2{"isRun2", 0, "if Run2: demand TPC refit"};

  /// DCA Selections for V0
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.05, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"}; // Pre-selection
  Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"}; // Pre-selection
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  Configurable<double> cMinV0Radius{"cMinV0Radius", 5.0, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};

  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hGoodTrackIndices", "hGoodTrackIndices", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
                                             {"hGoodV0Indices", "hGoodV0Indices", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
                                           },
                               OutputObjHandlingPolicy::QAObject};

  // Pre-filters for efficient process
  // Filter tofPIDFilter = aod::track::tofExpMom < 0.f || ((aod::track::tofExpMom > 0.f) && ((nabs(aod::pidtof::tofNSigmaPi) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < pidnSigmaPreSelectionCut))); // TOF
  Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta;                                                                                                                                                    // Eta cut
  Filter trackCutFilter = requireGlobalTrackInFilter();                                                                                                                                                      // Global track cuts

  // Filter for all tracks
  template <bool isMC, typename CollisionType, typename TrackType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // Dummy lists for V0s
    // Loop over tracks
    for (auto& track : tracks) {
      // Tracks are already filtered by the pre-filters
      qaRegistry.fill(HIST("hGoodTrackIndices"), 0.5);

      // Add PID selection criteria here
      uint8_t tpcPIDselections = 0;
      uint8_t tofPIDselections = 0;
      // TPC PID
      if (std::abs(track.tpcNSigmaPi()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kPion;
      if (std::abs(track.tpcNSigmaKa()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kKaon;
      if (std::abs(track.tpcNSigmaPr()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kProton;
      // TOF PID
      if (track.hasTOF()) {
        tofPIDselections |= aod::resodaughter::PDGtype::kHasTOF;
        if (std::abs(track.tofNSigmaPi()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kPion;
        if (std::abs(track.tofNSigmaKa()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kKaon;
        if (std::abs(track.tofNSigmaPr()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kProton;
      }

      reso2tracks(resoCollisions.lastIndex(),
                  track.pt(),
                  track.px(),
                  track.py(),
                  track.pz(),
                  track.eta(),
                  track.phi(),
                  aod::resodaughter::DaughterType::kTrack,
                  track.dcaXY(),
                  childIDs,
                  track.sign(),
                  (uint8_t)track.tpcNClsCrossedRows(),
                  track.dcaXY(),
                  track.dcaZ(),
                  track.x(),
                  track.alpha(),
                  tpcPIDselections,
                  tofPIDselections,
                  track.tpcNSigmaPi(),
                  track.tpcNSigmaKa(),
                  track.tpcNSigmaPr(),
                  track.tofNSigmaPi(),
                  track.tofNSigmaKa(),
                  track.tofNSigmaPr(),
                  0, 0, 0,
                  0, 0, 0, 0);
    }
  }

  // Filter for all V0s
  template <bool isMC, typename CollisionType, typename TrackType>
  void fillV0s(CollisionType const& collision, aod::V0Datas const& v0s, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto& v0 : v0s) {
      qaRegistry.fill(HIST("hGoodV0Indices"), 0.5);
      auto postrack = v0.posTrack_as<TrackType>();
      auto negtrack = v0.negTrack_as<TrackType>();

      if (postrack.tpcNClsCrossedRows() < mincrossedrows)
        continue;
      if (negtrack.tpcNClsCrossedRows() < mincrossedrows)
        continue;
      qaRegistry.fill(HIST("hGoodV0Indices"), 1.5);

      if (fabs(postrack.dcaXY()) < cMinV0PosDCArToPVcut)
        continue;
      if (fabs(negtrack.dcaXY()) < cMinV0NegDCArToPVcut)
        continue;
      qaRegistry.fill(HIST("hGoodV0Indices"), 2.5);

      if ((v0.v0radius() > cMaxV0Radius) || (v0.v0radius() < cMinV0Radius))
        continue;
      qaRegistry.fill(HIST("hGoodV0Indices"), 3.5);
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinV0CosPA)
        continue;
      qaRegistry.fill(HIST("hGoodV0Indices"), 4.5);
      childIDs[0] = v0.posTrackId();
      childIDs[1] = v0.negTrackId();
      reso2tracks(resoCollisions.lastIndex(),
                  v0.pt(),
                  v0.px(),
                  v0.py(),
                  v0.pz(),
                  v0.eta(),
                  v0.phi(),
                  aod::resodaughter::DaughterType::kV0,
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                  childIDs,
                  0,
                  0,
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                  0,
                  v0.x(), 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  v0.dcaV0daughters(), v0.mLambda(), v0.mAntiLambda(),
                  v0.v0radius(), v0.x(), v0.y(), v0.z());
    }
  }

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);
  }

  void process(const soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator& collision,
               soa::Filtered<aod::Reso2TracksPIDExt> const& tracks, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    // Loop over tracks
    fillTracks<false>(collision, tracks);
    /// V0s
    if (ConfStoreV0)
      fillV0s<false>(collision, V0s, tracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reso2initializer>(cfgc, TaskName{"lf-reso2initializer"}),
  };
}
