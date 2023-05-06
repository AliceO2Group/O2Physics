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
/// \file  TrackSelectionDefaults.cxx
/// \brief Class for the definition of standard track selection objects
/// \since 20-10-2020
///

#ifndef COMMON_CORE_TRACKSELECTIONDEFAULTS_H_
#define COMMON_CORE_TRACKSELECTIONDEFAULTS_H_

#include "Framework/DataTypes.h"
#include "Common/Core/TrackSelection.h"
#include "TrackSelectionDefaults.h"

// Default track selection requiring one hit in the SPD
TrackSelection getGlobalTrackSelection()
{
  TrackSelection selectedTracks;
  selectedTracks.SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
  selectedTracks.SetPtRange(0.1f, 1e10f);
  selectedTracks.SetEtaRange(-0.8f, 0.8f);
  selectedTracks.SetRequireITSRefit(true);
  selectedTracks.SetRequireTPCRefit(true);
  selectedTracks.SetRequireGoldenChi2(true);
  selectedTracks.SetMinNCrossedRowsTPC(70);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  selectedTracks.SetMaxChi2PerClusterTPC(4.f);
  selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
  selectedTracks.SetMaxChi2PerClusterITS(36.f);
  selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
  selectedTracks.SetMaxDcaZ(2.f);
  return selectedTracks;
}

// Default track selection requiring a particular Run 3 ITS matching
TrackSelection getGlobalTrackSelectionRun3ITSMatch(int matching)
{
  TrackSelection selectedTracks = getGlobalTrackSelection();
  selectedTracks.SetTrackType(o2::aod::track::TrackTypeEnum::Track); // Requiring that this is a Run 3 track
  selectedTracks.ResetITSRequirements();
  switch (matching) {
    case TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny:
      selectedTracks.SetRequireHitsInITSLayers(1, {0, 1, 2});
      break;
    case TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo:
      selectedTracks.SetRequireHitsInITSLayers(2, {0, 1, 2});
      break;
    case TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny:
      selectedTracks.SetRequireHitsInITSLayers(1, {0, 1, 2, 3, 4, 5, 6});
      break;
    case TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers:
      selectedTracks.SetRequireHitsInITSLayers(7, {0, 1, 2, 3, 4, 5, 6});
      break;
    default:
      LOG(fatal) << "getGlobalTrackSelectionRun3ITSMatch with undefined ITS matching";
      break;
  }
  return selectedTracks;
}

// Default track selection requiring no hit in the SPD and one in the innermost
// SDD -> complementary tracks to global selection
TrackSelection getGlobalTrackSelectionSDD()
{
  TrackSelection selectedTracks = getGlobalTrackSelection();
  selectedTracks.ResetITSRequirements();
  selectedTracks.SetRequireNoHitsInITSLayers({0, 1}); // no hit in SPD layers
  selectedTracks.SetRequireHitsInITSLayers(1, {2});   // one hit in first SDD layer
  return selectedTracks;
}

// Default track selection for nuclei analysis in run3 (STILL JUST A PLACEHOLDER)
TrackSelection getGlobalTrackSelectionRun3Nuclei()
{
  return getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo);
}

// Default track selection for HF analysis (global tracks, with its points, but no tight selection for primary) in run3 (STILL JUST A PLACEHOLDER)
TrackSelection getGlobalTrackSelectionRun3HF()
{
  TrackSelection selectedTracks;
  selectedTracks.SetTrackType(o2::aod::track::TrackTypeEnum::Track); // Run 3 track asked by default
  selectedTracks.SetPtRange(0.1f, 1e10f);
  selectedTracks.SetEtaRange(-0.8f, 0.8f);
  selectedTracks.SetRequireITSRefit(true);
  selectedTracks.SetRequireTPCRefit(true);
  selectedTracks.SetRequireGoldenChi2(true);
  selectedTracks.SetMinNCrossedRowsTPC(70);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  selectedTracks.SetMaxChi2PerClusterTPC(4.f);
  selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
  selectedTracks.SetMaxChi2PerClusterITS(36.f);
  // selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
  selectedTracks.SetMaxDcaZ(2.f);
  selectedTracks.SetMaxDcaXY(0.25);

  return selectedTracks;
}

// Reduced default track selection for jet validation based on hybrid cuts for converted (based on ESD's from run 2) A02D's
TrackSelection getJEGlobalTrackSelectionRun2()
{
  TrackSelection selectedTracks = getGlobalTrackSelection();
  selectedTracks.SetRequireGoldenChi2(false);
  selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 1e+10; });
  selectedTracks.SetEtaRange(-0.9f, 0.9f);
  selectedTracks.SetMaxDcaXY(2.4f);
  selectedTracks.SetMaxDcaZ(3.2f);
  return selectedTracks;
}

#endif
