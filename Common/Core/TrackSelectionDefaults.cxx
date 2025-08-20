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

#include "TrackSelectionDefaults.h"

#include "Common/Core/TrackSelection.h"

#include <Framework/DataTypes.h>
#include <Framework/Logger.h>

#include <cmath>

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
TrackSelection getGlobalTrackSelectionRun3ITSMatch(int matching, int passFlag)
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
    case TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibFirst:
      selectedTracks.SetRequireHitsInITSLayers(1, {0});
      break;
    default:
      LOG(fatal) << "getGlobalTrackSelectionRun3ITSMatch with undefined ITS matching";
      break;
  }
  switch (passFlag) {
    case TrackSelection::GlobalTrackRun3DCAxyCut::Default:
      break;
    case TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3:                            // Pass3 pp parameters
      selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; }); // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
      break;
    default:
      LOG(fatal) << "getGlobalTrackSelectionRun3ITSMatch with undefined DCA cut";
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
  TrackSelection selectedTracks;

  // These track selections are the same as the global track selections as of Jan 2025. Implemented seperately to prevent future
  // global track selection changes from affecting the Run 2 hybrid track selections.
  selectedTracks.SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
  selectedTracks.SetRequireITSRefit(true);
  selectedTracks.SetRequireTPCRefit(true);
  selectedTracks.SetRequireGoldenChi2(true);
  selectedTracks.SetMinNCrossedRowsTPC(70);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  selectedTracks.SetMaxChi2PerClusterTPC(4.f);
  selectedTracks.SetMaxChi2PerClusterITS(36.f);

  // These track selections are different to the global track selections as of Jan 2025.
  selectedTracks.SetPtRange(0.15f, 1000.f);
  selectedTracks.SetEtaRange(-0.9f, 0.9f);
  selectedTracks.SetMaxDcaXYPtDep([](float /*pt*/) { return 1e+10; });
  selectedTracks.SetMaxDcaXY(2.4f);
  selectedTracks.SetMaxDcaZ(3.2f);
  selectedTracks.SetRequireHitsInITSLayers(0, {0, 1}); // no minimum required number of hits in any SPD layer
  selectedTracks.SetMaxTPCFractionSharedCls(0.4f);     // This cut machinery was added since it's used in hybrid tracks Run 2

  return selectedTracks;
}
