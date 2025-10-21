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
/// \file   trackselection.cxx
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
///
/// \brief Task performing basic track selection.
///

#include "Common/Core/TrackSelection.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//****************************************************************************************
/**
 * Produce track filter table.
 */
//****************************************************************************************
struct TrackSelectionTask {
  // FIXME: this will be removed once we can get this via meta data
  Configurable<bool> isRun3{"isRun3", false, "temp option to enable run3 mode"};
  Configurable<int> produceTable{"produceTable", -1, "option to produce the standard table table with the track selection. -1 autosetup, 0 dislabled, 1 enabled"};
  Configurable<int> produceFBextendedTable{"produceFBextendedTable", -1, "option to produce table with FB selection information. -1 autosetup, 0 dislabled, 1 enabled"};
  Configurable<bool> compatibilityIU{"compatibilityIU", false, "compatibility option to allow the processing of tracks before the introduction of IU tracks"};
  Configurable<int> itsMatching{"itsMatching", 0, "condition for ITS matching (0: Run2 SPD kAny, 1: Run3ITSibAny, 2: Run3ITSallAny, 3: Run3ITSall7Layers, 4: Run3ITSibFirst)"};
  Configurable<int> dcaSetup{"dcaSetup", 0, "dca setup: (0: default, 1: ppPass3)"};
  Configurable<float> ptMin{"ptMin", 0.1f, "Lower cut on pt for the track selected"};
  Configurable<float> ptMax{"ptMax", 1e10f, "Upper cut on pt for the track selected"};
  Configurable<float> etaMin{"etaMin", -0.8, "Lower cut on eta for the track selected"};
  Configurable<float> etaMax{"etaMax", 0.8, "Upper cut on eta for the track selected"};

  Produces<aod::TrackSelection> filterTable;
  Produces<aod::TrackSelectionExtension> filterTableDetail;
  TrackSelection globalTracks;
  TrackSelection globalTracksSDD;
  TrackSelection filtBit1;
  TrackSelection filtBit2;
  TrackSelection filtBit3;
  TrackSelection filtBit4;
  TrackSelection filtBit5;

  void init(InitContext& initContext)
  {
    // Check which tables are used
    enableFlagIfTableRequired(initContext, "TrackSelection", produceTable);
    enableFlagIfTableRequired(initContext, "TrackSelectionExtension", produceFBextendedTable);

    // Set up the track cuts
    switch (itsMatching) {
      case 0:
        // Run 2 SPD kAny
        if (!isRun3) {
          LOG(info) << "setting up globalTracks = getGlobalTrackSelection();";
          globalTracks = getGlobalTrackSelection();
          break;
        }
        LOG(warning) << "isRun3 == true and itsMatching == 0: not setting globalTracks = getGlobalTrackSelection();, but going to itsMatching == 1 and set getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny)";
        [[fallthrough]];
      case 1:
        // Run 3 kAny on 3 IB layers of ITS
        if (isRun3) {
          LOG(info) << "setting up getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, " << dcaSetup.value << ");";
          globalTracks = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, dcaSetup.value);
          break;
        }
        [[fallthrough]];
      case 2:
        // Run 3 kAny on all 7 layers of ITS
        if (isRun3) {
          LOG(info) << "setting up getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, " << dcaSetup.value << ");";
          globalTracks = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, dcaSetup.value);
          break;
        }
        [[fallthrough]];
      case 3:
        // Run 3 kAll on all 7 layers of ITS
        if (isRun3) {
          LOG(info) << "setting up getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, " << dcaSetup.value << ");";
          globalTracks = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, dcaSetup.value);
          break;
        }
        [[fallthrough]];
      case 4:
        // Run 3 kFirst, i.e. 1 hit in first layer of ITS
        if (isRun3) {
          LOG(info) << "setting up getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibFirst, " << dcaSetup.value << ");";
          globalTracks = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibFirst, dcaSetup.value);
          break;
        }
        [[fallthrough]];
      default:
        LOG(fatal) << "TrackSelectionTask with undefined cuts. Fix it!";
        break;
    }
    globalTracks.SetPtRange(ptMin, ptMax);
    globalTracks.SetEtaRange(etaMin, etaMax);
    if (isRun3) {
      globalTracks.SetTrackType(o2::aod::track::TrackTypeEnum::Track); // Requiring that this is a Run 3 track
      if (compatibilityIU.value) {                                     // If in compatibility mode tracks are asked to be IU tracks
        globalTracks.SetTrackType(o2::aod::track::TrackTypeEnum::TrackIU);
      }
    }
    globalTracks.print();
    // Extra requirement on the ITS -> Run 2: asking for 1 hit SDD and no hit in SPD
    globalTracksSDD = getGlobalTrackSelectionSDD();
    globalTracksSDD.SetPtRange(ptMin, ptMax);
    globalTracksSDD.SetEtaRange(etaMin, etaMax);

    filtBit1 = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny);

    filtBit2 = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo);

    filtBit3 = getGlobalTrackSelectionRun3HF();

    filtBit4 = getGlobalTrackSelectionRun3Nuclei();

    LOG(info) << "setting up filtBit5 = getJEGlobalTrackSelectionRun2();";
    filtBit5 = getJEGlobalTrackSelectionRun2(); // Jet validation requires reduced set of cuts
  }

  void process(soa::Join<aod::FullTracks, aod::TracksDCA> const& tracks)
  {
    if (produceTable == 1) {
      filterTable.reserve(tracks.size());
    }
    if (produceFBextendedTable == 1) {
      filterTableDetail.reserve(tracks.size());
    }
    if (produceTable == 0 && produceFBextendedTable == 0) {
      return;
    }
    if (isRun3) {
      for (const auto& track : tracks) {

        if (produceTable == 1) {
          filterTable((uint8_t)0,
                      globalTracks.IsSelectedMask(track),
                      filtBit1.IsSelected(track),
                      filtBit2.IsSelected(track),
                      filtBit3.IsSelected(track),
                      filtBit4.IsSelected(track),
                      filtBit5.IsSelected(track));
        }
        if (produceFBextendedTable == 1) {
          o2::aod::track::TrackSelectionFlags::flagtype trackflagGlob = globalTracks.IsSelectedMask(track);
          o2::aod::track::TrackSelectionFlags::flagtype trackflagFB1 = filtBit1.IsSelectedMask(track);
          o2::aod::track::TrackSelectionFlags::flagtype trackflagFB2 = filtBit2.IsSelectedMask(track);
          // o2::aod::track::TrackSelectionFlags::flagtype trackflagFB3 = filtBit3.IsSelectedMask(track); // only temporarily commented, will be used
          // o2::aod::track::TrackSelectionFlags::flagtype trackflagFB4 = filtBit4.IsSelectedMask(track);
          // o2::aod::track::TrackSelectionFlags::flagtype trackflagFB5 = filtBit5.IsSelectedMask(track);

          filterTableDetail(o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTrackType),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kPtRange),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kEtaRange),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCNCls),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCCrossedRows),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCCrossedRowsOverNCls),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCChi2NDF),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCRefit),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSNCls),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSChi2NDF),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSRefit),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSHits),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kGoldenChi2),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kDCAxy),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kDCAz),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagFB1, o2::aod::track::TrackSelectionFlags::kITSHits),
                            o2::aod::track::TrackSelectionFlags::checkFlag(trackflagFB2, o2::aod::track::TrackSelectionFlags::kITSHits));
        }
      }
      return;
    }

    for (const auto& track : tracks) {
      o2::aod::track::TrackSelectionFlags::flagtype trackflagGlob = globalTracks.IsSelectedMask(track);
      if (produceTable == 1) {
        filterTable((uint8_t)globalTracksSDD.IsSelected(track),
                    globalTracks.IsSelectedMask(track),
                    filtBit1.IsSelected(track),
                    filtBit2.IsSelected(track),
                    filtBit3.IsSelected(track),
                    filtBit4.IsSelected(track),
                    filtBit5.IsSelected(track));
      }
      if (produceFBextendedTable == 1) {
        filterTableDetail(o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTrackType),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kPtRange),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kEtaRange),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCNCls),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCCrossedRows),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCCrossedRowsOverNCls),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCChi2NDF),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kTPCRefit),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSNCls),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSChi2NDF),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSRefit),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kITSHits),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kGoldenChi2),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kDCAxy),
                          o2::aod::track::TrackSelectionFlags::checkFlag(trackflagGlob, o2::aod::track::TrackSelectionFlags::kDCAz),
                          0,
                          0);
      }
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackSelectionTask>(cfgc, TaskName{"track-selection"})};
  return workflow;
}
