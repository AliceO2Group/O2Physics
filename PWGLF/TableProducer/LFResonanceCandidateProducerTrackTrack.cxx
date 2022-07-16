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

/// \file LFResonanceCandidateProducerTrackTrack.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
/// Inspired by lambdakzerofinder.cxx,
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Reconstruction of track-track decay resonance candidates
struct reso2trktrkbuilder {
  Produces<aod::StoredReso2TrackTrackDatas> reso2trktrkdata;

  HistogramRegistry registry{
    "registry",
    {
      {"hReso2trktrckCandidate", "hReso2trktrckCandidate", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hCandidateDCAd", "hCandidateDCAd", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hCandidateCPA", "hCandidateCPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}}},
      {"hCandidateRadius", "hCandidateRadius", {HistType::kTH1F, {{200, 0.0f, 200.0f}}}},
    },
    OutputObjHandlingPolicy::QAObject
  };

  // Configurables
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"}; // Choose if running on converted data or pilot beam
  /// Selection criteria
  Configurable<int> selectUnLikeSignOnly{"selectUnLikeSignOnly", 1, "Select only unlike sign pair for resonance"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 5.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  /// Partition for firstTrack
  Partition<aod::ResoDaughters> parts1 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack))
                                        && (nabs(aod::resodaughter::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut))
                                        && (nabs(aod::resodaughter::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) 
                                        && (nabs(aod::resodaughter::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts  
  Partition<aod::ResoDaughters> parts2 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack))
                                        && (nabs(aod::resodaughter::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut))
                                        && (nabs(aod::resodaughter::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) 
                                        && (nabs(aod::resodaughter::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  
  void process(aod::ResoCollision& collision,
               aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    // LOGF(info, "event id: %d", collision.bcId());
    auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
    auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());

    for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {
      registry.fill(HIST("hReso2trktrckCandidate"), 0.5);
      // Un-like sign pair only
      if (selectUnLikeSignOnly && (trk1.sign() * trk2.sign() > 0))
        continue;
      registry.fill(HIST("hReso2trktrckCandidate"), 1.5);
      reso2trktrkdata(
        trk1.globalIndex(),
        trk2.globalIndex(),
        trk1.resoCollisionId(),
        trk1.sign(), trk2.sign(),
        trk1.px(), trk1.py(), trk1.pz(),
        trk2.px(), trk2.py(), trk2.pz(),
        trk1.dcaXY(), trk2.dcaXY(),
        trk1.dcaZ(), trk2.dcaZ());
    }
  }
};

/// Extends the v0data table with expression columns
struct reso2trktrkinitializer {
  Spawns<aod::Reso2TrackTrackDatas> reso2tracktrackdatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reso2trktrkbuilder>(cfgc, TaskName{"lf-reso2trktrkbuilder"}),
    adaptAnalysisTask<reso2trktrkinitializer>(cfgc, TaskName{"lf-reso2trktrkinitializer"})};
}