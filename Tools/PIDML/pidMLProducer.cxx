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
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace pidtracks
{
DECLARE_SOA_COLUMN(Px, px, float);                                 //! Non-dynamic column with track x-momentum
DECLARE_SOA_COLUMN(Py, py, float);                                 //! Non-dynamic column with track y-momentum
DECLARE_SOA_COLUMN(Pz, pz, float);                                 //! Non-dynamic column with track z-momentum
DECLARE_SOA_COLUMN(Sign, sign, float);                             //! Non-dynamic column with track sign
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, uint8_t); //!
} // namespace pidtracks
DECLARE_SOA_TABLE(PidTracksReal, "AOD", "PIDTRACKSREAL", //! Real tracks for prediction and domain adaptation
                  aod::track::TPCSignal,
                  aod::pidtofsignal::TOFSignal,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ);
DECLARE_SOA_TABLE(PidTracksMc, "AOD", "PIDTRACKSMC", //! MC tracks for training
                  aod::track::TPCSignal,
                  aod::pidtofsignal::TOFSignal,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ,
                  aod::mcparticle::PdgCode,
                  pidtracks::IsPhysicalPrimary);
} // namespace o2::aod

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill PID train table with MC data."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

struct CreateTableMc {
  Produces<aod::PidTracksMc> pidTracksTable;

  Filter trackFilter = aod::track::isGlobalTrack == (uint8_t) true;
  using BigTracksMC = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksExtended, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;

  void process(BigTracksMC const& tracks, aod::McParticles const& mctracks)
  {
    for (const auto& track : tracks) {
      const auto mcParticle = track.mcParticle();
      uint8_t isPrimary = (uint8_t)mcParticle.isPhysicalPrimary();
      pidTracksTable(track.tpcSignal(), track.tofSignal(),
                     track.px(), track.py(), track.pz(),
                     track.sign(),
                     track.x(), track.y(), track.z(),
                     track.alpha(),
                     track.trackType(),
                     track.tpcNClsShared(),
                     track.dcaXY(), track.dcaZ(),
                     mcParticle.pdgCode(),
                     isPrimary);
    }
  }
};

struct CreateTableReal {
  Produces<aod::PidTracksReal> pidTracksTable;

  Filter trackFilter = aod::track::isGlobalTrack == (uint8_t) true;
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksExtended, aod::TrackSelection, aod::TOFSignal>>;

  void process(BigTracks const& tracks)
  {
    for (const auto& track : tracks) {
      pidTracksTable(track.tpcSignal(), track.tofSignal(),
                     track.px(), track.py(), track.pz(),
                     track.sign(),
                     track.x(), track.y(), track.z(),
                     track.alpha(),
                     track.trackType(),
                     track.tpcNClsShared(),
                     track.dcaXY(), track.dcaZ());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    return WorkflowSpec{
      adaptAnalysisTask<CreateTableMc>(cfgc)};
  } else {
    return WorkflowSpec{
      adaptAnalysisTask<CreateTableReal>(cfgc)};
  }
}
