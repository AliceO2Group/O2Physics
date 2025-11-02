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
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since November 27, 2023

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>

#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullAl, aod::pidTOFFullTr>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullAl, aod::pidTOFFullTr, aod::McTrackLabels>;

struct helium_flow {

  // QC Histograms
  HistogramRegistry registryQC{
    "registryQC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: MC
  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Global Parameters
  Configurable<int> particle_of_interest{"particle_of_interest", 0, "0=antihelium3, 1=antitriton, 2=antihelium4"};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 4, "minimum number of found ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of found TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> min_pt{"min_pt", 1.0f, "minimum pt of the tracks"};
  Configurable<float> min_eta{"min_eta", -0.8f, "minimum_eta"};
  Configurable<float> max_eta{"max_eta", +0.8f, "maximum_eta"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.1f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.1f, "Maximum DCAz"};
  Configurable<float> min_nsigmaTPC{"min_nsigmaTPC", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> max_nsigmaTPC{"max_nsigmaTPC", +3.0f, "Maximum nsigma TPC"};
  Configurable<bool> require_primVtx_contributor{"require_primVtx_contributor", true, "require that the track is a PV contributor"};

  // List of Particles
  enum particle { antihelium3,
                  antitriton,
                  antihelium4 };

  void init(InitContext const&)
  {
    // Global Properties
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{4, 0, 4, "1 = all, 2 = selected, 3 = events with particle of interest"}});
  }

  // Single-Track Selection for the Particle of Interest
  template <typename T1>
  bool passedTrackSelection(const T1& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < min_ITS_nClusters)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < min_TPC_nClusters)
      return false;
    if (track.tpcNClsCrossedRows() < min_TPC_nCrossedRows)
      return false;
    if (track.tpcChi2NCl() > max_chi2_TPC)
      return false;
    if (track.itsChi2NCl() > max_chi2_ITS)
      return false;
    if (track.eta() < min_eta || track.eta() > max_eta)
      return false;
    if (track.pt() < min_pt)
      return false;
    if (TMath::Abs(track.dcaXY()) > max_dcaxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > max_dcaz)
      return false;

    return true;
  }

  template <typename T2>
  bool isParticleOfInterest(const T2& track)
  {
    // Variables
    float nsigmaTPCHe3 = track.tpcNSigmaHe();
    float nsigmaTPCH3 = track.tpcNSigmaTr();
    float nsigmaTOFH3 = track.tofNSigmaTr();
    float nsigmaTPCHe4 = track.tpcNSigmaAl();
    float pt = track.pt();

    // Antimatter Only
    if (track.sign() > 0)
      return false;

    // Antihelium3 ID
    if (particle_of_interest == particle::antihelium3) {
      if (pt >= 1.0 && TMath::Abs(nsigmaTPCHe3) < 8.0)
        return true;
      return false;
    }

    // Antitriton ID
    if (particle_of_interest == particle::antitriton) {
      if (pt >= 1.0 && pt < 2.0 && TMath::Abs(nsigmaTPCH3) < 8.0)
        return true;
      if (pt >= 2.0 && pt < 4.0 && track.hasTOF() && nsigmaTPCH3 > min_nsigmaTPC && nsigmaTPCH3 < max_nsigmaTPC && TMath::Abs(nsigmaTOFH3) < 8.0)
        return true;
      return false;
    }

    // Antihelium4 ID
    if (particle_of_interest == particle::antihelium4) {
      if (pt >= 1.0 && TMath::Abs(nsigmaTPCHe4) < 8.0)
        return true;
      return false;
    }
    return false;
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, FullTracks const& tracks)
  {
    // Event Counter (before event selection)
    registryQC.fill(HIST("number_of_events_data"), 0.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter (after event selection)
    registryQC.fill(HIST("number_of_events_data"), 1.5);

    // Reduced Event
    std::vector<int> particle_ID;
    bool containsParticleOfInterest(false);

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      // Track Selection
      if (!passedTrackSelection(track))
        continue;
      if (!track.passedITSRefit())
        continue;
      if (!track.passedTPCRefit())
        continue;

      // Track Index
      int i = track.globalIndex();

      // Trigger: Particle of Interest
      if (isParticleOfInterest(track))
        containsParticleOfInterest = true;

      // Store Array Element
      particle_ID.push_back(i);
    }

    // Skip Events with no Particle of Interest
    if (!containsParticleOfInterest)
      return;

    // Event Counter (events with particle of interest)
    registryQC.fill(HIST("number_of_events_data"), 3.5);

  } // end processData
  PROCESS_SWITCH(helium_flow, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<helium_flow>(cfgc)};
}
