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
/// \since August 22, 2024

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
#include <TLorentzVector.h>
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
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::McTrackLabels>;

struct nuclei_in_toward_transv_regions {

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
  Configurable<double> min_pt_leading{"min_pt_leading", 5.0, "Minimum pt of leading particle"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 5, "minimum number of ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<double> max_chi2_TPC{"max_chi2_TPC", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> max_chi2_ITS{"max_chi2_ITS", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> min_pt{"min_pt", 0.3, "minimum pt of the tracks"};
  Configurable<double> min_eta{"min_eta", -0.8, "minimum eta"};
  Configurable<double> max_eta{"max_eta", +0.8, "maximum eta"};
  Configurable<double> min_y{"min_y", -0.5, "minimum y"};
  Configurable<double> max_y{"max_y", +0.5, "maximum y"};
  Configurable<double> max_dcaxy{"max_dcaxy", 0.1, "Maximum DCAxy"};
  Configurable<double> max_dcaz{"max_dcaz", 0.1, "Maximum DCAz"};
  Configurable<double> min_nsigmaTPC{"min_nsigmaTPC", -3.0, "Minimum nsigma TPC"};
  Configurable<double> max_nsigmaTPC{"max_nsigmaTPC", +3.0, "Maximum nsigma TPC"};
  Configurable<double> min_nsigmaTOF{"min_nsigmaTOF", -3.0, "Minimum nsigma TOF"};
  Configurable<double> max_nsigmaTOF{"max_nsigmaTOF", +3.5, "Maximum nsigma TOF"};
  Configurable<bool> require_PV_contributor{"require_PV_contributor", true, "require that the track is a PV contributor"};

  void init(InitContext const&)
  {
    // Event Counters
    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{10, 0, 10, "counter"}});
    registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});

    // Data
    registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

    // MC
    registryMC.add("antiproton_prim_jet", "antiproton_prim_jet", HistType::kTH1F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_jet", "antiproton_all_jet", HistType::kTH1F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_prim_ue", "antiproton_prim_ue", HistType::kTH1F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_ue", "antiproton_all_ue", HistType::kTH1F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}});
  }

  // Single-Track Selection for Particles inside Jets
  template <typename JetTrackType>
  bool passedTrackSelectionForJetReconstruction(const JetTrackType& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < 3)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcChi2NCl() > 4)
      return false;
    if (track.itsChi2NCl() > 36)
      return false;
    if (track.eta() < -0.8 || track.eta() > 0.8)
      return false;
    if (track.pt() < 0.1)
      return false;
    if (TMath::Abs(track.dcaXY()) > (0.004f + 0.013f / track.pt()))
      return false;
    if (TMath::Abs(track.dcaZ()) > 2.0)
      return false;

    return true;
  }

  // Single-Track Selection
  template <typename TrackType>
  bool passedTrackSelection(const TrackType& track)
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

    return true;
  }

  // Rapidity
  double get_rapidity(double px, double py, double pz, double mass)
  {
    double rap(0);
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(px, py, pz, mass);
    rap = lorentzVect.Rapidity();
    return rap;
  }

  template <typename T1, typename T2>
  bool isTrackInTowardRegion(const T1& track, const T2& leading_track)
  {
    // Initialization
    bool isInTowardRegion = false;

    // DeltaPhi
    double phi_ref = TVector2::Phi_0_2pi(leading_track.phi());
    double phi_trk = TVector2::Phi_0_2pi(track.phi());
    double delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi_trk - phi_ref);
    if (delta_phi >= 0.0 && delta_phi < 60.0)
      isInTowardRegion = true;
    if (delta_phi >= 300.0 && delta_phi <= 360.0)
      isInTowardRegion = true;

    return isInTowardRegion;
  }

  template <typename T3, typename T4>
  bool isTrackInTransverseRegion(const T3& track, const T4& leading_track)
  {
    // Initialization
    bool isInTransverseRegion = false;

    // DeltaPhi
    double phi_ref = TVector2::Phi_0_2pi(leading_track.phi());
    double phi_trk = TVector2::Phi_0_2pi(track.phi());
    double delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi_trk - phi_ref);
    if (delta_phi >= 60.0 && delta_phi < 120.0)
      isInTransverseRegion = true;
    if (delta_phi >= 240.0 && delta_phi < 300.0)
      isInTransverseRegion = true;

    return isInTransverseRegion;
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, FullTracks const& tracks)
  {
    // Event Counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter: after event selection sel8
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Cut on z-vertex
    if (abs(collision.posZ()) > zVtx)
      return;

    // Event Counter: after z-vertex cut
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // Leading Track
    int leading_ID(0);
    double pt_max(0);

    // Track Index
    int i = -1;

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      i++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      if (track.pt() > pt_max) {
        leading_ID = i;
        pt_max = track.pt();
      }
    }
    // Event Counter: Skip Events with pt<pt_leading_min
    if (pt_max < min_pt_leading)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Momentum of the Leading Particle
    auto const& leading_track = tracks.iteratorAt(leading_ID);

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      // Track Selection
      if (!passedTrackSelection(track))
        continue;
      if (require_PV_contributor && !(track.isPVContributor()))
        continue;
      if (track.sign() > 0)
        continue;
      if (TMath::Abs(track.dcaXY()) > max_dcaxy)
        continue;
      if (TMath::Abs(track.dcaZ()) > max_dcaz)
        continue;

      // Variables
      double nsigmaTPCPr = track.tpcNSigmaPr();
      double nsigmaTOFPr = track.tofNSigmaPr();
      double y_proton = get_rapidity(track.px(), track.py(), track.pz(), 0.93827208816);
      if (y_proton < min_y || y_proton > max_y)
        continue;

      // Jet
      if (isTrackInTowardRegion(track, leading_track)) {
        if (track.pt() < 1.0)
          registryData.fill(HIST("antiproton_jet_tpc"), track.pt(), nsigmaTPCPr);
        if (track.pt() >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
          registryData.fill(HIST("antiproton_jet_tof"), track.pt(), nsigmaTOFPr);
      }

      // UE
      if (isTrackInTransverseRegion(track, leading_track)) {
        if (track.pt() < 1.0)
          registryData.fill(HIST("antiproton_ue_tpc"), track.pt(), nsigmaTPCPr);
        if (track.pt() >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
          registryData.fill(HIST("antiproton_ue_tof"), track.pt(), nsigmaTOFPr);
      }
    }
  }

  Preslice<MCTracks> perCollision = o2::aod::track::collisionId;

  void processSecAntiprotons(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {
    for (const auto& collision : collisions) {

      registryMC.fill(HIST("number_of_events_mc"), 0.5);

      // Event Selection
      if (!collision.sel8())
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 1.5);

      if (abs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 2.5);

      auto tracks_per_coll = mcTracks.sliceBy(perCollision, collision.globalIndex());

      int leading_ID(0);
      double pt_max(0);
      int i = -1;

      // Loop over Reconstructed Tracks
      for (auto track : tracks_per_coll) {

        i++;
        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        if (track.pt() > pt_max) {
          leading_ID = i;
          pt_max = track.pt();
        }
      }
      if (pt_max < min_pt_leading)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 3.5);

      // Momentum of the Leading Particle
      auto const& leading_track = tracks_per_coll.iteratorAt(leading_ID);

      // Loop over Reconstructed Tracks
      for (auto track : tracks_per_coll) {

        if (!passedTrackSelection(track))
          continue;
        if (require_PV_contributor && !(track.isPVContributor()))
          continue;
        if (track.sign() > 0)
          continue;
        if (TMath::Abs(track.dcaXY()) > max_dcaxy)
          continue;
        if (TMath::Abs(track.dcaZ()) > max_dcaz)
          continue;
        double y_proton = get_rapidity(track.px(), track.py(), track.pz(), 0.93827208816);
        if (y_proton < min_y || y_proton > max_y)
          continue;

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;
        const auto particle = track.mcParticle();
        if (particle.pdgCode() != -2212)
          continue;

        // Jet
        if (isTrackInTowardRegion(track, leading_track)) {
          registryMC.fill(HIST("antiproton_all_jet"), track.pt());
          if (particle.isPhysicalPrimary()) {
            registryMC.fill(HIST("antiproton_prim_jet"), track.pt());
          }
        }
        // UE
        if (isTrackInTransverseRegion(track, leading_track)) {
          registryMC.fill(HIST("antiproton_all_ue"), track.pt());
          if (particle.isPhysicalPrimary()) {
            registryMC.fill(HIST("antiproton_prim_ue"), track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_toward_transv_regions, processData, "Process Data", true);
  PROCESS_SWITCH(nuclei_in_toward_transv_regions, processSecAntiprotons, "Process sec antip", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nuclei_in_toward_transv_regions>(cfgc)};
}
