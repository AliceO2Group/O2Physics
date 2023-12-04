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
/// \since November 22, 2023

#include <vector>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TDatabasePDG.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/DataTypes.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::McTrackLabels>;

struct nuclei_in_jets {

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
  Configurable<float> min_pt_leading{"min_pt_leading", 5.0f, "minimum pt of leading particle"};
  Configurable<float> max_jet_radius{"max_jet_radius", 0.1f, "maximum jet resolution parameter R"};
  Configurable<int> particle_of_interest{"particle_of_interest", 0, "0=antiproton, 1=antideuteron, 2=antihelium3"};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 4, "minimum number of found ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of found TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> min_pt{"min_pt", 0.2f, "minimum pt of the tracks"};
  Configurable<float> min_eta{"min_eta", -0.8f, "minimum_eta"};
  Configurable<float> max_eta{"max_eta", +0.8f, "maximum_eta"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.1f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.1f, "Maximum DCAz"};
  Configurable<float> min_nsigmaTPC{"min_nsigmaTPC", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> max_nsigmaTPC{"max_nsigmaTPC", +3.0f, "Maximum nsigma TPC"};
  Configurable<bool> require_primVtx_contributor{"require_primVtx_contributor", true, "require that the track is a PV contributor"};

  // List of Particles
  enum particle { proton,
                  deuteron,
                  helium };

  void init(InitContext const&)
  {
    // Global Properties and QC
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{5, 0, 5, "1 = all events, 2 = selected events, 3 = events with pt>pt_threshold, 4 = events with pt>pt_threshold and particle of interest"}});
    registryQC.add("jet_plus_ue_multiplicity", "jet + underlying-event multiplicity", HistType::kTH1F, {{300, 0, 300, "#it{N}_{ch}"}});
    registryQC.add("jet_multiplicity", "jet multiplicity", HistType::kTH1F, {{300, 0, 300, "#it{N}_{ch}"}});
    registryQC.add("ue_multiplicity", "underlying-event multiplicity", HistType::kTH1F, {{300, 0, 300, "#it{N}_{ch}"}});
    registryQC.add("pt_leading", "pt leading", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("eta_phi_jet", "DeltaEta DeltaPhi jet", HistType::kTH2F, {{500, -0.5, 0.5, "#Delta#eta"}, {500, 0.0, TMath::Pi(), "#Delta#phi"}});
    registryQC.add("eta_phi_ue", "DeltaEta DeltaPhi UE", HistType::kTH2F, {{500, -0.5, 0.5, "#Delta#eta"}, {500, 0.0, TMath::Pi(), "#Delta#phi"}});
    registryQC.add("r_jet", "R jet", HistType::kTH1F, {{400, 0.0, 0.8, "#it{R}"}});
    registryQC.add("r_ue", "R ue", HistType::kTH1F, {{400, 0.0, 0.8, "#it{R}"}});

    // Antiprotons
    registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH3F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH3F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TOF}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH3F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH3F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TOF}"}, {10, 0, 100, "#it{N}_{ch}"}});

    // Antideuterons
    registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH3F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH3F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TOF}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH3F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH3F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TOF}"}, {10, 0, 100, "#it{N}_{ch}"}});

    // Antihelium-3
    registryData.add("antihelium3_jet_tpc", "antihelium3_jet_tpc", HistType::kTH3F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antihelium3_ue_tpc", "antihelium3_ue_tpc", HistType::kTH3F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {10, 0, 100, "#it{N}_{ch}"}});
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

  // Single-Track Selection for Particles inside Jets
  template <typename T2>
  bool passedMinimalTrackSelection(const T2& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < 2)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < 70)
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcChi2NCl() > 4)
      return false;
    if (track.itsChi2NCl() > 36)
      return false;
    if (track.eta() < -0.8 || track.eta() > 0.8)
      return false;
    if (track.pt() < 0.15)
      return false;
    if (TMath::Abs(track.dcaXY()) > 0.2)
      return false;
    if (TMath::Abs(track.dcaZ()) > 0.2)
      return false;

    return true;
  }

  template <typename T3>
  bool isParticleOfInterest(const T3& track)
  {
    // Variables
    float nsigmaTPCPr = track.tpcNSigmaPr();
    float nsigmaTOFPr = track.tofNSigmaPr();
    float nsigmaTPCDe = track.tpcNSigmaDe();
    float nsigmaTOFDe = track.tofNSigmaDe();
    float nsigmaTPCHe = track.tpcNSigmaHe();
    float pt = track.pt();

    // Antimatter Only
    if (track.sign() > 0)
      return false;

    // Proton ID
    if (particle_of_interest == particle::proton) {
      if (pt < 1.0 && TMath::Abs(nsigmaTPCPr) < 8.0)
        return true;
      if (pt >= 1.0 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 10.0)
        return true;
      return false;
    }

    // Deuteron ID
    if (particle_of_interest == particle::deuteron) {
      if (pt < 1.0 && TMath::Abs(nsigmaTPCDe) < 8.0)
        return true;
      if (pt >= 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF() && TMath::Abs(nsigmaTOFDe) < 10.0)
        return true;
      return false;
    }

    // Helium-3 ID
    if (particle_of_interest == particle::helium) {
      if ((0.5 * pt) >= 1.0 && TMath::Abs(nsigmaTPCHe) < 8.0)
        return true;
      return false;
    }

    return false;
  }

  // Minimum
  float Minimum(float x1, float x2)
  {
    float x_min(x1);
    if (x1 < x2)
      x_min = x1;
    if (x1 >= x2)
      x_min = x2;

    return x_min;
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
    int leading_ID;
    bool containsParticleOfInterest(false);
    float pt_max(0);

    // Track Index Initialization
    int i = -1;

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      // Track Index
      i++;

      // Track Selection for Jet
      if (!passedMinimalTrackSelection(track))
        continue;
      if (!track.passedITSRefit())
        continue;
      if (!track.passedTPCRefit())
        continue;

      // Trigger: Particle of Interest
      if (isParticleOfInterest(track))
        containsParticleOfInterest = true;

      // Find pt Leading
      if (track.pt() > pt_max) {
        leading_ID = i;
        pt_max = track.pt();
      }

      // Store Array Element
      particle_ID.push_back(i);
    }

    // Skip Events with no trigger Particle
    if (pt_max == 0)
      return;

    // Histogram with pt_leading
    registryQC.fill(HIST("pt_leading"), pt_max);

    // Number of Stored Particles
    int nParticles = static_cast<int>(particle_ID.size());

    // Selection of Events with pt > pt_leading
    if (nParticles < 2)
      return;
    if (pt_max < min_pt_leading)
      return;

    // Event Counter (after pt > pt_max selection)
    registryQC.fill(HIST("number_of_events_data"), 2.5);

    // Skip Events with no Particle of Interest
    if (!containsParticleOfInterest)
      return;

    // Event Counter (events with pt > pt_max that contain particle of interest)
    registryQC.fill(HIST("number_of_events_data"), 3.5);

    // Momentum of the Leading Particle
    auto const& leading_track = tracks.iteratorAt(leading_ID);
    TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());

    // Instruction to be removed
    registryQC.fill(HIST("number_of_events_data"), 4.5);

    // Array of Particles inside Jet
    std::vector<int> jet_particle_ID;
    jet_particle_ID.push_back(leading_ID);

    // Labels
    Int_t exit(0);
    Int_t nPartAssociated(0);

    // Jet Finder
    do {
      // Initialization
      float distance_jet_min(1e+08);
      float distance_bkg_min(1e+08);
      int label_jet_particle(0);
      int i_jet_particle(0);

      for (int i = 0; i < nParticles; i++) {

        // Skip Leading Particle & Elements already associated to the Jet
        if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
          continue;

        // Get Particle Momentum
        auto stored_track = tracks.iteratorAt(particle_ID[i]);
        TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

        // Variables
        float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
        float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
        float deltaEta = p_particle.Eta() - p_leading.Eta();
        float deltaPhi = TVector2::Phi_0_2pi(p_particle.Phi() - p_leading.Phi());
        float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
        float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

        // Distances
        float distance_jet = min * Delta2 / (max_jet_radius * max_jet_radius);
        float distance_bkg = one_over_pt2_part;

        // Find Minimum Distance Jet
        if (distance_jet < distance_jet_min) {
          distance_jet_min = distance_jet;
          label_jet_particle = particle_ID[i];
          i_jet_particle = i;
        }

        // Find Minimum Distance Bkg
        if (distance_bkg < distance_bkg_min) {
          distance_bkg_min = distance_bkg;
        }
      }

      if (distance_jet_min <= distance_bkg_min) {

        // Add Particle to Jet
        jet_particle_ID.push_back(label_jet_particle);

        // Update Momentum of Leading Particle
        auto jet_track = tracks.iteratorAt(label_jet_particle);
        TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());
        p_leading = p_leading + p_i;

        // Remove Element
        particle_ID[i_jet_particle] = -1;
        nPartAssociated++;
      }

      if (nPartAssociated >= (nParticles - 1))
        exit = 1;
      if (distance_jet_min > distance_bkg_min)
        exit = 2;

    } while (exit == 0);

    // Multiplicity inside Jet + UE
    int nParticlesJetUE = static_cast<int>(jet_particle_ID.size());

    // Fill Jet Multiplicity
    registryQC.fill(HIST("jet_plus_ue_multiplicity"), static_cast<float>(jet_particle_ID.size()));

    // Perpendicular Cones for UE Estimate
    TVector3 z_positive(0.0, 0.0, 1.0);
    TVector3 z_negative(0.0, 0.0, -1.0);
    TVector3 v1 = (z_positive.Cross(p_leading)).Unit();
    TVector3 v2 = (z_negative.Cross(p_leading)).Unit();
    TVector3 v3 = (p_leading.Cross(v1)).Unit();
    TVector3 v4 = (p_leading.Cross(v2)).Unit();

    // Store UE
    std::vector<int> ue_particle_ID;

    for (int i = 0; i < nParticles; i++) {

      // Skip Leading Particle & Elements already associated to the Jet
      if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
        continue;

      // Get UE Track
      const auto& ue_track = tracks.iteratorAt(particle_ID[i]);

      // Variables
      float deltaEta1 = ue_track.eta() - v1.Eta();
      float deltaEta2 = ue_track.eta() - v2.Eta();
      float deltaEta3 = ue_track.eta() - v3.Eta();
      float deltaEta4 = ue_track.eta() - v4.Eta();

      float deltaPhi1 = TVector2::Phi_0_2pi(ue_track.phi() - v1.Phi());
      float deltaPhi2 = TVector2::Phi_0_2pi(ue_track.phi() - v2.Phi());
      float deltaPhi3 = TVector2::Phi_0_2pi(ue_track.phi() - v3.Phi());
      float deltaPhi4 = TVector2::Phi_0_2pi(ue_track.phi() - v4.Phi());

      float dr1 = TMath::Sqrt(deltaEta1 * deltaEta1 + deltaPhi1 * deltaPhi1);
      float dr2 = TMath::Sqrt(deltaEta2 * deltaEta2 + deltaPhi2 * deltaPhi2);
      float dr3 = TMath::Sqrt(deltaEta3 * deltaEta3 + deltaPhi3 * deltaPhi3);
      float dr4 = TMath::Sqrt(deltaEta4 * deltaEta4 + deltaPhi4 * deltaPhi4);

      registryQC.fill(HIST("eta_phi_ue"), deltaEta1, deltaPhi1);
      registryQC.fill(HIST("eta_phi_ue"), deltaEta2, deltaPhi2);
      registryQC.fill(HIST("eta_phi_ue"), deltaEta3, deltaPhi3);
      registryQC.fill(HIST("eta_phi_ue"), deltaEta4, deltaPhi4);
      registryQC.fill(HIST("r_ue"), TMath::Sqrt(deltaEta1 * deltaEta1 + deltaPhi1 * deltaPhi1));
      registryQC.fill(HIST("r_ue"), TMath::Sqrt(deltaEta2 * deltaEta2 + deltaPhi2 * deltaPhi2));
      registryQC.fill(HIST("r_ue"), TMath::Sqrt(deltaEta3 * deltaEta3 + deltaPhi3 * deltaPhi3));
      registryQC.fill(HIST("r_ue"), TMath::Sqrt(deltaEta4 * deltaEta4 + deltaPhi4 * deltaPhi4));

      // Store Particles in the UE
      if (dr1 < max_jet_radius || dr2 < max_jet_radius || dr3 < max_jet_radius || dr4 < max_jet_radius) {
        ue_particle_ID.push_back(particle_ID[i]);
      }
    }

    // UE Multiplicity
    int nParticlesUE = static_cast<int>(ue_particle_ID.size());

    registryQC.fill(HIST("ue_multiplicity"), static_cast<float>(ue_particle_ID.size()) / 4.0);

    // Jet Multiplicity
    float jet_Nch = static_cast<float>(jet_particle_ID.size()) - static_cast<float>(ue_particle_ID.size()) / 4.0;
    registryQC.fill(HIST("jet_multiplicity"), jet_Nch);

    // Loop over particles inside Jet
    for (int i = 0; i < nParticlesJetUE; i++) {

      const auto& jet_track = tracks.iteratorAt(jet_particle_ID[i]);
      TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());

      float deltaEta = p_i.Eta() - p_leading.Eta();
      float deltaPhi = TVector2::Phi_0_2pi(p_i.Phi() - p_leading.Phi());
      registryQC.fill(HIST("eta_phi_jet"), deltaEta, deltaPhi);
      registryQC.fill(HIST("r_jet"), TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi));

      // Track Selection
      if (!passedTrackSelection(jet_track))
        continue;
      if (require_primVtx_contributor && !(jet_track.isPVContributor()))
        continue;

      // Variables
      float nsigmaTPCPr = jet_track.tpcNSigmaPr();
      float nsigmaTOFPr = jet_track.tofNSigmaPr();
      float nsigmaTPCDe = jet_track.tpcNSigmaDe();
      float nsigmaTOFDe = jet_track.tofNSigmaDe();
      float nsigmaTPCHe = jet_track.tpcNSigmaHe();
      float pt = jet_track.pt();

      // Antiproton
      if (particle_of_interest == particle::proton) {
        if (pt < 1.0)
          registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr, jet_Nch);
        if (pt >= 1.0 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && jet_track.hasTOF())
          registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr, jet_Nch);
      }

      // Antideuteron
      if (particle_of_interest == particle::deuteron) {
        if (pt < 1.0)
          registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe, jet_Nch);
        if (pt >= 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && jet_track.hasTOF())
          registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe, jet_Nch);
      }

      // Antihelium3
      if (particle_of_interest == particle::helium) {
        registryData.fill(HIST("antihelium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe, jet_Nch);
      }
    }

    // Loop over particles inside UE
    for (int i = 0; i < nParticlesUE; i++) {

      const auto& ue_track = tracks.iteratorAt(ue_particle_ID[i]);

      // Track Selection
      if (!passedTrackSelection(ue_track))
        continue;
      if (require_primVtx_contributor && !(ue_track.isPVContributor()))
        continue;

      // Variables
      float nsigmaTPCPr = ue_track.tpcNSigmaPr();
      float nsigmaTOFPr = ue_track.tofNSigmaPr();
      float nsigmaTPCDe = ue_track.tpcNSigmaDe();
      float nsigmaTOFDe = ue_track.tofNSigmaDe();
      float nsigmaTPCHe = ue_track.tpcNSigmaHe();
      float pt = ue_track.pt();

      // Antiproton
      if (particle_of_interest == particle::proton) {
        if (pt < 1.0)
          registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr, jet_Nch);
        if (pt >= 1.0 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && ue_track.hasTOF())
          registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr, jet_Nch);
      }

      // Antideuteron
      if (particle_of_interest == particle::deuteron) {
        if (pt < 1.0)
          registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe, jet_Nch);
        if (pt >= 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && ue_track.hasTOF())
          registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe, jet_Nch);
      }

      // Antihelium3
      if (particle_of_interest == particle::helium) {
        registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe, jet_Nch);
      }
    }

  } // end processData
  PROCESS_SWITCH(nuclei_in_jets, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nuclei_in_jets>(cfgc)};
}
