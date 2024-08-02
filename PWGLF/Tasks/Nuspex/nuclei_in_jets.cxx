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
#include <TLorentzVector.h>
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
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

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
  Configurable<double> min_pt_leading{"min_pt_leading", 5.0, "Minimum pt of leading particle"};
  Configurable<double> min_jet_pt{"min_jet_pt", 10.0, "Minimum pt of the jet"};
  Configurable<double> Rjet{"Rjet", 0.3, "Jet resolution parameter R"};
  Configurable<double> Rmax{"Rmax", 0.3, "Maximum radius for jet and UE regions"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<int> min_nPartInJet{"min_nPartInJet", 2, "Minimum number of particles inside jet"};

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
    // QC Histograms
    registryQC.add("multiplicityJetPlusUE", "multiplicityJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("multiplicityJet", "multiplicityJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("multiplicityUE", "multiplicityUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("ptLeading", "ptLeading", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("etaLeading", "etaLeading", HistType::kTH1F, {{100, -0.8, 0.8, "#eta"}});
    registryQC.add("phiLeading", "phiLeading", HistType::kTH1F, {{100, 0, TMath::TwoPi(), "#phi"}});
    registryQC.add("rJet", "rJet", HistType::kTH1F, {{100, 0.0, 0.5, "#it{R}"}});
    registryQC.add("rUE", "rUE", HistType::kTH1F, {{100, 0.0, 0.5, "#it{R}"}});
    registryQC.add("angle_jet_leading_track", "angle_jet_leading_track", HistType::kTH1F, {{200, 0.0, 50.0, "#theta"}});
    registryQC.add("ptJetPlusUE", "ptJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("ptJet", "ptJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("ptUE", "ptUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("deltaEtadeltaPhiJet", "deltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhiUE", "deltaEtadeltaPhiUE", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhi_leading_jet", "deltaEtadeltaPhi_leading_jet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});

    // Event Counters
    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{10, 0, 10, "counter"}});
    registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});

    // Antiprotons
    registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH2F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH2F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH2F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH2F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antiproton_dca_jet", "antiproton_dca_jet", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});
    registryData.add("antiproton_dca_ue", "antiproton_dca_ue", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});

    // Antideuterons
    registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH2F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH2F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH2F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH2F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

    // Antihelium-3
    registryData.add("antihelium3_jet_tpc", "antihelium3_jet_tpc", HistType::kTH2F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antihelium3_ue_tpc", "antihelium3_ue_tpc", HistType::kTH2F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});

    // Generated
    registryMC.add("antiproton_jet_gen", "antiproton_jet_gen", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_gen", "antideuteron_jet_gen", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_jet_gen", "antihelium3_jet_gen", HistType::kTH1F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_gen", "antiproton_ue_gen", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_gen", "antideuteron_ue_gen", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_ue_gen", "antihelium3_ue_gen", HistType::kTH1F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Reconstructed TPC
    registryMC.add("antiproton_jet_rec_tpc", "antiproton_jet_rec_tpc", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_rec_tpc", "antideuteron_jet_rec_tpc", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_jet_rec_tpc", "antihelium3_jet_rec_tpc", HistType::kTH1F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_rec_tpc", "antiproton_ue_rec_tpc", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_rec_tpc", "antideuteron_ue_rec_tpc", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_ue_rec_tpc", "antihelium3_ue_rec_tpc", HistType::kTH1F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Reconstructed TOF
    registryMC.add("antiproton_jet_rec_tof", "antiproton_jet_rec_tof", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_rec_tof", "antideuteron_jet_rec_tof", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_rec_tof", "antiproton_ue_rec_tof", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_rec_tof", "antideuteron_ue_rec_tof", HistType::kTH1F, {{50, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});

    // DCA Templates
    registryMC.add("antiproton_dca_prim", "antiproton_dca_prim", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});
    registryMC.add("antiproton_dca_sec", "antiproton_dca_sec", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});

    // Fraction of Primary Antiprotons from MC
    registryMC.add("antiproton_prim", "antiproton_prim", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all", "antiproton_all", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_prim_jet", "antiproton_prim_jet", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_jet", "antiproton_all_jet", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_prim_ue", "antiproton_prim_ue", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_ue", "antiproton_all_ue", HistType::kTH1F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
  }

  // Single-Track Selection for Particles inside Jets
  template <typename T1>
  bool passedTrackSelectionForJetReconstruction(const T1& track)
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
    if (TMath::Abs(track.dcaXY()) > (0.0105 * 0.035 / TMath::Power(track.pt(), 1.1)))
      return false;
    if (TMath::Abs(track.dcaZ()) > 2.0)
      return false;

    return true;
  }

  // Single-Track Selection
  template <typename T2>
  bool passedTrackSelection(const T2& track)
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

  template <typename T3>
  bool isHighPurityAntiproton(const T3& track)
  {
    // Variables
    double nsigmaTPCPr = track.tpcNSigmaPr();
    double nsigmaTOFPr = track.tofNSigmaPr();
    double pt = track.pt();

    if (pt < 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0)
      return true;
    if (pt >= 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0 && track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 2.0)
      return true;
    return false;
  }

  // Minimum
  double Minimum(double x1, double x2)
  {
    double x_min(x1);
    if (x1 < x2)
      x_min = x1;
    if (x1 >= x2)
      x_min = x2;

    return x_min;
  }

  // Deltaphi
  double GetDeltaPhi(double a1, double a2)
  {
    double delta_phi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = TMath::Abs(phi1 - phi2);

    if (diff <= TMath::Pi())
      delta_phi = diff;
    if (diff > TMath::Pi())
      delta_phi = TMath::TwoPi() - diff;

    return delta_phi;
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

  void get_perpendicular_axis(TVector3 p, TVector3& u, double sign)
  {
    // Initialization
    double ux(0), uy(0), uz(0);

    // Components of Vector p
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    // Protection 1
    if (px == 0 && py != 0) {

      uy = -(pz * pz) / py;
      ux = sign * sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Protection 2
    if (py == 0 && px != 0) {

      ux = -(pz * pz) / px;
      uy = sign * sqrt(px * px - (pz * pz * pz * pz) / (px * px));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Equation Parameters
    double a = px * px + py * py;
    double b = 2.0 * px * pz * pz;
    double c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    double delta = b * b - 4.0 * a * c;

    // Protection agains delta<0
    if (delta < 0) {
      return;
    }

    // Solutions
    ux = (-b + sign * sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
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

    // Indices of Reduced Event
    std::vector<int> particle_ID;

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
      particle_ID.push_back(i);
    }

    // Momentum of the Leading Particle
    auto const& leading_track = tracks.iteratorAt(leading_ID);
    TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());
    TVector3 p_jet(leading_track.px(), leading_track.py(), leading_track.pz());

    // QC: pt, eta, and phi Distributions of Leading Track
    registryQC.fill(HIST("ptLeading"), p_leading.Pt());
    registryQC.fill(HIST("etaLeading"), p_leading.Eta());
    registryQC.fill(HIST("phiLeading"), p_leading.Phi());

    // Event Counter: Skip Events with pt<pt_leading_min
    if (pt_max < min_pt_leading)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Labels
    int exit(0);
    int nPartAssociated(0);
    int nParticles = static_cast<int>(particle_ID.size());
    std::vector<int> jet_particle_ID;
    jet_particle_ID.push_back(leading_ID);

    // Jet Finder
    do {
      // Initialization
      double distance_jet_min(1e+08);
      double distance_bkg_min(1e+08);
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
        double one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
        double one_over_pt2_lead = 1.0 / (p_jet.Pt() * p_jet.Pt());
        double deltaEta = p_particle.Eta() - p_jet.Eta();
        double deltaPhi = GetDeltaPhi(p_particle.Phi(), p_jet.Phi());
        double min = Minimum(one_over_pt2_part, one_over_pt2_lead);
        double Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

        // Distances
        double distance_jet = min * Delta2 / (Rjet * Rjet);
        double distance_bkg = one_over_pt2_part;

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
        p_jet = p_jet + p_i;

        // Remove Element
        particle_ID[i_jet_particle] = -1;
        nPartAssociated++;
      }

      if (nPartAssociated >= (nParticles - 1))
        exit = 1;
      if (distance_jet_min > distance_bkg_min)
        exit = 2;

    } while (exit == 0);

    // Event Counter: Skip Events with jet not fully inside acceptance
    if ((TMath::Abs(p_jet.Eta()) + Rmax) > max_eta)
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Perpendicular Cones for the UE Estimate
    TVector3 ue_axis1(0.0, 0.0, 0.0);
    TVector3 ue_axis2(0.0, 0.0, 0.0);
    get_perpendicular_axis(p_jet, ue_axis1, +1.0);
    get_perpendicular_axis(p_jet, ue_axis2, -1.0);

    // Protection against delta<0
    if (ue_axis1.Mag() == 0 || ue_axis2.Mag() == 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    // QA Plots
    double deltaEta = p_leading.Eta() - p_jet.Eta();
    double deltaPhi = GetDeltaPhi(p_leading.Phi(), p_jet.Phi());
    registryQC.fill(HIST("angle_jet_leading_track"), (180.0 / TMath::Pi()) * p_leading.Angle(p_jet));
    registryQC.fill(HIST("deltaEtadeltaPhi_leading_jet"), deltaEta, deltaPhi);

    double NchJetPlusUE(0);
    double NchJet(0);
    double NchUE(0);
    double ptJetPlusUE(0);
    double ptJet(0);
    double ptUE(0);

    for (int i = 0; i < nParticles; i++) {

      auto track = tracks.iteratorAt(particle_ID[i]);

      TVector3 particle_dir(track.px(), track.py(), track.pz());
      double deltaEta_jet = particle_dir.Eta() - p_jet.Eta();
      double deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), p_jet.Phi());
      double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
      double deltaEta_ue1 = particle_dir.Eta() - ue_axis1.Eta();
      double deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue_axis1.Phi());
      double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
      double deltaEta_ue2 = particle_dir.Eta() - ue_axis2.Eta();
      double deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue_axis2.Phi());
      double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

      if (deltaR_jet < Rmax) {
        registryQC.fill(HIST("deltaEtadeltaPhiJet"), deltaEta_jet, deltaPhi_jet);
        registryQC.fill(HIST("rJet"), deltaR_jet);
        NchJetPlusUE++;
        ptJetPlusUE = ptJetPlusUE + track.pt();
      }
      if (deltaR_ue1 < Rmax) {
        registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEta_ue1, deltaPhi_ue1);
        registryQC.fill(HIST("rUE"), deltaR_ue1);
        NchUE++;
        ptUE = ptUE + track.pt();
      }
      if (deltaR_ue2 < Rmax) {
        registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEta_ue2, deltaPhi_ue2);
        registryQC.fill(HIST("rUE"), deltaR_ue2);
        NchUE++;
        ptUE = ptUE + track.pt();
      }
    }

    NchJet = NchJetPlusUE - 0.5 * NchUE;
    ptJet = ptJetPlusUE - 0.5 * ptUE;
    registryQC.fill(HIST("multiplicityJetPlusUE"), NchJetPlusUE);
    registryQC.fill(HIST("multiplicityJet"), NchJet);
    registryQC.fill(HIST("multiplicityUE"), 0.5 * NchUE);
    registryQC.fill(HIST("ptJetPlusUE"), ptJetPlusUE);
    registryQC.fill(HIST("ptJet"), jetPt);
    registryQC.fill(HIST("ptUE"), 0.5 * ptUE);

    // Event Counter: Skip Events with n. particles in jet less than given value
    if (NchJetPlusUE < min_nPartInJet)
      return;
    registryData.fill(HIST("number_of_events_data"), 6.5);

    // Event Counter: Skip Events with Jet Pt lower than threshold
    if (ptJet < min_jet_pt)
      return;
    registryData.fill(HIST("number_of_events_data"), 7.5);

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      if (!passedTrackSelection(track))
        continue;
      if (require_PV_contributor && !(track.isPVContributor()))
        continue;
      if (track.sign() > 0)
        continue;

      // Variables
      double nsigmaTPCPr = track.tpcNSigmaPr();
      double nsigmaTOFPr = track.tofNSigmaPr();
      double nsigmaTPCDe = track.tpcNSigmaDe();
      double nsigmaTOFDe = track.tofNSigmaDe();
      double nsigmaTPCHe = track.tpcNSigmaHe();
      double pt = track.pt();
      double dcaxy = track.dcaXY();
      double dcaz = track.dcaZ();
      double y_proton = get_rapidity(track.px(), track.py(), track.pz(), 0.93827208816);
      double y_deuteron = get_rapidity(track.px(), track.py(), track.pz(), 1.87561294257);
      double y_helium3 = get_rapidity(2.0 * track.px(), 2.0 * track.py(), 2.0 * track.pz(), 2.80839160743);

      TVector3 particle_dir(track.px(), track.py(), track.pz());
      double deltaEta_jet = particle_dir.Eta() - p_jet.Eta();
      double deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), p_jet.Phi());
      double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
      double deltaEta_ue1 = particle_dir.Eta() - ue_axis1.Eta();
      double deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue_axis1.Phi());
      double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
      double deltaEta_ue2 = particle_dir.Eta() - ue_axis2.Eta();
      double deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue_axis2.Phi());
      double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

      // DCAxy Distributions of Antiprotons
      if (isHighPurityAntiproton(track) && TMath::Abs(dcaz) < max_dcaz) {
        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("antiproton_dca_jet"), pt, dcaxy);
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("antiproton_dca_ue"), pt, dcaxy);
        }
      }

      // DCA Cuts
      if (TMath::Abs(dcaxy) > max_dcaxy)
        continue;
      if (TMath::Abs(dcaz) > max_dcaz)
        continue;

      // Jet
      if (deltaR_jet < Rmax) {

        // Antiproton
        if (y_proton > min_y && y_proton < max_y) {
          if (pt < 1.0)
            registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr);
          if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
            registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr);
        }

        // Antideuteron
        if (y_deuteron > min_y && y_deuteron < max_y) {
          if (pt < 1.0)
            registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe);
          if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
            registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe);
        }

        // Antihelium3
        if (y_helium3 > min_y && y_helium3 < max_y) {
          registryData.fill(HIST("antihelium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
        }
      }

      // UE
      if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {

        // Antiproton
        if (y_proton > min_y && y_proton < max_y) {
          if (pt < 1.0)
            registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr);
          if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
            registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr);
        }

        // Antideuteron
        if (y_deuteron > min_y && y_deuteron < max_y) {
          if (pt < 1.0)
            registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe);
          if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
            registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe);
        }

        // Antihelium3
        if (y_helium3 > min_y && y_helium3 < max_y) {
          registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
        }
      }
    }
  }

  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<MCTracks> perCollision = o2::aod::track::collisionId;

  void processMC(o2::aod::McCollisions const& mcCollisions, SimCollisions const& collisions, MCTracks const& mcTracks, aod::McParticles const& mcParticles)
  {
    // Generated Events
    for (const auto& mccollision : mcCollisions) {

      registryMC.fill(HIST("number_of_events_mc"), 0.5);
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      for (auto& particle : mcParticles_per_coll) {

        if (!particle.isPhysicalPrimary())
          continue;
        if ((particle.pdgCode() != -2212) && (particle.pdgCode() != -1000010020) && (particle.pdgCode() != -1000020030))
          continue;
        if (particle.y() < min_y || particle.y() > max_y)
          continue;

        if (particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_jet_gen"), particle.pt());
          registryMC.fill(HIST("antiproton_ue_gen"), particle.pt());
        }
        if (particle.pdgCode() == -1000010020) {
          registryMC.fill(HIST("antideuteron_jet_gen"), particle.pt());
          registryMC.fill(HIST("antideuteron_ue_gen"), particle.pt());
        }
        if (particle.pdgCode() == -1000020030) {
          registryMC.fill(HIST("antihelium3_jet_gen"), particle.pt());
          registryMC.fill(HIST("antihelium3_ue_gen"), particle.pt());
        }
      }
    }

    // Reconstructed Events
    for (const auto& collision : collisions) {

      registryMC.fill(HIST("number_of_events_mc"), 1.5);

      // Event Selection
      if (!collision.sel8())
        continue;

      if (abs(collision.posZ()) > 10)
        continue;

      // Event Counter (after event sel)
      registryMC.fill(HIST("number_of_events_mc"), 2.5);

      auto tracks_per_coll = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // Reconstructed Tracks
      for (auto track : tracks_per_coll) {

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;

        const auto particle = track.mcParticle();
        if ((particle.pdgCode() != -2212) && (particle.pdgCode() != -1000010020) && (particle.pdgCode() != -1000020030))
          continue;

        // Track Selection
        if (!passedTrackSelection(track))
          continue;
        if (require_PV_contributor && !(track.isPVContributor()))
          continue;

        // Variables
        float nsigmaTPCPr = track.tpcNSigmaPr();
        float nsigmaTOFPr = track.tofNSigmaPr();
        float nsigmaTPCDe = track.tpcNSigmaDe();
        float nsigmaTOFDe = track.tofNSigmaDe();
        float nsigmaTPCHe = track.tpcNSigmaHe();
        float pt = track.pt();

        // DCA Templates
        if (particle.pdgCode() == -2212 && particle.isPhysicalPrimary() && TMath::Abs(track.dcaZ()) < max_dcaz)
          registryMC.fill(HIST("antiproton_dca_prim"), pt, track.dcaXY());

        if (particle.pdgCode() == -2212 && (!particle.isPhysicalPrimary()) && TMath::Abs(track.dcaZ()) < max_dcaz)
          registryMC.fill(HIST("antiproton_dca_sec"), pt, track.dcaXY());

        // DCA Cuts
        if (TMath::Abs(track.dcaXY()) > max_dcaxy)
          continue;
        if (TMath::Abs(track.dcaZ()) > max_dcaz)
          continue;

        // Fraction of Primary Antiprotons
        if (particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_all"), pt);
          if (particle.isPhysicalPrimary()) {
            registryMC.fill(HIST("antiproton_prim"), pt);
          }
        }

        if (!particle.isPhysicalPrimary())
          continue;

        // Antiproton
        if (particle.pdgCode() == -2212) {
          if (pt < 1.0 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC) {
            registryMC.fill(HIST("antiproton_jet_rec_tpc"), pt);
            registryMC.fill(HIST("antiproton_ue_rec_tpc"), pt);
          }
          if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF() && nsigmaTOFPr > min_nsigmaTOF && nsigmaTOFPr < max_nsigmaTOF) {
            registryMC.fill(HIST("antiproton_jet_rec_tof"), pt);
            registryMC.fill(HIST("antiproton_ue_rec_tof"), pt);
          }
        }

        // Antideuteron
        if (particle.pdgCode() == -1000010020) {
          if (pt < 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC) {
            registryMC.fill(HIST("antideuteron_jet_rec_tpc"), pt);
            registryMC.fill(HIST("antideuteron_ue_rec_tpc"), pt);
          }
          if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF() && nsigmaTOFDe > min_nsigmaTOF && nsigmaTOFDe < max_nsigmaTOF) {
            registryMC.fill(HIST("antideuteron_jet_rec_tof"), pt);
            registryMC.fill(HIST("antideuteron_ue_rec_tof"), pt);
          }
        }

        // Antihelium-3
        if (particle.pdgCode() == -1000020030) {
          if (nsigmaTPCHe > min_nsigmaTPC && nsigmaTPCHe < max_nsigmaTPC) {
            registryMC.fill(HIST("antihelium3_jet_rec_tpc"), 2.0 * pt);
            registryMC.fill(HIST("antihelium3_ue_rec_tpc"), 2.0 * pt);
          }
        }
      }
    }
  }

  void processSecAntiprotons(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {
    for (const auto& collision : collisions) {

      // Event Selection
      if (!collision.sel8())
        continue;

      if (abs(collision.posZ()) > zVtx)
        continue;

      auto tracks_per_coll = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // Reduced Event
      std::vector<int> particle_ID;
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
        particle_ID.push_back(i);
      }
      if (pt_max < min_pt_leading)
        continue;

      // Number of Stored Particles
      int nParticles = static_cast<int>(particle_ID.size());

      // Momentum of the Leading Particle
      auto const& leading_track = tracks_per_coll.iteratorAt(leading_ID);
      TVector3 p_jet(leading_track.px(), leading_track.py(), leading_track.pz());

      // Labels
      Int_t exit(0);
      Int_t nPartAssociated(0);
      std::vector<int> jet_particle_ID;
      Double_t jetPt(0);

      // Jet Finder
      do {
        // Initialization
        double distance_jet_min(1e+08);
        double distance_bkg_min(1e+08);
        int label_jet_particle(0);
        int i_jet_particle(0);

        for (int i = 0; i < nParticles; i++) {

          // Skip Leading Particle & Elements already associated to the Jet
          if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
            continue;

          // Get Particle Momentum
          auto stored_track = tracks_per_coll.iteratorAt(particle_ID[i]);
          TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

          // Variables
          double one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
          double one_over_pt2_lead = 1.0 / (p_jet.Pt() * p_jet.Pt());
          double deltaEta = p_particle.Eta() - p_jet.Eta();
          double deltaPhi = GetDeltaPhi(p_particle.Phi(), p_jet.Phi());
          double min = Minimum(one_over_pt2_part, one_over_pt2_lead);
          double Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

          // Distances
          double distance_jet = min * Delta2 / (Rjet * Rjet);
          double distance_bkg = one_over_pt2_part;

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
          auto jet_track = tracks_per_coll.iteratorAt(label_jet_particle);
          TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());
          p_jet = p_jet + p_i;
          jetPt = jetPt + jet_track.pt();

          // Remove Element
          particle_ID[i_jet_particle] = -1;
          nPartAssociated++;
        }

        if (nPartAssociated >= (nParticles - 1))
          exit = 1;
        if (distance_jet_min > distance_bkg_min)
          exit = 2;

      } while (exit == 0);

      // Cut on Jet Pt
      if (jetPt < min_jet_pt)
        return;

      // Skip Events with n. particles in jet less than given value
      int nClusteredParticles = static_cast<int>(jet_particle_ID.size());
      if (nClusteredParticles < min_nPartInJet)
        return;

      if ((TMath::Abs(p_jet.Eta()) + Rmax) > max_eta)
        continue;

      // Perpendicular Cones for UE Estimate
      TVector3 ue_axis1(0.0, 0.0, 0.0);
      TVector3 ue_axis2(0.0, 0.0, 0.0);
      get_perpendicular_axis(p_jet, ue_axis1, +1.0);
      get_perpendicular_axis(p_jet, ue_axis2, -1.0);

      // Protection against delta<0
      if (ue_axis1.Mag() == 0 || ue_axis2.Mag() == 0)
        continue;

      for (auto track : tracks_per_coll) {

        if (!passedTrackSelection(track))
          continue;
        if (require_PV_contributor && !(track.isPVContributor()))
          continue;
        if (TMath::Abs(track.dcaXY()) > max_dcaxy)
          continue;
        if (TMath::Abs(track.dcaZ()) > max_dcaz)
          continue;

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;
        const auto particle = track.mcParticle();
        if (particle.pdgCode() != -2212)
          continue;

        TVector3 particle_dir(track.px(), track.py(), track.pz());
        float deltaEta_jet = particle_dir.Eta() - p_jet.Eta();
        float deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), p_jet.Phi());
        float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        float deltaEta_ue1 = particle_dir.Eta() - ue_axis1.Eta();
        float deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue_axis1.Phi());
        float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        float deltaEta_ue2 = particle_dir.Eta() - ue_axis2.Eta();
        float deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue_axis2.Phi());
        float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

        if (deltaR_jet < Rmax) {
          registryMC.fill(HIST("antiproton_all_jet"), track.pt());
          if (particle.isPhysicalPrimary())
            registryMC.fill(HIST("antiproton_prim_jet"), track.pt());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryMC.fill(HIST("antiproton_all_ue"), track.pt());
          if (particle.isPhysicalPrimary())
            registryMC.fill(HIST("antiproton_prim_ue"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_jets, processData, "Process Data", true);
  PROCESS_SWITCH(nuclei_in_jets, processMC, "process MC", false);
  PROCESS_SWITCH(nuclei_in_jets, processSecAntiprotons, "Process sec antip", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nuclei_in_jets>(cfgc)};
}
