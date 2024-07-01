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
  Configurable<float> min_pt_leading{"min_pt_leading", 5.0f, "minimum pt of leading particle"};
  Configurable<float> Rparameter_jet{"Rparameter_jet", 0.3f, "jet resolution parameter R"};
  Configurable<float> Rmax_jet_ue{"Rmax_jet_ue", 0.3f, "Maximum radius"};
  Configurable<int> particle_of_interest{"particle_of_interest", 0, "0=antiproton, 1=antideuteron, 2=antihelium3"};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 4, "minimum number of found ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of found TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> min_pt{"min_pt", 0.2f, "minimum pt of the tracks"};
  Configurable<float> min_eta{"min_eta", -0.8f, "minimum eta"};
  Configurable<float> max_eta{"max_eta", +0.8f, "maximum eta"};
  Configurable<float> min_y{"min_y", -0.5f, "minimum y"};
  Configurable<float> max_y{"max_y", +0.5f, "maximum y"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.1f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.1f, "Maximum DCAz"};
  Configurable<float> min_nsigmaTPC{"min_nsigmaTPC", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> max_nsigmaTPC{"max_nsigmaTPC", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> min_nsigmaTOF{"min_nsigmaTOF", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> max_nsigmaTOF{"max_nsigmaTOF", +3.5f, "Maximum nsigma TOF"};
  Configurable<float> coalMom{"coalMom", +0.3f, "Coalescence Momentum"};
  Configurable<bool> require_primVtx_contributor{"require_primVtx_contributor", true, "require that the track is a PV contributor"};
  Configurable<bool> applyReweighting{"applyReweighting", true, "apply Reweighting"};
  Configurable<std::vector<float>> param_proton_ref{"param_proton_ref", {0.000000000151, 984.796940000000, 0.458560000000, 0.000360000000}, "Parameters of Levi-Tsallis fit of Protons from pythia"};
  Configurable<std::vector<float>> param_proton_jet{"param_proton_jet", {0.025285984480, 4.271090000000, 0.511730000000, 2.261760000000}, "Parameters for reweighting protons in jets"};
  Configurable<std::vector<float>> param_deuteron_jet{"param_deuteron_jet", {20.121453090900, 2.096620000000, 1.000000000000, 1.875610000000}, "Parameters for reweighting deuterons in jets"};
  Configurable<std::vector<float>> param_helium3_jet{"param_helium3_jet", {0.00026, 2.09662, 1.00000, 1.87561}, "Parameters for reweighting helium3 in jets"};
  Configurable<std::vector<float>> param_proton_ue{"param_proton_ue", {0.000000977546, 88.480170000000, 0.539520000000, 0.062120000000}, "Parameters for reweighting protons in ue"};
  Configurable<std::vector<float>> param_deuteron_ue{"param_deuteron_ue", {0.000294563800, 25.000000000000, 0.514970000000, 1.000000000000}, "Parameters for reweighting deuterons in ue"};
  Configurable<std::vector<float>> param_helium3_ue{"param_helium3_ue", {0.00000, 25.00000, 0.51497, 1.00000}, "Parameters for reweighting helium3 in ue"};

  // List of Particles
  enum nucleus { proton,
                 deuteron,
                 helium };

  enum region { jet,
                underlying_event };

  void init(InitContext const&)
  {
    // Global Properties and QC
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{15, 0, 15, "counter"}});
    registryQC.add("number_of_events_mc", "number of events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});
    registryQC.add("event_counter_mc", "event_counter_mc", HistType::kTH1F, {{10, 0, 10, "counter"}});
    registryQC.add("jet_plus_ue_multiplicity", "jet + underlying-event multiplicity", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("jet_multiplicity", "jet multiplicity", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("ue_multiplicity", "underlying-event multiplicity", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("pt_leading", "pt leading", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("eta_phi_jet", "DeltaEta DeltaPhi jet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("eta_phi_ue", "DeltaEta DeltaPhi UE", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("r_max_jet", "R Max jet", HistType::kTH1F, {{200, 0.0, 6.0, "#it{R}_{max}"}});
    registryQC.add("r_jet", "R jet", HistType::kTH1F, {{200, 0.0, 1.0, "#it{R}"}});
    registryQC.add("r_ue", "R ue", HistType::kTH1F, {{200, 0.0, 1.0, "#it{R}"}});
    registryQC.add("eta_leading", "eta_leading", HistType::kTH1F, {{100, -1, 1, "#eta"}});
    registryQC.add("phi_leading", "phi_leading", HistType::kTH1F, {{100, -TMath::Pi(), TMath::Pi(), "#phi"}});
    registryQC.add("angle_jet_leading_track", "angle_jet_leading_track", HistType::kTH1F, {{200, 0.0, TMath::Pi() / 4.0, "#theta"}});
    registryQC.add("deltaPt", "deltaPt", HistType::kTH1F, {{200, -2, 2, "#Delta p_{T}"}});

    // Antiprotons
    registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH3F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH3F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH3F, {{20, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH3F, {{90, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}, {2, 0, 100, "#it{N}_{ch}"}});

    // Eta - pt Maps Antiprotons
    registryData.add("antip_tpc_jet_eta_pt", "antip_tpc_jet_eta_pt", HistType::kTH2F, {{60, 0.0, 3.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("antip_tpc_ue_eta_pt", "antip_tpc_ue_eta_pt", HistType::kTH2F, {{60, 0.0, 3.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("antip_tof_jet_eta_pt", "antip_tof_jet_eta_pt", HistType::kTH2F, {{60, 0.0, 3.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("antip_tof_ue_eta_pt", "antip_tof_ue_eta_pt", HistType::kTH2F, {{60, 0.0, 3.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});

    // DCA Distributions
    registryData.add("antiproton_dca_jet", "antiproton_dca_jet", HistType::kTH3F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antiproton_dca_ue", "antiproton_dca_ue", HistType::kTH3F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}, {2, 0, 100, "#it{N}_{ch}"}});

    // Antideuterons
    registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH3F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH3F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH3F, {{10, 0.0, 1.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH3F, {{45, 0.5, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}, {2, 0, 100, "#it{N}_{ch}"}});

    // Antihelium-3
    registryData.add("antihelium3_jet_tpc", "antihelium3_jet_tpc", HistType::kTH3F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});
    registryData.add("antihelium3_ue_tpc", "antihelium3_ue_tpc", HistType::kTH3F, {{40, 1.0, 7.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10.0, 10.0, "n#sigma_{TPC}"}, {2, 0, 100, "#it{N}_{ch}"}});

    // Input Antiproton Distribution
    registryMC.add("antiproton_input", "antiproton_input", HistType::kTH1F, {{1000, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_weighted_jet", "antiproton_weighted_jet", HistType::kTH1F, {{1000, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_weighted_ue", "antiproton_weighted_ue", HistType::kTH1F, {{1000, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});

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

    // Histograms for reweighting
    registryMC.add("antiproton_eta_pt_pythia", "antiproton_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
    registryMC.add("antiproton_eta_pt_jet", "antiproton_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
    registryMC.add("antiproton_eta_pt_ue", "antiproton_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
    registryMC.add("antideuteron_eta_pt_injected", "antideuteron_eta_pt_injected", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
    registryMC.add("antideuteron_eta_pt_jet", "antideuteron_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
    registryMC.add("antideuteron_eta_pt_ue", "antideuteron_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {80, -0.8, 0.8, "#eta"}});
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

    // Rapidity Cut
    double mass(0);
    if (particle_of_interest == nucleus::proton)
      mass = 0.93827208816; // Proton
    if (particle_of_interest == nucleus::deuteron)
      mass = 1.87561294257; // Deuteron
    if (particle_of_interest == nucleus::helium)
      mass = 2.80839160743; // Helium-3

    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(track.px(), track.py(), track.pz(), mass);
    if (lorentzVect.Rapidity() < min_y || lorentzVect.Rapidity() > max_y)
      return false;

    return true;
  }

  // Single-Track Selection for Particles inside Jets
  template <typename T2>
  bool passedTrackSelectionJets(const T2& track)
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
    if (TMath::Abs(track.dcaXY()) > 0.5)
      return false;
    if (TMath::Abs(track.dcaZ()) > 0.5)
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
    if (particle_of_interest == nucleus::proton) {
      if (pt < 1.0 && TMath::Abs(nsigmaTPCPr) < 8.0)
        return true;
      if (pt >= 1.0 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 10.0)
        return true;
      return false;
    }

    // Deuteron ID
    if (particle_of_interest == nucleus::deuteron) {
      if (pt < 1.0 && TMath::Abs(nsigmaTPCDe) < 8.0)
        return true;
      if (pt >= 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF() && TMath::Abs(nsigmaTOFDe) < 10.0)
        return true;
      return false;
    }

    // Helium-3 ID
    if (particle_of_interest == nucleus::helium) {
      if ((0.5 * pt) >= 1.0 && TMath::Abs(nsigmaTPCHe) < 8.0)
        return true;
      return false;
    }

    return false;
  }

  template <typename T4>
  bool isHighPurityAntiproton(const T4& track)
  {
    // Variables
    float nsigmaTPCPr = track.tpcNSigmaPr();
    float nsigmaTOFPr = track.tofNSigmaPr();
    float pt = track.pt();

    if (pt < 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0)
      return true;
    if (pt >= 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0 && track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 2.0)
      return true;
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

  float Weight(float pt, int event_region, int nucleus_of_interest)
  {

    if (!applyReweighting)
      return 1;

    auto par_proton_ref = static_cast<std::vector<float>>(param_proton_ref);
    auto par_proton_jet = static_cast<std::vector<float>>(param_proton_jet);
    auto par_deuteron_jet = static_cast<std::vector<float>>(param_deuteron_jet);
    auto par_helium3_jet = static_cast<std::vector<float>>(param_helium3_jet);
    auto par_proton_ue = static_cast<std::vector<float>>(param_proton_ue);
    auto par_deuteron_ue = static_cast<std::vector<float>>(param_deuteron_ue);
    auto par_helium3_ue = static_cast<std::vector<float>>(param_helium3_ue);

    float dNdpt_proton_ref = GetTsallis(par_proton_ref[0], par_proton_ref[1], par_proton_ref[2], par_proton_ref[3], pt);
    float dNdpt_proton_jet = GetTsallis(par_proton_jet[0], par_proton_jet[1], par_proton_jet[2], par_proton_jet[3], pt);
    float dNdpt_proton_ue = GetTsallis(par_proton_ue[0], par_proton_ue[1], par_proton_ue[2], par_proton_ue[3], pt);
    float dNdpt_deuteron_jet = GetTsallis(par_deuteron_jet[0], par_deuteron_jet[1], par_deuteron_jet[2], par_deuteron_jet[3], pt);
    float dNdpt_deuteron_ue = GetTsallis(par_deuteron_ue[0], par_deuteron_ue[1], par_deuteron_ue[2], par_deuteron_ue[3], pt);
    float dNdpt_helium3_jet = GetTsallis(par_helium3_jet[0], par_helium3_jet[1], par_helium3_jet[2], par_helium3_jet[3], pt);
    float dNdpt_helium3_ue = GetTsallis(par_helium3_ue[0], par_helium3_ue[1], par_helium3_ue[2], par_helium3_ue[3], pt);

    if (nucleus_of_interest == nucleus::proton && event_region == region::jet)
      return dNdpt_proton_jet / dNdpt_proton_ref;
    if (nucleus_of_interest == nucleus::proton && event_region == region::underlying_event)
      return dNdpt_proton_ue / dNdpt_proton_ref;
    if (nucleus_of_interest == nucleus::deuteron && event_region == region::jet)
      return dNdpt_deuteron_jet;
    if (nucleus_of_interest == nucleus::deuteron && event_region == region::underlying_event)
      return dNdpt_deuteron_ue;
    if (nucleus_of_interest == nucleus::helium && event_region == region::jet)
      return dNdpt_helium3_jet;
    if (nucleus_of_interest == nucleus::helium && event_region == region::underlying_event)
      return dNdpt_helium3_ue;

    return 1;
  }

  float GetTsallis(float p0, float p1, float p2, float p3, float pt)
  {

    float dNdpt(1);
    float part1 = p0 * pt * (p1 - 1.0) * (p1 - 2.0);
    float part2 = part1 * TMath::Power(1.0 + (TMath::Sqrt(p3 * p3 + pt * pt) - p3) / (p1 * p2), -p1);
    float part3 = part2 / TMath::TwoPi() * (p1 * p2 * (p1 * p2 + p3 * (p1 - 2.0)));
    dNdpt = part3;

    return dNdpt;
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

  // Lorentz Transformation
  TLorentzVector LorentzTransform(TLorentzVector R, TVector3 beta_vect)
  {

    // Inizialization
    TLorentzVector R_prime(0.0, 0.0, 0.0, 0.0);

    // Beta Components
    Double_t Bx = beta_vect.X();
    Double_t By = beta_vect.Y();
    Double_t Bz = beta_vect.Z();

    // Beta & Gamma
    Double_t beta = TMath::Sqrt(Bx * Bx + By * By + Bz * Bz);
    if (beta >= 1.0) {
      return R_prime;
    }
    Double_t gamma = 1.0 / TMath::Sqrt(1.0 - (beta * beta));

    // Coordinates in the Lab System
    Double_t t = R.T();
    Double_t x = R.X();
    Double_t y = R.Y();
    Double_t z = R.Z();

    // Coordinates in the Deuteron Frame
    Double_t t_prime = gamma * t - gamma * Bx * x - gamma * By * y - gamma * Bz * z;
    Double_t x_prime = -gamma * Bx * t + (1.0 + (gamma - 1.0) * Bx * Bx / (beta * beta)) * x + (gamma - 1.0) * (Bx * By / (beta * beta)) * y + (gamma - 1.0) * (Bx * Bz / (beta * beta)) * z;
    Double_t y_prime = -gamma * By * t + (gamma - 1.0) * (Bx * By / (beta * beta)) * x + (1.0 + (gamma - 1.0) * By * By / (beta * beta)) * y + (gamma - 1.0) * (By * Bz / (beta * beta)) * z;
    Double_t z_prime = -gamma * Bz * t + (gamma - 1.0) * (Bx * Bz / (beta * beta)) * x + (gamma - 1.0) * (By * Bz / (beta * beta)) * y + (1.0 + (gamma - 1.0) * Bz * Bz / (beta * beta)) * z;

    // Set Coordinates
    R_prime.SetXYZT(x_prime, y_prime, z_prime, t_prime);

    return R_prime;
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, FullTracks const& tracks)
  {

    // Event Counter: before event selection
    registryQC.fill(HIST("number_of_events_data"), 0.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter: after event selection sel8
    registryQC.fill(HIST("number_of_events_data"), 1.5);

    // Cut on Zvertex
    if (abs(collision.posZ()) > 10.0)
      return;

    // Event Counter: after |z|<10 cm cut
    registryQC.fill(HIST("number_of_events_data"), 2.5);

    // Reduced Event
    std::vector<int> particle_ID;
    int leading_ID = 0;
    bool containsParticleOfInterest(false);
    float pt_max(0);
    int i = -1;

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      // Global Track
      i++;

      // Track Selection for Jet
      if (!passedTrackSelectionJets(track))
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

    // Event Counter: Skip Events with no trigger Particle (pmax=0)
    if (pt_max == 0)
      return;
    registryQC.fill(HIST("number_of_events_data"), 3.5);

    // Histogram with pt_leading
    registryQC.fill(HIST("pt_leading"), pt_max);

    // Event Counter: Skip Events with pt<pt_leading_min
    if (pt_max < min_pt_leading)
      return;
    registryQC.fill(HIST("number_of_events_data"), 4.5);

    // Number of Stored Particles
    int nParticles = static_cast<int>(particle_ID.size());

    // Event Counter: Skip Events with 0 Particles
    if (nParticles < 1)
      return;
    registryQC.fill(HIST("number_of_events_data"), 5.5);

    // Momentum of the Leading Particle
    auto const& leading_track = tracks.iteratorAt(leading_ID);
    TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());
    TVector3 p_highest_pt_track(leading_track.px(), leading_track.py(), leading_track.pz());

    // Event Counter: Skip Events with no Particle of Interest
    if (!containsParticleOfInterest)
      return;
    registryQC.fill(HIST("number_of_events_data"), 6.5);

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
        float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
        float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
        float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

        // Distances
        float distance_jet = min * Delta2 / (Rparameter_jet * Rparameter_jet);
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

    // QA Plots
    registryQC.fill(HIST("eta_leading"), p_leading.Eta());
    registryQC.fill(HIST("phi_leading"), TVector2::Phi_0_2pi(p_leading.Phi()));
    registryQC.fill(HIST("angle_jet_leading_track"), p_leading.Angle(p_highest_pt_track));

    // Rmax and Area Cut
    float Rmax(0);
    int nParticlesJetAndUE(0);
    for (int i = 0; i < nParticlesJetUE; i++) {

      const auto& jet_track = tracks.iteratorAt(jet_particle_ID[i]);
      TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());

      float deltaEta = p_i.Eta() - p_leading.Eta();
      float deltaPhi = GetDeltaPhi(p_i.Phi(), p_leading.Phi());
      float R = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
      if (R < Rmax_jet_ue)
        nParticlesJetAndUE++;
      if (R > Rmax)
        Rmax = R;
    }

    // Rmax Distribution
    registryQC.fill(HIST("r_max_jet"), Rmax);

    // Event Counter: Skip Events with jet not fully inside acceptance
    float eta_jet_axis = p_leading.Eta();
    if ((TMath::Abs(eta_jet_axis) + Rmax_jet_ue) > max_eta)
      return;
    registryQC.fill(HIST("number_of_events_data"), 7.5);

    // Fill Jet Multiplicity
    registryQC.fill(HIST("jet_plus_ue_multiplicity"), nParticlesJetAndUE);

    // Perpendicular Cones for UE Estimate
    TVector3 ue_axis1(0.0, 0.0, 0.0);
    TVector3 ue_axis2(0.0, 0.0, 0.0);
    get_perpendicular_axis(p_leading, ue_axis1, +1.0);
    get_perpendicular_axis(p_leading, ue_axis2, -1.0);

    // Protection against delta<0
    if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
      return;
    if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
      return;

    registryQC.fill(HIST("number_of_events_data"), 8.5);

    // Store UE
    std::vector<int> ue_particle_ID;

    for (int i = 0; i < nParticles; i++) {

      // Skip Leading Particle & Elements already associated to the Jet
      if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
        continue;

      // Get UE Track
      const auto& ue_track = tracks.iteratorAt(particle_ID[i]);

      // Variables
      float deltaEta1 = ue_track.eta() - ue_axis1.Eta();
      float deltaPhi1 = GetDeltaPhi(ue_track.phi(), ue_axis1.Phi());
      float dr1 = TMath::Sqrt(deltaEta1 * deltaEta1 + deltaPhi1 * deltaPhi1);
      float deltaEta2 = ue_track.eta() - ue_axis2.Eta();
      float deltaPhi2 = GetDeltaPhi(ue_track.phi(), ue_axis2.Phi());
      float dr2 = TMath::Sqrt(deltaEta2 * deltaEta2 + deltaPhi2 * deltaPhi2);

      // Store Particles in the UE
      if (dr1 < Rmax_jet_ue) {
        registryQC.fill(HIST("eta_phi_ue"), deltaEta1, deltaPhi1);
        registryQC.fill(HIST("r_ue"), dr1);
        ue_particle_ID.push_back(particle_ID[i]);
      }

      if (dr2 < Rmax_jet_ue) {
        registryQC.fill(HIST("eta_phi_ue"), deltaEta2, deltaPhi2);
        registryQC.fill(HIST("r_ue"), dr2);
        ue_particle_ID.push_back(particle_ID[i]);
      }
    }

    // UE Multiplicity
    int nParticlesUE = static_cast<int>(ue_particle_ID.size());
    registryQC.fill(HIST("ue_multiplicity"), static_cast<float>(nParticlesUE) / 2.0);

    // Jet Multiplicity
    float jet_Nch = static_cast<float>(nParticlesJetAndUE) - static_cast<float>(nParticlesUE) / 2.0;
    registryQC.fill(HIST("jet_multiplicity"), jet_Nch);

    // Loop over particles inside Jet
    for (int i = 0; i < nParticlesJetUE; i++) {

      const auto& jet_track = tracks.iteratorAt(jet_particle_ID[i]);
      TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());

      float deltaEta = p_i.Eta() - p_leading.Eta();
      float deltaPhi = GetDeltaPhi(p_i.Phi(), p_leading.Phi());
      float R = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
      if (R > Rmax_jet_ue)
        continue;

      if (deltaEta != 0 && deltaPhi != 0) {
        registryQC.fill(HIST("eta_phi_jet"), deltaEta, deltaPhi);
        registryQC.fill(HIST("r_jet"), TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi));
      }

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
      float eta = jet_track.eta();

      // DCA
      if (isHighPurityAntiproton(jet_track) && TMath::Abs(jet_track.dcaZ()) < max_dcaz) {
        registryData.fill(HIST("antiproton_dca_jet"), pt, jet_track.dcaXY(), jet_Nch);
      }

      // DCA Cuts
      if (TMath::Abs(jet_track.dcaXY()) > max_dcaxy)
        continue;
      if (TMath::Abs(jet_track.dcaZ()) > max_dcaz)
        continue;

      // Antiproton
      if (particle_of_interest == nucleus::proton) {
        if (pt < 1.0)
          registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr, jet_Nch);
        if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && jet_track.hasTOF())
          registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr, jet_Nch);
        if (pt < 1.0 && TMath::Abs(nsigmaTPCPr) < 2.0)
          registryData.fill(HIST("antip_tpc_jet_eta_pt"), pt, eta);
        if (pt >= 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0 && jet_track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 2.0)
          registryData.fill(HIST("antip_tof_jet_eta_pt"), pt, eta);
      }

      // Antideuteron
      if (particle_of_interest == nucleus::deuteron) {
        if (pt < 1.0)
          registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe, jet_Nch);
        if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && jet_track.hasTOF())
          registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe, jet_Nch);
      }

      // Antihelium3
      if (particle_of_interest == nucleus::helium) {
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
      float eta = ue_track.eta();

      // DCA
      if (isHighPurityAntiproton(ue_track) && TMath::Abs(ue_track.dcaZ()) < max_dcaz) {
        registryData.fill(HIST("antiproton_dca_ue"), pt, ue_track.dcaXY(), jet_Nch);
      }

      // DCA Cuts
      if (TMath::Abs(ue_track.dcaXY()) > max_dcaxy)
        continue;
      if (TMath::Abs(ue_track.dcaZ()) > max_dcaz)
        continue;

      // Antiproton
      if (particle_of_interest == nucleus::proton) {
        if (pt < 1.0)
          registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr, jet_Nch);
        if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && ue_track.hasTOF())
          registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr, jet_Nch);
        if (pt < 1.0 && TMath::Abs(nsigmaTPCPr) < 2.0)
          registryData.fill(HIST("antip_tpc_ue_eta_pt"), pt, eta);
        if (pt >= 0.5 && TMath::Abs(nsigmaTPCPr) < 2.0 && ue_track.hasTOF() && TMath::Abs(nsigmaTOFPr) < 2.0)
          registryData.fill(HIST("antip_tof_ue_eta_pt"), pt, eta);
      }

      // Antideuteron
      if (particle_of_interest == nucleus::deuteron) {
        if (pt < 1.0)
          registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe, jet_Nch);
        if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && ue_track.hasTOF())
          registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe, jet_Nch);
      }

      // Antihelium3
      if (particle_of_interest == nucleus::helium) {
        registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe, jet_Nch);
      }
    }
  }

  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<MCTracks> perCollision = o2::aod::track::collisionId;

  void processGen(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mcCollisions) {

      // Event Counter (before event sel)
      registryQC.fill(HIST("number_of_events_mc"), 0.5);

      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // Generated Particles
      for (auto& particle : mcParticles_per_coll) {

        if (!particle.isPhysicalPrimary())
          continue;
        if ((particle.pdgCode() != -2212) && (particle.pdgCode() != -1000010020) && (particle.pdgCode() != -1000020030))
          continue;
        if (particle.y() < min_y || particle.y() > max_y)
          continue;

        // Reweighting
        float wpr_jet = Weight(particle.pt(), region::jet, nucleus::proton);
        float wpr_ue = Weight(particle.pt(), region::underlying_event, nucleus::proton);
        float wde_jet = Weight(particle.pt(), region::jet, nucleus::deuteron);
        float wde_ue = Weight(particle.pt(), region::underlying_event, nucleus::deuteron);
        float whe_jet = Weight(particle.pt(), region::jet, nucleus::helium);
        float whe_ue = Weight(particle.pt(), region::underlying_event, nucleus::helium);

        // Fill Histograms
        if (particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_input"), particle.pt());
          registryMC.fill(HIST("antiproton_weighted_jet"), particle.pt(), wpr_jet);
          registryMC.fill(HIST("antiproton_weighted_ue"), particle.pt(), wpr_ue);
          registryMC.fill(HIST("antiproton_jet_gen"), particle.pt(), wpr_jet);
          registryMC.fill(HIST("antiproton_ue_gen"), particle.pt(), wpr_ue);
          registryMC.fill(HIST("antiproton_eta_pt_pythia"), particle.pt(), particle.eta());
        }
        if (particle.pdgCode() == -1000010020) {
          registryMC.fill(HIST("antideuteron_jet_gen"), particle.pt(), wde_jet);
          registryMC.fill(HIST("antideuteron_ue_gen"), particle.pt(), wde_ue);
          registryMC.fill(HIST("antideuteron_eta_pt_injected"), particle.pt(), particle.eta());
        }
        if (particle.pdgCode() == -1000020030) {
          registryMC.fill(HIST("antihelium3_jet_gen"), particle.pt(), whe_jet);
          registryMC.fill(HIST("antihelium3_ue_gen"), particle.pt(), whe_ue);
        }
      }
    }
  }

  void processRec(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {

    for (const auto& collision : collisions) {

      // Event Counter (before event sel)
      registryQC.fill(HIST("number_of_events_mc"), 1.5);

      // Event Selection
      if (!collision.sel8())
        continue;

      if (abs(collision.posZ()) > 10)
        continue;

      // Event Counter (after event sel)
      registryQC.fill(HIST("number_of_events_mc"), 2.5);

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
        if (require_primVtx_contributor && !(track.isPVContributor()))
          continue;

        // Reweighting
        float wpr_jet = Weight(particle.pt(), region::jet, nucleus::proton);
        float wpr_ue = Weight(particle.pt(), region::underlying_event, nucleus::proton);
        float wde_jet = Weight(particle.pt(), region::jet, nucleus::deuteron);
        float wde_ue = Weight(particle.pt(), region::underlying_event, nucleus::deuteron);
        float whe_jet = Weight(particle.pt(), region::jet, nucleus::helium);
        float whe_ue = Weight(particle.pt(), region::underlying_event, nucleus::helium);

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
            registryMC.fill(HIST("antiproton_jet_rec_tpc"), pt, wpr_jet);
            registryMC.fill(HIST("antiproton_ue_rec_tpc"), pt, wpr_ue);
          }
          if (pt >= 0.5 && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF() && nsigmaTOFPr > min_nsigmaTOF && nsigmaTOFPr < max_nsigmaTOF) {
            registryMC.fill(HIST("antiproton_jet_rec_tof"), pt, wpr_jet);
            registryMC.fill(HIST("antiproton_ue_rec_tof"), pt, wpr_ue);
          }
        }

        // Antideuteron
        if (particle.pdgCode() == -1000010020) {
          if (pt < 1.0 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC) {
            registryMC.fill(HIST("antideuteron_jet_rec_tpc"), pt, wde_jet);
            registryMC.fill(HIST("antideuteron_ue_rec_tpc"), pt, wde_ue);
          }
          if (pt >= 0.5 && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF() && nsigmaTOFDe > min_nsigmaTOF && nsigmaTOFDe < max_nsigmaTOF) {
            registryMC.fill(HIST("antideuteron_jet_rec_tof"), pt, wde_jet);
            registryMC.fill(HIST("antideuteron_ue_rec_tof"), pt, wde_ue);
          }
        }

        // Antihelium-3
        if (particle.pdgCode() == -1000020030) {
          if (nsigmaTPCHe > min_nsigmaTPC && nsigmaTPCHe < max_nsigmaTPC) {
            registryMC.fill(HIST("antihelium3_jet_rec_tpc"), 2.0 * pt, whe_jet);
            registryMC.fill(HIST("antihelium3_ue_rec_tpc"), 2.0 * pt, whe_ue);
          }
        }
      }
    }
  }

  void processAntiNucleiJet(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {

    // Loop over MC Collisions
    for (const auto& mccollision : mcCollisions) {

      registryQC.fill(HIST("event_counter_mc"), 0.5);

      // Selection on z_{vertex}
      if (abs(mccollision.posZ()) > 10)
        continue;
      registryQC.fill(HIST("event_counter_mc"), 1.5);

      // MC Particles per Collision
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // Reduced Event
      std::vector<int> particle_ID;
      int leading_ID = 0;
      float pt_max(0);

      // Generated Particles
      for (auto& particle : mcParticles_per_coll) {

        // Global Index
        int i = particle.globalIndex();

        // Select Primary Particles
        float deltaX = particle.vx() - mccollision.posX();
        float deltaY = particle.vy() - mccollision.posY();
        float deltaZ = particle.vz() - mccollision.posZ();
        float deltaR = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        if (deltaR > 0.1)
          continue;

        // Pseudorapidity Selection
        if (abs(particle.eta()) > 0.8)
          continue;

        // PDG Selection
        int pdg = abs(particle.pdgCode());
        if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212) && (pdg != 1000010020))
          continue;

        // Find pt Leading
        if (particle.pt() > pt_max) {
          leading_ID = i;
          pt_max = particle.pt();
        }

        // Store Array Element
        particle_ID.push_back(i);
      }

      // Skip Events with pt<pt_leading_min
      if (pt_max < min_pt_leading)
        continue;
      registryQC.fill(HIST("event_counter_mc"), 2.5);

      // Number of Stored Particles
      int nParticles = static_cast<int>(particle_ID.size());

      // Momentum of the Leading Particle
      auto const& leading_track = mcParticles_per_coll.iteratorAt(leading_ID);
      TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());
      registryQC.fill(HIST("deltaPt"), leading_track.pt() - pt_max);

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
          auto stored_track = mcParticles_per_coll.iteratorAt(particle_ID[i]);
          TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

          // Variables
          float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
          float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
          float deltaEta = p_particle.Eta() - p_leading.Eta();
          float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
          float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
          float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

          // Distances
          float distance_jet = min * Delta2 / (Rparameter_jet * Rparameter_jet);
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
          auto jet_track = mcParticles_per_coll.iteratorAt(label_jet_particle);
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

      // Event Counter: Skip Events with jet not fully inside acceptance
      float eta_jet_axis = p_leading.Eta();
      if ((TMath::Abs(eta_jet_axis) + Rmax_jet_ue) > max_eta)
        continue;
      registryQC.fill(HIST("event_counter_mc"), 3.5);

      // Perpendicular Cones for UE Estimate
      TVector3 ue_axis1(0.0, 0.0, 0.0);
      TVector3 ue_axis2(0.0, 0.0, 0.0);
      get_perpendicular_axis(p_leading, ue_axis1, +1.0);
      get_perpendicular_axis(p_leading, ue_axis2, -1.0);

      // Protection against delta<0
      if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
        continue;
      if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
        continue;
      registryQC.fill(HIST("event_counter_mc"), 4.5);

      // Store UE
      std::vector<int> ue_particle_ID;

      for (int i = 0; i < nParticles; i++) {

        // Skip Leading Particle & Elements already associated to the Jet
        if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
          continue;

        // Get UE Track
        const auto& ue_track = mcParticles_per_coll.iteratorAt(particle_ID[i]);

        // Variables
        float deltaEta1 = ue_track.eta() - ue_axis1.Eta();
        float deltaPhi1 = GetDeltaPhi(ue_track.phi(), ue_axis1.Phi());
        float dr1 = TMath::Sqrt(deltaEta1 * deltaEta1 + deltaPhi1 * deltaPhi1);
        float deltaEta2 = ue_track.eta() - ue_axis2.Eta();
        float deltaPhi2 = GetDeltaPhi(ue_track.phi(), ue_axis2.Phi());
        float dr2 = TMath::Sqrt(deltaEta2 * deltaEta2 + deltaPhi2 * deltaPhi2);

        // Store Particles in the UE 1
        if (dr1 < Rmax_jet_ue) {
          ue_particle_ID.push_back(particle_ID[i]);
        }

        // Store Particles in the UE 2
        if (dr2 < Rmax_jet_ue) {
          ue_particle_ID.push_back(particle_ID[i]);
        }
      }

      // UE Multiplicity
      int nParticlesUE = static_cast<int>(ue_particle_ID.size());

      // Loop over particles inside Jet
      for (int i = 0; i < nParticlesJetUE; i++) {

        const auto& jet_track = mcParticles_per_coll.iteratorAt(jet_particle_ID[i]);

        float deltaEta = jet_track.eta() - p_leading.Eta();
        float deltaPhi = GetDeltaPhi(jet_track.phi(), p_leading.Phi());
        float R = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
        if (R > Rmax_jet_ue)
          continue;

        // Rapidity Cut
        if (jet_track.y() < min_y || jet_track.y() > max_y)
          continue;

        if (jet_track.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_eta_pt_jet"), jet_track.pt(),
                          jet_track.eta());
        }

        if (jet_track.pdgCode() == -1000010020) {
          registryMC.fill(HIST("antideuteron_eta_pt_jet"), jet_track.pt(),
                          jet_track.eta());
        }
      }

      // Loop over particles inside UE
      for (int i = 0; i < nParticlesUE; i++) {

        const auto& ue_track = mcParticles_per_coll.iteratorAt(ue_particle_ID[i]);

        // Rapidity Cut
        if (ue_track.y() < min_y || ue_track.y() > max_y)
          continue;

        if (ue_track.pdgCode() != -2212) {
          registryMC.fill(HIST("antiproton_eta_pt_ue"), ue_track.pt(),
                          ue_track.eta());
        }

        if (ue_track.pdgCode() != -1000010020) {
          registryMC.fill(HIST("antideuteron_eta_pt_ue"), ue_track.pt(),
                          ue_track.eta());
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

      if (abs(collision.posZ()) > 10)
        continue;

      auto tracks_per_coll = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // Reduced Event
      std::vector<int> particle_ID;
      int leading_ID = 0;
      float pt_max(0);
      int i = -1;

      // Loop over Reconstructed Tracks
      for (auto track : tracks_per_coll) {

        // Global Track Index
        i++;

        // Track Selection for Jet
        if (!passedTrackSelectionJets(track))
          continue;

        // Find pt Leading
        if (track.pt() > pt_max) {
          leading_ID = i;
          pt_max = track.pt();
        }

        // Store Array Element
        particle_ID.push_back(i);
      }

      if (pt_max < min_pt_leading)
        continue;

      // Number of Stored Particles
      int nParticles = static_cast<int>(particle_ID.size());

      // Event Counter: Skip Events with 0 Particles
      if (nParticles < 1)
        continue;

      // Momentum of the Leading Particle
      auto const& leading_track = tracks_per_coll.iteratorAt(leading_ID);
      TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());

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
          auto stored_track = tracks_per_coll.iteratorAt(particle_ID[i]);
          TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

          // Variables
          float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
          float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
          float deltaEta = p_particle.Eta() - p_leading.Eta();
          float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
          float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
          float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

          // Distances
          float distance_jet = min * Delta2 / (Rparameter_jet * Rparameter_jet);
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

          // Update Momentum of Leading Particle
          auto jet_track = tracks_per_coll.iteratorAt(label_jet_particle);
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

      float eta_jet_axis = p_leading.Eta();
      if ((TMath::Abs(eta_jet_axis) + Rmax_jet_ue) > max_eta)
        continue;
      TVector3 jet_axis(p_leading.X(), p_leading.Y(), p_leading.Z());

      // Perpendicular Cones for UE Estimate
      TVector3 ue_axis1(0.0, 0.0, 0.0);
      TVector3 ue_axis2(0.0, 0.0, 0.0);
      get_perpendicular_axis(jet_axis, ue_axis1, +1.0);
      get_perpendicular_axis(jet_axis, ue_axis2, -1.0);

      // Protection against delta<0
      if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
        continue;
      if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
        continue;

      for (auto track : tracks_per_coll) {

        if (!passedTrackSelection(track))
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
        float deltaEta_jet = particle_dir.Eta() - jet_axis.Eta();
        float deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), jet_axis.Phi());
        float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        float deltaEta_ue1 = particle_dir.Eta() - ue_axis1.Eta();
        float deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue_axis1.Phi());
        float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        float deltaEta_ue2 = particle_dir.Eta() - ue_axis2.Eta();
        float deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue_axis2.Phi());
        float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

        if (deltaR_jet < Rmax_jet_ue) {
          registryMC.fill(HIST("antiproton_all_jet"), track.pt());
          if (particle.isPhysicalPrimary())
            registryMC.fill(HIST("antiproton_prim_jet"), track.pt());
        }
        if (deltaR_ue1 < Rmax_jet_ue || deltaR_ue2 < Rmax_jet_ue) {
          registryMC.fill(HIST("antiproton_all_ue"), track.pt());
          if (particle.isPhysicalPrimary())
            registryMC.fill(HIST("antiproton_prim_ue"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_jets, processSecAntiprotons, "Process sec antip", false);

  /*
  void processAntidJet(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {

    // Particle Masses [GeV/c^2]
    Double_t mp = 0.93827208816;
    Double_t mn = 0.93956542052;
    Double_t md = 1.87561294257;

    for (const auto& mccollision : mcCollisions) {

      // Selection on z_{vertex}
      if (abs(mccollision.posZ()) > 10)
        continue;

      // MC Particles per Collision
      auto mcParticles_per_coll =
        mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // Reduced Event
      std::vector<int> particle_ID;
      int leading_ID = 0;
      float pt_max(0);

      // Specific Flags
      bool containsProton(false);
      bool containsNeutron(false);

      // Generated Particles
      for (auto& particle : mcParticles_per_coll) {

        // Global Index
        int i = particle.globalIndex();

        // Select Primary Particles
        float deltaX = particle.vx() - mccollision.posX();
        float deltaY = particle.vy() - mccollision.posY();
        float deltaZ = particle.vz() - mccollision.posZ();
        float deltaR = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        if (deltaR > 0.1)
          continue;

        // Pseudorapidity Selection
        if (abs(particle.eta()) > 0.8)
          continue;

        // PDG Selection
        int pdg = abs(particle.pdgCode());
        if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212) && (pdg != 2112))
          continue;

        // Nucleon Selection
        if (particle.pdgCode() == -2212)
          containsProton = true;
        if (particle.pdgCode() == -2112)
          containsNeutron = true;

        // Find pt Leading
        if (particle.pt() > pt_max && (pdg != 2112)) {
          leading_ID = i;
          pt_max = particle.pt();
        }

        // Store Array Element
        particle_ID.push_back(i);
      }

      // Skip Events with pt<pt_leading_min
      if (pt_max < min_pt_leading)
          continue;

      // Number of Stored Particles
      int nParticles = static_cast<int>(particle_ID.size());

      if (nParticles < 2)
        continue;
      if (!containsProton)
        continue;
      if (!containsNeutron)
        continue;

      // Momentum of the Leading Particle
      auto const& leading_track = mcParticles_per_coll.iteratorAt(leading_ID);
      TVector3 p_leading(leading_track.px(), leading_track.py(),
                         leading_track.pz());

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
          auto stored_track = mcParticles_per_coll.iteratorAt(particle_ID[i]);
          TVector3 p_particle(stored_track.px(), stored_track.py(),
                              stored_track.pz());

          // Variables
          float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
          float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
          float deltaEta = p_particle.Eta() - p_leading.Eta();
          float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
          float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
          float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

          // Distances
          float distance_jet = min * Delta2 / (Rparameter_jet * Rparameter_jet);
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
          auto jet_track = mcParticles_per_coll.iteratorAt(label_jet_particle);
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

      // Event Counter: Skip Events with jet not fully inside acceptance
      float eta_jet_axis = p_leading.Eta();
      if ((TMath::Abs(eta_jet_axis) + Rmax_jet_ue) > max_eta)
          continue;

      // Perpendicular Cones for UE Estimate
      TVector3 ue_axis1(0.0, 0.0, 0.0);
      TVector3 ue_axis2(0.0, 0.0, 0.0);
      get_perpendicular_axis(p_leading, ue_axis1, +1.0);
      get_perpendicular_axis(p_leading, ue_axis2, -1.0);

      // Protection against delta<0
      if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
          continue;
      if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
          continue;

      // Store UE
      std::vector<int> ue_particle_ID;

      for (int i = 0; i < nParticles; i++) {

        // Skip Leading Particle & Elements already associated to the Jet
        if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
          continue;

        // Get UE Track
        const auto& ue_track = mcParticles_per_coll.iteratorAt(particle_ID[i]);

        // Variables
        float deltaEta1 = ue_track.eta() - ue_axis1.Eta();
        float deltaPhi1 = GetDeltaPhi(ue_track.phi(), ue_axis1.Phi());
        float dr1 = TMath::Sqrt(deltaEta1 * deltaEta1 + deltaPhi1 * deltaPhi1);
        float deltaEta2 = ue_track.eta() - ue_axis2.Eta();
        float deltaPhi2 = GetDeltaPhi(ue_track.phi(), ue_axis2.Phi());
        float dr2 = TMath::Sqrt(deltaEta2 * deltaEta2 + deltaPhi2 * deltaPhi2);

        // Store Particles in the UE
        if (dr1 < Rmax_jet_ue) {
          ue_particle_ID.push_back(particle_ID[i]);
        }

        if (dr2 < Rmax_jet_ue) {
          ue_particle_ID.push_back(particle_ID[i]);
        }
      }

      // UE Multiplicity
      int nParticlesUE = static_cast<int>(ue_particle_ID.size());

      vector<int> neutron_status_jet;
      for (int i = 0; i < nParticlesJetUE; i++) {

        neutron_status_jet.push_back(0);
      }

      // Loop over particles inside Jet
      for (int i = 0; i < nParticlesJetUE; i++) {

        // Proton Selection
        const auto& proton_jet = mcParticles_per_coll.iteratorAt(jet_particle_ID[i]);
        if (proton_jet.pdgCode() != -2212)
          continue;

        // Proton 4-Momentum
        TLorentzVector p_proton;
        p_proton.SetXYZM(proton_jet.px(), proton_jet.py(), proton_jet.pz(), mp);

        for (int j = 0; j < nParticlesJetUE; j++) {

          // Neutron Selection
          const auto& neutron_jet = mcParticles_per_coll.iteratorAt(jet_particle_ID[j]);
          if (neutron_jet.pdgCode() != -2112)
            continue;

          // Neutron 4-Momentum
          TLorentzVector p_neutron;
          p_neutron.SetXYZM(neutron_jet.px(), neutron_jet.py(), neutron_jet.pz(), mn);

          // Deuteron 4-Momentum
          TLorentzVector p_deuteron;
          p_deuteron.SetXYZM(p_proton.X() + p_neutron.X(), p_proton.Y() + p_neutron.Y(), p_proton.Z() + p_neutron.Z(), md);
          Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
          Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
          Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
          TVector3 beta(beta_x, beta_y, beta_z);

          // Lorentz Transformation
          TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
          TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

          // Variables
          Double_t deltaP = (p_proton_prime - p_neutron_prime).P();

          // Skip already used Neutrons
          if (neutron_status_jet[j] == 1)
            continue;

          // Coalescence Condition
          if (deltaP < coalMom) {

            neutron_status_jet[j] = 1;
            Double_t y = p_deuteron.Rapidity();

            float deltaEta = p_deuteron.Eta() - p_leading.Eta();
            float deltaPhi = GetDeltaPhi(p_deuteron.Phi(), p_leading.Phi());
            float R = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

            if (R < Rmax_jet_ue && TMath::Abs(y) < 0.5) {

              registryMC.fill(HIST("antideuteron_eta_pt_jet"), p_deuteron.Pt(), p_deuteron.Eta());
            }
            break;
          }
        }
      }

      vector<int> neutron_status_ue;
      for (int i = 0; i < nParticlesUE; i++) {

        neutron_status_ue.push_back(0);
      }

      // Loop over particles inside UE
      for (int i = 0; i < nParticlesUE; i++) {

        // Proton Selection
        const auto& proton_ue = mcParticles_per_coll.iteratorAt(ue_particle_ID[i]);
        if (proton_ue.pdgCode() != -2212)
          continue;

        // Proton 4-Momentum
        TLorentzVector p_proton;
        p_proton.SetXYZM(proton_ue.px(), proton_ue.py(), proton_ue.pz(), mp);

        for (int j = 0; j < nParticlesUE; j++) {

          // Neutron Selection
          const auto& neutron_ue = mcParticles_per_coll.iteratorAt(ue_particle_ID[j]);
          if (neutron_ue.pdgCode() != -2112)
            continue;

          // Neutron 4-Momentum
          TLorentzVector p_neutron;
          p_neutron.SetXYZM(neutron_ue.px(), neutron_ue.py(), neutron_ue.pz(), mn);

          // Deuteron 4-Momentum
          TLorentzVector p_deuteron;
          p_deuteron.SetXYZM(p_proton.X() + p_neutron.X(), p_proton.Y() + p_neutron.Y(), p_proton.Z() + p_neutron.Z(), md);
          Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
          Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
          Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
          TVector3 beta(beta_x, beta_y, beta_z);

          // Lorentz Transformation
          TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
          TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

          // Variables
          Double_t deltaP = (p_proton_prime - p_neutron_prime).P();

          // Skip already used Neutrons
          if (neutron_status_ue[j] == 1)
            continue;

          // Coalescence Condition
          if (deltaP < coalMom) {

            neutron_status_ue[j] = 1;
            Double_t y = p_deuteron.Rapidity();

            float deltaEta = p_deuteron.Eta() - p_leading.Eta();
            float deltaPhi = GetDeltaPhi(p_deuteron.Phi(), p_leading.Phi());
            float R = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

            if (R < Rmax_jet_ue && TMath::Abs(y) < 0.5) {

              registryMC.fill(HIST("antideuteron_eta_pt_ue"), p_deuteron.Pt(), p_deuteron.Eta());
            }
            break;
          }
        }
      }
    }
  }
  */

  PROCESS_SWITCH(nuclei_in_jets, processData, "Process data", true);
  PROCESS_SWITCH(nuclei_in_jets, processGen, "process Gen MC", false);
  PROCESS_SWITCH(nuclei_in_jets, processRec, "process Rec MC", false);
  PROCESS_SWITCH(nuclei_in_jets, processAntiNucleiJet, "process antinuclei in jet and UE (MC)", false);
  // PROCESS_SWITCH(nuclei_in_jets, processAntidJet, "process antideuterons in jet and UE (MC)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nuclei_in_jets>(cfgc)};
}
