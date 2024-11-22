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
#include <TList.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "TGrid.h"
#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
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
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

using FullNucleiTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe>;

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
  Configurable<double> min_jet_pt{"min_jet_pt", 10.0, "Minimum pt of the jet"};
  Configurable<double> Rjet{"Rjet", 0.3, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<int> min_nPartInJet{"min_nPartInJet", 2, "Minimum number of particles inside jet"};
  Configurable<int> n_jets_per_event_max{"n_jets_per_event_max", 1000, "Maximum number of jets per event"};
  Configurable<bool> requireNoOverlap{"requireNoOverlap", false, "require no overlap between jets and UE cones"};

  // Track Parameters
  Configurable<double> par0{"par0", 0.004, "par 0"};
  Configurable<double> par1{"par1", 0.013, "par 1"};
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 5, "minimum number of ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<double> max_chi2_TPC{"max_chi2_TPC", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> max_chi2_ITS{"max_chi2_ITS", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> min_pt{"min_pt", 0.3, "minimum pt of the tracks"};
  Configurable<double> min_eta{"min_eta", -0.8, "minimum eta"};
  Configurable<double> max_eta{"max_eta", +0.8, "maximum eta"};
  Configurable<double> max_dcaxy{"max_dcaxy", 0.05, "Maximum DCAxy"};
  Configurable<double> max_dcaz{"max_dcaz", 0.05, "Maximum DCAz"};
  Configurable<double> min_nsigmaTPC{"min_nsigmaTPC", -3.0, "Minimum nsigma TPC"};
  Configurable<double> max_nsigmaTPC{"max_nsigmaTPC", +3.0, "Maximum nsigma TPC"};
  Configurable<double> min_nsigmaTOF{"min_nsigmaTOF", -3.0, "Minimum nsigma TOF"};
  Configurable<double> max_nsigmaTOF{"max_nsigmaTOF", +3.5, "Maximum nsigma TOF"};
  Configurable<double> max_pt_for_nsigmaTPC{"max_pt_for_nsigmaTPC", 2.0, "Maximum pt for TPC analysis"};
  Configurable<double> min_pt_for_nsigmaTOF{"min_pt_for_nsigmaTOF", 0.5, "Minimum pt for TOF analysis"};
  Configurable<bool> require_PV_contributor{"require_PV_contributor", true, "require that the track is a PV contributor"};
  Configurable<bool> setDCAselectionPtDep{"setDCAselectionPtDep", true, "require pt dependent selection"};
  Configurable<bool> applyReweighting{"applyReweighting", true, "apply reweighting"};

  Configurable<std::string> url_to_ccdb{"url_to_ccdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> path_to_file{"path_to_file", "", "path to file with reweighting"};
  Configurable<std::string> histo_name_weight_antip_jet{"histo_name_weight_antip_jet", "", "reweighting histogram: antip in jet"};
  Configurable<std::string> histo_name_weight_antip_ue{"histo_name_weight_antip_ue", "", "reweighting histogram: antip in ue"};

  TH2F* twod_weights_antip_jet;
  TH2F* twod_weights_antip_ue;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  void init(InitContext const&)
  {
    ccdb->setURL(url_to_ccdb.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (applyReweighting) {
      GetReweightingHistograms(ccdb, TString(path_to_file), TString(histo_name_weight_antip_jet), TString(histo_name_weight_antip_ue));
    } else {
      twod_weights_antip_jet = nullptr;
      twod_weights_antip_ue = nullptr;
    }

    // QC Histograms
    registryQC.add("deltaEtadeltaPhi_jet", "deltaEtadeltaPhi_jet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhi_ue", "deltaEtadeltaPhi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQC.add("NchJetPlusUE", "NchJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchJet", "NchJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchUE", "NchUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("sumPtJetPlusUE", "sumPtJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtJet", "sumPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtUE", "sumPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("nJets_found", "nJets_found", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("nJets_selected", "nJets_selected", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("dcaxy_vs_pt", "dcaxy_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    registryQC.add("dcaz_vs_pt", "dcaz_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    registryQC.add("jet_ue_overlaps", "jet_ue_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});

    // Event Counters
    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{10, 0, 10, "counter"}});
    registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});

    // Binning
    double min = 0.0;
    double max = 6.0;
    int nbins = 120;

    // Antiprotons
    registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antiproton_dca_jet", "antiproton_dca_jet", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});
    registryData.add("antiproton_dca_ue", "antiproton_dca_ue", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});

    // Antideuterons
    registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

    // Deuterons
    registryData.add("deuteron_jet_tof", "deuteron_jet_tof", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    registryData.add("deuteron_ue_tof", "deuteron_ue_tof", HistType::kTH2F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

    // Antihelium-3
    registryData.add("antihelium3_jet_tpc", "antihelium3_jet_tpc", HistType::kTH2F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("antihelium3_ue_tpc", "antihelium3_ue_tpc", HistType::kTH2F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});

    // Helium-3
    registryData.add("helium3_jet_tpc", "helium3_jet_tpc", HistType::kTH2F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    registryData.add("helium3_ue_tpc", "helium3_ue_tpc", HistType::kTH2F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});

    // Generated
    registryMC.add("antiproton_jet_gen", "antiproton_jet_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_gen", "antideuteron_jet_gen", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_jet_gen", "antihelium3_jet_gen", HistType::kTH1F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_gen", "antiproton_ue_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_gen", "antideuteron_ue_gen", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_ue_gen", "antihelium3_ue_gen", HistType::kTH1F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}});

    // Reconstructed TPC
    registryMC.add("antiproton_jet_rec_tpc", "antiproton_jet_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_rec_tpc", "antideuteron_jet_rec_tpc", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_jet_rec_tpc", "antihelium3_jet_rec_tpc", HistType::kTH1F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_rec_tpc", "antiproton_ue_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_rec_tpc", "antideuteron_ue_rec_tpc", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antihelium3_ue_rec_tpc", "antihelium3_ue_rec_tpc", HistType::kTH1F, {{nbins, min * 3, max * 3, "#it{p}_{T} (GeV/#it{c})"}});

    // Reconstructed TOF
    registryMC.add("antiproton_jet_rec_tof", "antiproton_jet_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_jet_rec_tof", "antideuteron_jet_rec_tof", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_ue_rec_tof", "antiproton_ue_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antideuteron_ue_rec_tof", "antideuteron_ue_rec_tof", HistType::kTH1F, {{nbins, min * 2, max * 2, "#it{p}_{T} (GeV/#it{c})"}});

    // DCA Templates
    registryMC.add("antiproton_dca_prim", "antiproton_dca_prim", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});
    registryMC.add("antiproton_dca_sec", "antiproton_dca_sec", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -0.5, 0.5, "DCA_{xy} (cm)"}});

    // Fraction of Primary Antiprotons from MC
    registryMC.add("antiproton_prim", "antiproton_prim", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all", "antiproton_all", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_prim_jet", "antiproton_prim_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_jet", "antiproton_all_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_prim_ue", "antiproton_prim_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("antiproton_all_ue", "antiproton_all_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

    // Antiproton Reweighting
    registryMC.add("antiproton_eta_pt_pythia", "antiproton_eta_pt_pythia", HistType::kTH2F, {{200, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {20, -1.0, 1.0, "#it{#eta}"}});
    registryMC.add("antiproton_eta_pt_jet", "antiproton_eta_pt_jet", HistType::kTH2F, {{200, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {20, -1.0, 1.0, "#it{#eta}"}});
    registryMC.add("antiproton_eta_pt_ue", "antiproton_eta_pt_ue", HistType::kTH2F, {{200, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {20, -1.0, 1.0, "#it{#eta}"}});
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
    if (track.pt() < 0.15)
      return false;

    // pt-dependent selection
    if (setDCAselectionPtDep) {
      if (TMath::Abs(track.dcaXY()) > (par0 + par1 / track.pt()))
        return false;
      if (TMath::Abs(track.dcaZ()) > (par0 + par1 / track.pt()))
        return false;
    }

    // standard selection
    if (!setDCAselectionPtDep) {
      if (TMath::Abs(track.dcaXY()) > max_dcaxy)
        return false;
      if (TMath::Abs(track.dcaZ()) > max_dcaz)
        return false;
    }

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

    // Protection against delta<0
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

  double calculate_dij(TVector3 t1, TVector3 t2, double R)
  {
    double distance_jet(0);
    double x1 = 1.0 / (t1.Pt() * t1.Pt());
    double x2 = 1.0 / (t2.Pt() * t2.Pt());
    double deltaEta = t1.Eta() - t2.Eta();
    double deltaPhi = GetDeltaPhi(t1.Phi(), t2.Phi());
    double min = Minimum(x1, x2);
    double Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;
    distance_jet = min * Delta2 / (R * R);
    return distance_jet;
  }

  bool overlap(TVector3 v1, TVector3 v2, double R)
  {
    double dx = v1.Eta() - v2.Eta();
    double dy = GetDeltaPhi(v1.Phi(), v2.Phi());
    double d = sqrt(dx * dx + dy * dy);
    if (d < 2.0 * R)
      return true;
    return false;
  }

  void GetReweightingHistograms(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, TString filepath, TString histname_antip_jet, TString histname_antip_ue)
  {
    TList* l = ccdbObj->get<TList>(filepath.Data());
    if (!l) {
      LOGP(error, "Could not open the file {}", Form("%s", filepath.Data()));
      return;
    }
    twod_weights_antip_jet = static_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname_antip_jet.Data())));
    if (!twod_weights_antip_jet) {
      LOGP(error, "Could not open histogram {}", Form("%s_antiproton", histname_antip_jet.Data()));
      return;
    }
    twod_weights_antip_ue = static_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname_antip_ue.Data())));
    if (!twod_weights_antip_ue) {
      LOGP(error, "Could not open histogram {}", Form("%s_antiproton", histname_antip_ue.Data()));
      return;
    }
    LOGP(info, "Opened histogram {}", Form("%s_antiproton", histname_antip_jet.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_antiproton", histname_antip_ue.Data()));
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, FullNucleiTracks const& tracks)
  {
    // Event Counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter: after event selection sel8
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Cut on z-vertex
    if (TMath::Abs(collision.posZ()) > zVtx)
      return;

    // Event Counter: after z-vertex cut
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // List of Tracks
    std::vector<TVector3> trk;

    for (auto track : tracks) {

      if (!passedTrackSelectionForJetReconstruction(track))
        continue;
      registryQC.fill(HIST("dcaxy_vs_pt"), track.pt(), track.dcaXY());
      registryQC.fill(HIST("dcaz_vs_pt"), track.pt(), track.dcaZ());

      TVector3 momentum(track.px(), track.py(), track.pz());
      trk.push_back(momentum);
    }

    // Anti-kt Jet Finder
    int n_particles_removed(0);
    std::vector<TVector3> jet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    do {
      double dij_min(1e+06), diB_min(1e+06);
      int i_min(0), j_min(0), iB_min(0);
      for (int i = 0; i < static_cast<int>(trk.size()); i++) {
        if (trk[i].Mag() == 0)
          continue;
        double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
        if (diB < diB_min) {
          diB_min = diB;
          iB_min = i;
        }
        for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) {
          if (trk[j].Mag() == 0)
            continue;
          double dij = calculate_dij(trk[i], trk[j], Rjet);
          if (dij < dij_min) {
            dij_min = dij;
            i_min = i;
            j_min = j;
          }
        }
      }
      if (dij_min < diB_min) {
        trk[i_min] = trk[i_min] + trk[j_min];
        trk[j_min].SetXYZ(0, 0, 0);
        n_particles_removed++;
      }
      if (dij_min > diB_min) {
        jet.push_back(trk[iB_min]);
        trk[iB_min].SetXYZ(0, 0, 0);
        n_particles_removed++;
      }
    } while (n_particles_removed < static_cast<int>(trk.size()));

    registryQC.fill(HIST("nJets_found"), static_cast<int>(jet.size()));

    // Jet Selection
    std::vector<int> isSelected;
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {
      isSelected.push_back(0);
    }

    int n_jets_selected(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {

      if ((TMath::Abs(jet[i].Eta()) + Rjet) > max_eta)
        continue;

      // Perpendicular cones
      TVector3 ue_axis1(0, 0, 0);
      TVector3 ue_axis2(0, 0, 0);
      get_perpendicular_axis(jet[i], ue_axis1, +1);
      get_perpendicular_axis(jet[i], ue_axis2, -1);
      ue1.push_back(ue_axis1);
      ue2.push_back(ue_axis2);

      double nPartJetPlusUE(0);
      double nPartJet(0);
      double nPartUE(0);
      double ptJetPlusUE(0);
      double ptJet(0);
      double ptUE(0);

      for (auto track : tracks) {

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        TVector3 sel_track(track.px(), track.py(), track.pz());

        double deltaEta_jet = sel_track.Eta() - jet[i].Eta();
        double deltaPhi_jet = GetDeltaPhi(sel_track.Phi(), jet[i].Phi());
        double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        double deltaEta_ue1 = sel_track.Eta() - ue_axis1.Eta();
        double deltaPhi_ue1 = GetDeltaPhi(sel_track.Phi(), ue_axis1.Phi());
        double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        double deltaEta_ue2 = sel_track.Eta() - ue_axis2.Eta();
        double deltaPhi_ue2 = GetDeltaPhi(sel_track.Phi(), ue_axis2.Phi());
        double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

        if (deltaR_jet < Rjet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_jet"), deltaEta_jet, deltaPhi_jet);
          nPartJetPlusUE++;
          ptJetPlusUE = ptJetPlusUE + sel_track.Pt();
        }
        if (deltaR_ue1 < Rjet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEta_ue1, deltaPhi_ue1);
          nPartUE++;
          ptUE = ptUE + sel_track.Pt();
        }
        if (deltaR_ue2 < Rjet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEta_ue2, deltaPhi_ue2);
          nPartUE++;
          ptUE = ptUE + sel_track.Pt();
        }
      }
      nPartJet = nPartJetPlusUE - 0.5 * nPartUE;
      ptJet = ptJetPlusUE - 0.5 * ptUE;
      registryQC.fill(HIST("NchJetPlusUE"), nPartJetPlusUE);
      registryQC.fill(HIST("NchJet"), nPartJet);
      registryQC.fill(HIST("NchUE"), nPartUE);
      registryQC.fill(HIST("sumPtJetPlusUE"), ptJetPlusUE);
      registryQC.fill(HIST("sumPtJet"), ptJet);
      registryQC.fill(HIST("sumPtUE"), ptUE);

      if (ptJet < min_jet_pt)
        continue;
      if (nPartJetPlusUE < min_nPartInJet)
        continue;
      n_jets_selected++;
      isSelected[i] = 1;
    }
    registryQC.fill(HIST("nJets_selected"), n_jets_selected);

    if (n_jets_selected == 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);
    //************************************************************************************************************************************

    // Overlaps
    int nOverlaps(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {
      if (isSelected[i] == 0)
        continue;

      for (int j = 0; j < static_cast<int>(jet.size()); j++) {
        if (isSelected[j] == 0 || i == j)
          continue;
        if (overlap(jet[i], ue1[j], Rjet) || overlap(jet[i], ue2[j], Rjet))
          nOverlaps++;
      }
    }
    registryQC.fill(HIST("jet_ue_overlaps"), n_jets_selected, nOverlaps);

    if (n_jets_selected > n_jets_per_event_max)
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    if (requireNoOverlap && nOverlaps > 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    //************************************************************************************************************************************

    for (int i = 0; i < static_cast<int>(jet.size()); i++) {

      if (isSelected[i] == 0)
        continue;

      for (auto track : tracks) {
        if (!passedTrackSelection(track))
          continue;
        if (require_PV_contributor && !(track.isPVContributor()))
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

        TVector3 particle_dir(track.px(), track.py(), track.pz());
        double deltaEta_jet = particle_dir.Eta() - jet[i].Eta();
        double deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), jet[i].Phi());
        double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        double deltaEta_ue1 = particle_dir.Eta() - ue1[i].Eta();
        double deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue1[i].Phi());
        double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        double deltaEta_ue2 = particle_dir.Eta() - ue2[i].Eta();
        double deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue2[i].Phi());
        double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

        // DCAxy Distributions of Antiprotons
        if (track.sign() < 0) { // only antiprotons
          if (isHighPurityAntiproton(track) && TMath::Abs(dcaz) < max_dcaz) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("antiproton_dca_jet"), pt, dcaxy);
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("antiproton_dca_ue"), pt, dcaxy);
            }
          }
        }
        // DCA Cuts
        if (TMath::Abs(dcaxy) > max_dcaxy)
          continue;
        if (TMath::Abs(dcaz) > max_dcaz)
          continue;

        // Jet
        if (deltaR_jet < Rjet) {

          if (track.sign() < 0) { // only antimatter
            // Antiproton
            if (pt < max_pt_for_nsigmaTPC)
              registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr);
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr);

            // Antideuteron
            if (pt < max_pt_for_nsigmaTPC)
              registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe);
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe);

            // Antihelium3
            registryData.fill(HIST("antihelium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
          }

          if (track.sign() > 0) { // only matter
            // Deuteron
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("deuteron_jet_tof"), pt, nsigmaTOFDe);

            // Helium3
            registryData.fill(HIST("helium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }

        // UE
        if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {

          if (track.sign() < 0) { // only antimatter
            // Antiproton
            if (pt < max_pt_for_nsigmaTPC)
              registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr);
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr);

            // Antideuteron
            if (pt < max_pt_for_nsigmaTPC)
              registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe);
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe);

            // Antihelium3
            registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
          }

          if (track.sign() > 0) { // only matter
            // Deuteron
            if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF())
              registryData.fill(HIST("deuteron_ue_tof"), pt, nsigmaTOFDe);

            // Helium3
            registryData.fill(HIST("helium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_jets, processData, "Process Data", true);

  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<MCTracks> perCollision = o2::aod::track::collisionId;

  void processEfficiency(o2::aod::McCollisions const& mcCollisions, SimCollisions const& collisions, MCTracks const& mcTracks, aod::McParticles const& mcParticles)
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
        if (particle.eta() < min_eta || particle.eta() > max_eta)
          continue;

        double w_antip_jet(1.0);
        double w_antip_ue(1.0);
        if (applyReweighting) {
          int ix = twod_weights_antip_jet->GetXaxis()->FindBin(particle.pt());
          int iy = twod_weights_antip_jet->GetYaxis()->FindBin(particle.eta());
          w_antip_jet = twod_weights_antip_jet->GetBinContent(ix, iy);
          w_antip_ue = twod_weights_antip_ue->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twod_weights_antip_jet->GetNbinsX()) {
            w_antip_jet = 1.0;
            w_antip_ue = 1.0;
          }
          if (iy == 0 || iy > twod_weights_antip_jet->GetNbinsY()) {
            w_antip_jet = 1.0;
            w_antip_ue = 1.0;
          }
        }

        if (particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_jet_gen"), particle.pt(), w_antip_jet);
          registryMC.fill(HIST("antiproton_ue_gen"), particle.pt(), w_antip_ue);
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

      if (TMath::Abs(collision.posZ()) > 10)
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

        double w_antip_jet(1.0);
        double w_antip_ue(1.0);
        if (applyReweighting) {
          int ix = twod_weights_antip_jet->GetXaxis()->FindBin(particle.pt());
          int iy = twod_weights_antip_jet->GetYaxis()->FindBin(particle.eta());
          w_antip_jet = twod_weights_antip_jet->GetBinContent(ix, iy);
          w_antip_ue = twod_weights_antip_ue->GetBinContent(ix, iy);

          // protection
          if (ix == 0 || ix > twod_weights_antip_jet->GetNbinsX()) {
            w_antip_jet = 1.0;
            w_antip_ue = 1.0;
          }
          if (iy == 0 || iy > twod_weights_antip_jet->GetNbinsY()) {
            w_antip_jet = 1.0;
            w_antip_ue = 1.0;
          }
        }

        // Antiproton
        if (particle.pdgCode() == -2212) {
          if (pt < max_pt_for_nsigmaTPC && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC) {
            registryMC.fill(HIST("antiproton_jet_rec_tpc"), pt, w_antip_jet);
            registryMC.fill(HIST("antiproton_ue_rec_tpc"), pt, w_antip_ue);
          }
          if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCPr > min_nsigmaTPC && nsigmaTPCPr < max_nsigmaTPC && track.hasTOF() && nsigmaTOFPr > min_nsigmaTOF && nsigmaTOFPr < max_nsigmaTOF) {
            registryMC.fill(HIST("antiproton_jet_rec_tof"), pt, w_antip_jet);
            registryMC.fill(HIST("antiproton_ue_rec_tof"), pt, w_antip_ue);
          }
        }

        // Antideuteron
        if (particle.pdgCode() == -1000010020) {
          if (pt < max_pt_for_nsigmaTPC && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC) {
            registryMC.fill(HIST("antideuteron_jet_rec_tpc"), pt);
            registryMC.fill(HIST("antideuteron_ue_rec_tpc"), pt);
          }
          if (pt >= min_pt_for_nsigmaTOF && nsigmaTPCDe > min_nsigmaTPC && nsigmaTPCDe < max_nsigmaTPC && track.hasTOF() && nsigmaTOFDe > min_nsigmaTOF && nsigmaTOFDe < max_nsigmaTOF) {
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
  PROCESS_SWITCH(nuclei_in_jets, processEfficiency, "process efficiency", false);

  void processSecondaryAntiprotons(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {
    for (const auto& collision : collisions) {

      registryMC.fill(HIST("number_of_events_mc"), 3.5);

      // Event Selection
      if (!collision.sel8())
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 4.5);

      if (TMath::Abs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 5.5);

      auto tracks_per_coll = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // List of Tracks
      std::vector<TVector3> trk;

      for (auto track : tracks_per_coll) {

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        TVector3 momentum(track.px(), track.py(), track.pz());
        trk.push_back(momentum);
      }

      // Anti-kt Jet Finder
      int n_particles_removed(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      do {
        double dij_min(1e+06), diB_min(1e+06);
        int i_min(0), j_min(0), iB_min(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) {
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diB_min) {
            diB_min = diB;
            iB_min = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) {
            if (trk[j].Mag() == 0)
              continue;
            double dij = calculate_dij(trk[i], trk[j], Rjet);
            if (dij < dij_min) {
              dij_min = dij;
              i_min = i;
              j_min = j;
            }
          }
        }
        if (dij_min < diB_min) {
          trk[i_min] = trk[i_min] + trk[j_min];
          trk[j_min].SetXYZ(0, 0, 0);
          n_particles_removed++;
        }
        if (dij_min > diB_min) {
          jet.push_back(trk[iB_min]);
          trk[iB_min].SetXYZ(0, 0, 0);
          n_particles_removed++;
        }
      } while (n_particles_removed < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {
        isSelected.push_back(0);
      }

      int n_jets_selected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {

        if ((TMath::Abs(jet[i].Eta()) + Rjet) > max_eta)
          continue;

        // Perpendicular cones
        TVector3 ue_axis1(0, 0, 0);
        TVector3 ue_axis2(0, 0, 0);
        get_perpendicular_axis(jet[i], ue_axis1, +1);
        get_perpendicular_axis(jet[i], ue_axis2, -1);
        ue1.push_back(ue_axis1);
        ue2.push_back(ue_axis2);

        double nPartJetPlusUE(0);
        double ptJetPlusUE(0);
        double ptJet(0);
        double ptUE(0);

        for (auto track : tracks_per_coll) {

          if (!passedTrackSelectionForJetReconstruction(track))
            continue;
          TVector3 sel_track(track.px(), track.py(), track.pz());

          double deltaEta_jet = sel_track.Eta() - jet[i].Eta();
          double deltaPhi_jet = GetDeltaPhi(sel_track.Phi(), jet[i].Phi());
          double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          double deltaEta_ue1 = sel_track.Eta() - ue_axis1.Eta();
          double deltaPhi_ue1 = GetDeltaPhi(sel_track.Phi(), ue_axis1.Phi());
          double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          double deltaEta_ue2 = sel_track.Eta() - ue_axis2.Eta();
          double deltaPhi_ue2 = GetDeltaPhi(sel_track.Phi(), ue_axis2.Phi());
          double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          if (deltaR_jet < Rjet) {
            nPartJetPlusUE++;
            ptJetPlusUE = ptJetPlusUE + sel_track.Pt();
          }
          if (deltaR_ue1 < Rjet) {
            ptUE = ptUE + sel_track.Pt();
          }
          if (deltaR_ue2 < Rjet) {
            ptUE = ptUE + sel_track.Pt();
          }
        }
        ptJet = ptJetPlusUE - 0.5 * ptUE;

        if (ptJet < min_jet_pt)
          continue;
        if (nPartJetPlusUE < min_nPartInJet)
          continue;
        n_jets_selected++;
        isSelected[i] = 1;
      }
      if (n_jets_selected == 0)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 6.5);

      for (int i = 0; i < static_cast<int>(jet.size()); i++) {

        if (isSelected[i] == 0)
          continue;

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
          if (!track.has_mcParticle())
            continue;
          const auto particle = track.mcParticle();
          if (particle.pdgCode() != -2212)
            continue;

          TVector3 particle_dir(track.px(), track.py(), track.pz());
          float deltaEta_jet = particle_dir.Eta() - jet[i].Eta();
          float deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), jet[i].Phi());
          float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          float deltaEta_ue1 = particle_dir.Eta() - ue1[i].Eta();
          float deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue1[i].Phi());
          float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          float deltaEta_ue2 = particle_dir.Eta() - ue2[i].Eta();
          float deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue2[i].Phi());
          float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          if (deltaR_jet < Rjet) {
            registryMC.fill(HIST("antiproton_all_jet"), track.pt());
            if (particle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_jet"), track.pt());
            }
          }
          if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
            registryMC.fill(HIST("antiproton_all_ue"), track.pt());
            if (particle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_ue"), track.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_jets, processSecondaryAntiprotons, "process secondary antiprotons", false);

  void processAntiprotonReweighting(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mcCollisions) {

      registryMC.fill(HIST("number_of_events_mc"), 7.5);

      // Selection on z_{vertex}
      if (TMath::Abs(mccollision.posZ()) > 10)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 8.5);

      // MC Particles per Collision
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // List of Tracks
      std::vector<TVector3> trk;

      for (auto& particle : mcParticles_per_coll) {
        if (particle.isPhysicalPrimary() && particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_eta_pt_pythia"), particle.pt(), particle.eta());
        }

        // Select Primary Particles
        double dx = particle.vx() - mccollision.posX();
        double dy = particle.vy() - mccollision.posY();
        double dz = particle.vz() - mccollision.posZ();
        double dcaxy = sqrt(dx * dx + dy * dy);
        double dcaz = TMath::Abs(dz);

        if (setDCAselectionPtDep) {
          if (dcaxy > (par0 + par1 / particle.pt()))
            continue;
          if (dcaz > (par0 + par1 / particle.pt()))
            continue;
        }
        if (!setDCAselectionPtDep) {
          if (dcaxy > max_dcaxy)
            continue;
          if (dcaz > max_dcaz)
            continue;
        }

        if (TMath::Abs(particle.eta()) > 0.8)
          continue;
        if (particle.pt() < 0.15)
          continue;

        // PDG Selection
        int pdg = TMath::Abs(particle.pdgCode());
        if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
          continue;

        TVector3 momentum(particle.px(), particle.py(), particle.pz());
        trk.push_back(momentum);
      }

      // Anti-kt Jet Finder
      int n_particles_removed(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      do {
        double dij_min(1e+06), diB_min(1e+06);
        int i_min(0), j_min(0), iB_min(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) {
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diB_min) {
            diB_min = diB;
            iB_min = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) {
            if (trk[j].Mag() == 0)
              continue;
            double dij = calculate_dij(trk[i], trk[j], Rjet);
            if (dij < dij_min) {
              dij_min = dij;
              i_min = i;
              j_min = j;
            }
          }
        }
        if (dij_min < diB_min) {
          trk[i_min] = trk[i_min] + trk[j_min];
          trk[j_min].SetXYZ(0, 0, 0);
          n_particles_removed++;
        }
        if (dij_min > diB_min) {
          jet.push_back(trk[iB_min]);
          trk[iB_min].SetXYZ(0, 0, 0);
          n_particles_removed++;
        }
      } while (n_particles_removed < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {
        isSelected.push_back(0);
      }

      int n_jets_selected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {

        if ((TMath::Abs(jet[i].Eta()) + Rjet) > max_eta)
          continue;

        // Perpendicular cones
        TVector3 ue_axis1(0, 0, 0);
        TVector3 ue_axis2(0, 0, 0);
        get_perpendicular_axis(jet[i], ue_axis1, +1);
        get_perpendicular_axis(jet[i], ue_axis2, -1);
        ue1.push_back(ue_axis1);
        ue2.push_back(ue_axis2);

        double nPartJetPlusUE(0);
        double ptJetPlusUE(0);
        double ptJet(0);
        double ptUE(0);

        for (auto& particle : mcParticles_per_coll) {

          // Select Primary Particles
          double dx = particle.vx() - mccollision.posX();
          double dy = particle.vy() - mccollision.posY();
          double dz = particle.vz() - mccollision.posZ();
          double dcaxy = sqrt(dx * dx + dy * dy);
          double dcaz = TMath::Abs(dz);

          if (setDCAselectionPtDep) {
            if (dcaxy > (par0 + par1 / particle.pt()))
              continue;
            if (dcaz > (par0 + par1 / particle.pt()))
              continue;
          }
          if (!setDCAselectionPtDep) {
            if (dcaxy > max_dcaxy)
              continue;
            if (dcaz > max_dcaz)
              continue;
          }

          if (TMath::Abs(particle.eta()) > 0.8)
            continue;
          if (particle.pt() < 0.15)
            continue;

          // PDG Selection
          int pdg = TMath::Abs(particle.pdgCode());
          if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
            continue;

          TVector3 sel_track(particle.px(), particle.py(), particle.pz());

          double deltaEta_jet = sel_track.Eta() - jet[i].Eta();
          double deltaPhi_jet = GetDeltaPhi(sel_track.Phi(), jet[i].Phi());
          double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          double deltaEta_ue1 = sel_track.Eta() - ue_axis1.Eta();
          double deltaPhi_ue1 = GetDeltaPhi(sel_track.Phi(), ue_axis1.Phi());
          double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          double deltaEta_ue2 = sel_track.Eta() - ue_axis2.Eta();
          double deltaPhi_ue2 = GetDeltaPhi(sel_track.Phi(), ue_axis2.Phi());
          double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          if (deltaR_jet < Rjet) {
            nPartJetPlusUE++;
            ptJetPlusUE = ptJetPlusUE + sel_track.Pt();
          }
          if (deltaR_ue1 < Rjet) {
            ptUE = ptUE + sel_track.Pt();
          }
          if (deltaR_ue2 < Rjet) {
            ptUE = ptUE + sel_track.Pt();
          }
        }
        ptJet = ptJetPlusUE - 0.5 * ptUE;

        if (ptJet < min_jet_pt)
          continue;
        if (nPartJetPlusUE < min_nPartInJet)
          continue;
        n_jets_selected++;
        isSelected[i] = 1;
      }
      if (n_jets_selected == 0)
        continue;

      for (int i = 0; i < static_cast<int>(jet.size()); i++) {

        if (isSelected[i] == 0)
          continue;

        // Generated Particles
        for (auto& particle : mcParticles_per_coll) {

          if (!particle.isPhysicalPrimary())
            continue;
          if (particle.pdgCode() != -2212)
            continue;

          TVector3 particle_dir(particle.px(), particle.py(), particle.pz());
          double deltaEta_jet = particle_dir.Eta() - jet[i].Eta();
          double deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), jet[i].Phi());
          double deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          double deltaEta_ue1 = particle_dir.Eta() - ue1[i].Eta();
          double deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue1[i].Phi());
          double deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          double deltaEta_ue2 = particle_dir.Eta() - ue2[i].Eta();
          double deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue2[i].Phi());
          double deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          if (deltaR_jet < Rjet) {
            registryMC.fill(HIST("antiproton_eta_pt_jet"), particle.pt(), particle.eta());
          }
          if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
            registryMC.fill(HIST("antiproton_eta_pt_ue"), particle.pt(), particle.eta());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nuclei_in_jets, processAntiprotonReweighting, "Process antiproton reweighting", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nuclei_in_jets>(cfgc)};
}
