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
/// \file nucleiInJets.cxx
///
/// \brief task for analysis of nuclei in jets
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since November 22, 2023

#include <vector>
#include <TList.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TGrid.h"
#include "TF1.h"
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
#include "Common/DataModel/PIDResponseITS.h"

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

using FullNucleiTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::McTrackLabels>;

struct NucleiInJets {

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
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet"};
  Configurable<double> rJet{"rJet", 0.3, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<int> minNparticlesInJet{"minNparticlesInJet", 2, "Minimum number of particles inside jet"};
  Configurable<int> nJetsPerEventMax{"nJetsPerEventMax", 1000, "Maximum number of jets per event"};
  Configurable<bool> requireNoOverlap{"requireNoOverlap", true, "require no overlap between jets and UE cones"};
  Configurable<int> nGhosts{"nGhosts", 1000, "number of ghost particles"};
  Configurable<double> alpha{"alpha", 1.0, "Alpha"};
  Configurable<double> averagePtUE{"averagePtUE", 0.1, "Average pt of UE"};
  Configurable<double> averagePtUEMC{"averagePtUEMC", 0.1, "Average pt of UE in MC"};

  // Track Parameters
  Configurable<double> par0{"par0", 0.00164, "par 0"};
  Configurable<double> par1{"par1", 0.00231, "par 1"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 80, "minimum number of TPC crossed pad rows"};
  Configurable<double> minTpcNcrossedRowsOverFindable{"minTpcNcrossedRowsOverFindable", 0.8, "crossed rows/findable"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.3, "minimum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  Configurable<double> minNsigmaTpc{"minNsigmaTpc", -3.0, "Minimum nsigma TPC"};
  Configurable<double> maxNsigmaTpc{"maxNsigmaTpc", +3.0, "Maximum nsigma TPC"};
  Configurable<double> minNsigmaTof{"minNsigmaTof", -3.0, "Minimum nsigma TOF"};
  Configurable<double> maxNsigmaTof{"maxNsigmaTof", +3.5, "Maximum nsigma TOF"};
  Configurable<double> maxPtForNsigmaTpc{"maxPtForNsigmaTpc", 2.0, "Maximum pt for TPC analysis"};
  Configurable<double> minPtForNsigmaTof{"minPtForNsigmaTof", 0.5, "Minimum pt for TOF analysis"};
  Configurable<bool> requirePvContributor{"requirePvContributor", true, "require that the track is a PV contributor"};
  Configurable<bool> setDCAselectionPtDep{"setDCAselectionPtDep", true, "require pt dependent selection"};
  Configurable<bool> applyReweighting{"applyReweighting", true, "apply reweighting"};
  Configurable<bool> applyItsPid{"applyItsPid", true, "apply ITS PID"};
  Configurable<double> ptMaxItsPid{"ptMaxItsPid", 1.0, "maximum pt for ITS PID"};
  Configurable<double> nSigmaItsMin{"nSigmaItsMin", -2.0, "nSigmaITS min"};
  Configurable<double> nSigmaItsMax{"nSigmaItsMax", +2.0, "nSigmaITS max"};
  Configurable<std::string> urlToCcdb{"urlToCcdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> pathToFile{"pathToFile", "", "path to file with reweighting"};
  Configurable<std::string> histoNameWeightAntipJet{"histoNameWeightAntipJet", "", "reweighting histogram: antip in jet"};
  Configurable<std::string> histoNameWeightAntipUe{"histoNameWeightAntipUe", "", "reweighting histogram: antip in ue"};

  TH2F* twoDweightsAntipJet;
  TH2F* twoDweightsAntipUe;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  void init(InitContext const&)
  {
    ccdb->setURL(urlToCcdb.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (applyReweighting) {
      getReweightingHistograms(ccdb, TString(pathToFile), TString(histoNameWeightAntipJet), TString(histoNameWeightAntipUe));
    } else {
      twoDweightsAntipJet = nullptr;
      twoDweightsAntipUe = nullptr;
    }

    // QC Histograms
    registryQC.add("deltaEtadeltaPhiJet", "deltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, PIHalf, "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhi_ue", "deltaEtadeltaPhi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, PIHalf, "#Delta#phi"}});
    registryQC.add("NchJetPlusUE", "NchJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchJet", "NchJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchUE", "NchUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("sumPtJetPlusUE", "sumPtJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtJet", "sumPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtUE", "sumPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtUE_MC", "sumPtUE_MC", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("nJets_found", "nJets_found", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("nJets_selected", "nJets_selected", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("dcaxy_vs_pt", "dcaxy_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    registryQC.add("dcaz_vs_pt", "dcaz_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    registryQC.add("jet_jet_overlaps", "jet_jet_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("jet_ue_overlaps", "jet_ue_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("ue_ue_overlaps", "ue_ue_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("tot_overlaps", "tot_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("hJetArea", "hJetArea", HistType::kTH1F, {{2000, 0, 2, "Area"}});
    registryQC.add("hError", "hError", HistType::kTH1F, {{5, 0, 5, "error"}});

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

    // Detector Response Matrix
    if (doprocessDetResponseMatrix) {
      registryMC.add("detectorResponseMatrix", "detectorResponseMatrix", HistType::kTH2F, {{5000, 0.0, 50.0, "#it{p}_{T}^{gen} (GeV/#it{c})"}, {5000, 0.0, 50.0, "#it{p}_{T}^{rec} (GeV/#it{c})"}});
    }
  }

  // ITS Hit
  template <typename T>
  bool hasITSHit(T const& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-Track Selection for Particles inside Jets
  template <typename T1>
  bool passedTrackSelectionForJetReconstruction(const T1& track)
  {
    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < 0.8)
      return false;
    if (track.tpcChi2NCl() > 4)
      return false;
    if (track.itsChi2NCl() > 36)
      return false;
    if (track.eta() < -0.8 || track.eta() > 0.8)
      return false;
    if (track.pt() < 0.15)
      return false;
    if (std::fabs(track.dcaXY()) > 0.25)
      return false;
    if (std::fabs(track.dcaZ()) > 2.0)
      return false;

    /*
    // pt-dependent selection
    if (setDCAselectionPtDep) {
      if (std::fabs(track.dcaXY()) > (par0 + par1 / track.pt()))
        return false;
      if (std::fabs(track.dcaZ()) > (par0 + par1 / track.pt()))
        return false;
    }

    // standard selection
    if (!setDCAselectionPtDep) {
      if (std::fabs(track.dcaXY()) > maxDcaxy)
        return false;
      if (std::fabs(track.dcaZ()) > maxDcaz)
        return false;
    }
    */

    return true;
  }

  // Single-Track Selection
  template <typename T2>
  bool passedTrackSelection(const T2& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minItsNclusters)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows)
      return false;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < minTpcNcrossedRowsOverFindable)
      return false;
    if (track.tpcChi2NCl() > maxChiSquareTpc)
      return false;
    if (track.itsChi2NCl() > maxChiSquareIts)
      return false;
    if (track.eta() < minEta || track.eta() > maxEta)
      return false;
    if (track.pt() < minPt)
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

    if (pt < 0.5 && std::fabs(nsigmaTPCPr) < 2.0)
      return true;
    if (pt >= 0.5 && std::fabs(nsigmaTPCPr) < 2.0 && track.hasTOF() && std::fabs(nsigmaTOFPr) < 2.0)
      return true;
    return false;
  }

  // Minimum
  double minimumValue(double x1, double x2)
  {
    double xMin(x1);
    if (x1 < x2)
      xMin = x1;
    if (x1 >= x2)
      xMin = x2;

    return xMin;
  }

  // Deltaphi
  double getDeltaPhi(double a1, double a2)
  {
    double deltaPhi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = std::fabs(phi1 - phi2);

    if (diff <= PI)
      deltaPhi = diff;
    if (diff > PI)
      deltaPhi = TwoPI - diff;

    return deltaPhi;
  }

  void getPerpendicularAxis(TVector3 p, TVector3& u, double sign)
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
      ux = sign * std::sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Protection 2
    if (py == 0 && px != 0) {

      ux = -(pz * pz) / px;
      uy = sign * std::sqrt(px * px - (pz * pz * pz * pz) / (px * px));
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
    ux = (-b + sign * std::sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }

  double calculateDij(TVector3 t1, TVector3 t2, double R)
  {
    double distanceJet(0);
    double x1 = 1.0 / (t1.Pt() * t1.Pt());
    double x2 = 1.0 / (t2.Pt() * t2.Pt());
    double deltaEta = t1.Eta() - t2.Eta();
    double deltaPhi = getDeltaPhi(t1.Phi(), t2.Phi());
    double min = minimumValue(x1, x2);
    double deltaSquare = deltaEta * deltaEta + deltaPhi * deltaPhi;
    distanceJet = min * deltaSquare / (R * R);
    return distanceJet;
  }

  bool overlap(TVector3 v1, TVector3 v2, double R)
  {
    double dx = v1.Eta() - v2.Eta();
    double dy = getDeltaPhi(v1.Phi(), v2.Phi());
    double d = std::sqrt(dx * dx + dy * dy);
    if (d < 2.0 * R)
      return true;
    return false;
  }

  double getCorrectedPt(double ptRec)
  {
    // to be developed
    return ptRec;
  }

  void getReweightingHistograms(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, TString filepath, TString histname_antip_jet, TString histname_antip_ue)
  {
    TList* l = ccdbObj->get<TList>(filepath.Data());
    if (!l) {
      LOGP(error, "Could not open the file {}", Form("%s", filepath.Data()));
      return;
    }
    twoDweightsAntipJet = static_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname_antip_jet.Data())));
    if (!twoDweightsAntipJet) {
      LOGP(error, "Could not open histogram {}", Form("%s_antiproton", histname_antip_jet.Data()));
      return;
    }
    twoDweightsAntipUe = static_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname_antip_ue.Data())));
    if (!twoDweightsAntipUe) {
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
    if (std::fabs(collision.posZ()) > zVtx)
      return;

    // Event Counter: after z-vertex cut
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // ITS Response
    o2::aod::ITSResponse itsResponse;

    // List of Tracks
    std::vector<TVector3> trk;
    std::vector<int> ntrk;

    for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]

      if (!passedTrackSelectionForJetReconstruction(track))
        continue;
      registryQC.fill(HIST("dcaxy_vs_pt"), track.pt(), track.dcaXY());
      registryQC.fill(HIST("dcaz_vs_pt"), track.pt(), track.dcaZ());

      TVector3 momentum(track.px(), track.py(), track.pz());
      trk.push_back(momentum);
      ntrk.push_back(1);
    }

    // Reject Empty Events
    if (static_cast<int>(trk.size()) < 1)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Anti-kt Jet Finder
    int nParticlesRemoved(0);
    std::vector<TVector3> jet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;
    std::vector<int> nParticlesInjet;

    do {
      double dijMin(1e+300), diBmin(1e+300);
      int iMin(0), jMin(0), iBmin(0);
      for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
        if (trk[i].Mag() == 0)
          continue;
        double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
        if (diB < diBmin) {
          diBmin = diB;
          iBmin = i;
        }
        for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[j].Mag() == 0)
            continue;
          double dij = calculateDij(trk[i], trk[j], rJet);
          if (dij < dijMin) {
            dijMin = dij;
            iMin = i;
            jMin = j;
          }
        }
      }
      if (dijMin < diBmin) {
        trk[iMin] = trk[iMin] + trk[jMin];
        ntrk[iMin] = ntrk[iMin] + ntrk[jMin];
        trk[jMin].SetXYZ(0, 0, 0);
        ntrk[jMin] = 0;
        nParticlesRemoved++;
      }
      if (dijMin > diBmin) {
        jet.push_back(trk[iBmin]);
        nParticlesInjet.push_back(ntrk[iBmin]);
        trk[iBmin].SetXYZ(0, 0, 0);
        nParticlesRemoved++;
      }
    } while (nParticlesRemoved < static_cast<int>(trk.size()));

    registryQC.fill(HIST("nJets_found"), static_cast<int>(jet.size()));

    // Jet Selection
    std::vector<int> isSelected;
    int nJetsSelected(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

      // Initialization
      isSelected.push_back(0);

      // Jet fully contained inside acceptance
      if ((std::fabs(jet[i].Eta()) + rJet) > (maxEta - 0.5))
        continue;
      if (nParticlesInjet[i] < minNparticlesInJet)
        continue;

      // Perpendicular cones
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jet[i], ueAxis1, +1);
      getPerpendicularAxis(jet[i], ueAxis2, -1);
      ue1.push_back(ueAxis1);
      ue2.push_back(ueAxis2);

      double ptUE(0);
      for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        TVector3 selectedTrack(track.px(), track.py(), track.pz());

        double deltaEtaUe1 = selectedTrack.Eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(selectedTrack.Phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = selectedTrack.Eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(selectedTrack.Phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        if (deltaRUe1 < alpha * rJet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEtaUe1, deltaPhiUe1);
          ptUE = ptUE + selectedTrack.Pt();
        }
        if (deltaRUe2 < alpha * rJet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEtaUe2, deltaPhiUe2);
          ptUE = ptUE + selectedTrack.Pt();
        }
      }
      registryQC.fill(HIST("sumPtUE"), 0.5 * ptUE);
      registryQC.fill(HIST("NchJetPlusUE"), nParticlesInjet[i]);

      double ptJetRec = jet[i].Pt() - averagePtUE;
      double ptJetCorr = getCorrectedPt(ptJetRec);

      if (ptJetCorr < minJetPt)
        continue;

      nJetsSelected++;
      isSelected[i] = 1;
    }
    registryQC.fill(HIST("nJets_selected"), nJetsSelected);

    if (nJetsSelected == 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Overlaps
    int nOverlapsJetJet(0);
    int nOverlapsJetUe(0);
    int nOverlapsUeUe(0);
    int nOverlapsTot(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {
      if (isSelected[i] == 0)
        continue;

      for (int j = (i + 1); j < static_cast<int>(jet.size()); j++) {
        if (isSelected[j] == 0)
          continue;
        if (overlap(jet[i], jet[j], rJet))
          nOverlapsJetJet++;
        if (overlap(jet[i], ue1[j], rJet) || overlap(jet[i], ue2[j], rJet))
          nOverlapsJetUe++;
        if (overlap(ue1[i], ue1[j], rJet) || overlap(ue1[i], ue2[j], rJet) || overlap(ue2[i], ue2[j], rJet))
          nOverlapsUeUe++;
      }
    }
    nOverlapsTot = nOverlapsJetJet + nOverlapsJetUe + nOverlapsUeUe;
    registryQC.fill(HIST("jet_jet_overlaps"), nJetsSelected, nOverlapsJetJet);
    registryQC.fill(HIST("jet_ue_overlaps"), nJetsSelected, nOverlapsJetUe);
    registryQC.fill(HIST("ue_ue_overlaps"), nJetsSelected, nOverlapsUeUe);
    registryQC.fill(HIST("tot_overlaps"), nJetsSelected, nOverlapsTot);

    if (nJetsSelected > nJetsPerEventMax)
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    if (requireNoOverlap && nOverlapsTot > 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 6.5);

    //**************************************************************************************************************

    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

      if (isSelected[i] == 0)
        continue;

      for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]
        if (!passedTrackSelection(track))
          continue;
        if (requirePvContributor && !(track.isPVContributor()))
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

        // Selection on <ITS Cluster size>
        bool passedItsPid = false;
        if (itsResponse.nSigmaITS<o2::track::PID::Proton>(track) > nSigmaItsMin && itsResponse.nSigmaITS<o2::track::PID::Proton>(track) < nSigmaItsMax) {
          passedItsPid = true;
        }

        bool passedItsPidSelection = true;
        if (applyItsPid && pt < ptMaxItsPid && (!passedItsPid))
          passedItsPidSelection = false;

        TVector3 particleDirection(track.px(), track.py(), track.pz());
        double deltaEtaJet = particleDirection.Eta() - jet[i].Eta();
        double deltaPhiJet = getDeltaPhi(particleDirection.Phi(), jet[i].Phi());
        double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
        double deltaEtaUe1 = particleDirection.Eta() - ue1[i].Eta();
        double deltaPhiUe1 = getDeltaPhi(particleDirection.Phi(), ue1[i].Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = particleDirection.Eta() - ue2[i].Eta();
        double deltaPhiUe2 = getDeltaPhi(particleDirection.Phi(), ue2[i].Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        // DCAxy Distributions of Antiprotons
        if (track.sign() < 0) { // only antiprotons
          if (isHighPurityAntiproton(track) && std::fabs(dcaz) < maxDcaz) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("antiproton_dca_jet"), pt, dcaxy);
            }
            if (deltaRUe1 < rJet || deltaRUe2 < rJet) {
              registryData.fill(HIST("antiproton_dca_ue"), pt, dcaxy);
            }
          }
        }
        // DCA Cuts
        if (std::fabs(dcaxy) > maxDcaxy)
          continue;
        if (std::fabs(dcaz) > maxDcaz)
          continue;

        // Jet
        if (deltaRjet < rJet) {

          if (track.sign() < 0) { // only antimatter
            // Antiproton
            if (passedItsPidSelection) {
              if (pt < maxPtForNsigmaTpc)
                registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr);
              if (pt >= minPtForNsigmaTof && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && track.hasTOF())
                registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr);
            }

            // Antideuteron
            if (pt < maxPtForNsigmaTpc)
              registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe);
            if (pt >= minPtForNsigmaTof && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe);

            // Antihelium3
            registryData.fill(HIST("antihelium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
          }

          if (track.sign() > 0) { // only matter
            // Deuteron
            if (pt >= minPtForNsigmaTof && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("deuteron_jet_tof"), pt, nsigmaTOFDe);

            // Helium3
            registryData.fill(HIST("helium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }

        // UE
        if (deltaRUe1 < rJet || deltaRUe2 < rJet) {

          if (track.sign() < 0) { // only antimatter
            // Antiproton
            if (passedItsPidSelection) {
              if (pt < maxPtForNsigmaTpc)
                registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr);
              if (pt >= minPtForNsigmaTof && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && track.hasTOF())
                registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr);
            }
            // Antideuteron
            if (pt < maxPtForNsigmaTpc)
              registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe);
            if (pt >= minPtForNsigmaTof && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe);

            // Antihelium3
            registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
          }

          if (track.sign() > 0) { // only matter
            // Deuteron
            if (pt >= minPtForNsigmaTof && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("deuteron_ue_tof"), pt, nsigmaTOFDe);

            // Helium3
            registryData.fill(HIST("helium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiInJets, processData, "Process Data", true);

  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<MCTracks> perCollision = o2::aod::track::collisionId;

  void processEfficiency(o2::aod::McCollisions const& mcCollisions, SimCollisions const& collisions, MCTracks const& mcTracks, aod::McParticles const& mcParticles)
  {
    // Generated Events
    for (const auto& mccollision : mcCollisions) { // o2-linter: disable=[const-ref-in-for-loop]

      registryMC.fill(HIST("number_of_events_mc"), 0.5);
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      for (auto& particle : mcParticlesPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!particle.isPhysicalPrimary())
          continue;
        if ((particle.pdgCode() != -2212) && (particle.pdgCode() != -1000010020) && (particle.pdgCode() != -1000020030))
          continue;
        if (particle.eta() < minEta || particle.eta() > maxEta)
          continue;

        double wAntipJet(1.0);
        double wAntipUe(1.0);
        if (applyReweighting) {
          int ix = twoDweightsAntipJet->GetXaxis()->FindBin(particle.pt());
          int iy = twoDweightsAntipJet->GetYaxis()->FindBin(particle.eta());
          wAntipJet = twoDweightsAntipJet->GetBinContent(ix, iy);
          wAntipUe = twoDweightsAntipUe->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twoDweightsAntipJet->GetNbinsX()) {
            wAntipJet = 1.0;
            wAntipUe = 1.0;
          }
          if (iy == 0 || iy > twoDweightsAntipJet->GetNbinsY()) {
            wAntipJet = 1.0;
            wAntipUe = 1.0;
          }
        }

        if (particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_jet_gen"), particle.pt(), wAntipJet);
          registryMC.fill(HIST("antiproton_ue_gen"), particle.pt(), wAntipUe);
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
    for (const auto& collision : collisions) { // o2-linter: disable=[const-ref-in-for-loop]

      registryMC.fill(HIST("number_of_events_mc"), 1.5);

      // Event Selection
      if (!collision.sel8())
        continue;

      if (std::fabs(collision.posZ()) > 10)
        continue;

      // Event Counter (after event sel)
      registryMC.fill(HIST("number_of_events_mc"), 2.5);

      auto tracksPerColl = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // Reconstructed Tracks
      for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;

        const auto particle = track.mcParticle();
        if ((particle.pdgCode() != -2212) && (particle.pdgCode() != -1000010020) && (particle.pdgCode() != -1000020030))
          continue;

        // Track Selection
        if (!passedTrackSelection(track))
          continue;
        if (requirePvContributor && !(track.isPVContributor()))
          continue;

        // Variables
        float nsigmaTPCPr = track.tpcNSigmaPr();
        float nsigmaTOFPr = track.tofNSigmaPr();
        float nsigmaTPCDe = track.tpcNSigmaDe();
        float nsigmaTOFDe = track.tofNSigmaDe();
        float nsigmaTPCHe = track.tpcNSigmaHe();
        float pt = track.pt();

        // DCA Templates
        if (particle.pdgCode() == -2212 && particle.isPhysicalPrimary() && std::fabs(track.dcaZ()) < maxDcaz)
          registryMC.fill(HIST("antiproton_dca_prim"), pt, track.dcaXY());

        if (particle.pdgCode() == -2212 && (!particle.isPhysicalPrimary()) && std::fabs(track.dcaZ()) < maxDcaz)
          registryMC.fill(HIST("antiproton_dca_sec"), pt, track.dcaXY());

        // DCA Cuts
        if (std::fabs(track.dcaXY()) > maxDcaxy)
          continue;
        if (std::fabs(track.dcaZ()) > maxDcaz)
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

        double wAntipJet(1.0);
        double wAntipUe(1.0);
        if (applyReweighting) {
          int ix = twoDweightsAntipJet->GetXaxis()->FindBin(particle.pt());
          int iy = twoDweightsAntipJet->GetYaxis()->FindBin(particle.eta());
          wAntipJet = twoDweightsAntipJet->GetBinContent(ix, iy);
          wAntipUe = twoDweightsAntipUe->GetBinContent(ix, iy);

          // protection
          if (ix == 0 || ix > twoDweightsAntipJet->GetNbinsX()) {
            wAntipJet = 1.0;
            wAntipUe = 1.0;
          }
          if (iy == 0 || iy > twoDweightsAntipJet->GetNbinsY()) {
            wAntipJet = 1.0;
            wAntipUe = 1.0;
          }
        }

        // Antiproton
        if (particle.pdgCode() == -2212) {
          if (pt < maxPtForNsigmaTpc && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_jet_rec_tpc"), pt, wAntipJet);
            registryMC.fill(HIST("antiproton_ue_rec_tpc"), pt, wAntipUe);
          }
          if (pt >= minPtForNsigmaTof && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
            registryMC.fill(HIST("antiproton_jet_rec_tof"), pt, wAntipJet);
            registryMC.fill(HIST("antiproton_ue_rec_tof"), pt, wAntipUe);
          }
        }

        // Antideuteron
        if (particle.pdgCode() == -1000010020) {
          if (pt < maxPtForNsigmaTpc && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc) {
            registryMC.fill(HIST("antideuteron_jet_rec_tpc"), pt);
            registryMC.fill(HIST("antideuteron_ue_rec_tpc"), pt);
          }
          if (pt >= minPtForNsigmaTof && nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF() && nsigmaTOFDe > minNsigmaTof && nsigmaTOFDe < maxNsigmaTof) {
            registryMC.fill(HIST("antideuteron_jet_rec_tof"), pt);
            registryMC.fill(HIST("antideuteron_ue_rec_tof"), pt);
          }
        }

        // Antihelium-3
        if (particle.pdgCode() == -1000020030) {
          if (nsigmaTPCHe > minNsigmaTpc && nsigmaTPCHe < maxNsigmaTpc) {
            registryMC.fill(HIST("antihelium3_jet_rec_tpc"), 2.0 * pt);
            registryMC.fill(HIST("antihelium3_ue_rec_tpc"), 2.0 * pt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiInJets, processEfficiency, "process efficiency", false);

  void processSecondaryAntiprotons(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {
    for (const auto& collision : collisions) { // o2-linter: disable=[const-ref-in-for-loop]

      registryMC.fill(HIST("number_of_events_mc"), 3.5);

      // Event Selection
      if (!collision.sel8())
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 4.5);

      if (std::fabs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 5.5);

      auto tracksPerColl = mcTracks.sliceBy(perCollision, collision.globalIndex());

      // List of Tracks
      std::vector<TVector3> trk;
      std::vector<TVector3> part;
      std::vector<int> ntrk;

      for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        if (!track.has_mcParticle())
          continue;

        const auto particle = track.mcParticle();
        TVector3 pRec(track.px(), track.py(), track.pz());
        TVector3 pGen(particle.px(), particle.py(), particle.pz());
        trk.push_back(pRec);
        part.push_back(pGen);
        ntrk.push_back(1);
      }

      // Reject Empty Events
      if (static_cast<int>(trk.size()) < 1)
        continue;

      // Anti-kt Jet Finder
      int nParticlesRemoved(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> jetGen;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;
      std::vector<int> nParticlesInjet;

      do {
        double dijMin(1e+300), diBmin(1e+300);
        int iMin(0), jMin(0), iBmin(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diBmin) {
            diBmin = diB;
            iBmin = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
            if (trk[j].Mag() == 0)
              continue;
            double dij = calculateDij(trk[i], trk[j], rJet);
            if (dij < dijMin) {
              dijMin = dij;
              iMin = i;
              jMin = j;
            }
          }
        }
        if (dijMin < diBmin) {
          trk[iMin] = trk[iMin] + trk[jMin];
          ntrk[iMin] = ntrk[iMin] + ntrk[jMin];
          part[iMin] = part[iMin] + part[jMin];
          trk[jMin].SetXYZ(0, 0, 0);
          part[jMin].SetXYZ(0, 0, 0);
          ntrk[jMin] = 0;
          nParticlesRemoved++;
        }
        if (dijMin > diBmin) {
          jet.push_back(trk[iBmin]);
          jetGen.push_back(part[iBmin]);
          nParticlesInjet.push_back(ntrk[iBmin]);
          trk[iBmin].SetXYZ(0, 0, 0);
          part[iBmin].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
      } while (nParticlesRemoved < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      int nJetsSelected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        // Initialization
        isSelected.push_back(0);

        // Jet fully contained inside acceptance
        if ((std::fabs(jet[i].Eta()) + rJet) > (maxEta - 0.5))
          continue;
        if (nParticlesInjet[i] < minNparticlesInJet)
          continue;

        // Perpendicular cones
        TVector3 ueAxis1(0, 0, 0);
        TVector3 ueAxis2(0, 0, 0);
        getPerpendicularAxis(jet[i], ueAxis1, +1);
        getPerpendicularAxis(jet[i], ueAxis2, -1);
        ue1.push_back(ueAxis1);
        ue2.push_back(ueAxis2);

        double ptUE(0);
        for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

          if (!passedTrackSelectionForJetReconstruction(track))
            continue;
          TVector3 selectedTrack(track.px(), track.py(), track.pz());
          if (!track.has_mcParticle())
            continue;

          const auto particle = track.mcParticle();
          double deltaEtaUe1 = selectedTrack.Eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(selectedTrack.Phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = selectedTrack.Eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(selectedTrack.Phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          if ((deltaRUe1 < alpha * rJet) || (deltaRUe2 < alpha * rJet)) {
            ptUE = ptUE + particle.pt();
          }
        }
        registryQC.fill(HIST("sumPtUE_MC"), 0.5 * ptUE);

        double ptJetCorr = jetGen[i].Pt() - averagePtUEMC;

        if (ptJetCorr < minJetPt)
          continue;

        nJetsSelected++;
        isSelected[i] = 1;
      }
      if (nJetsSelected == 0)
        continue;

      // Overlaps
      int nOverlapsJetJet(0);
      int nOverlapsJetUe(0);
      int nOverlapsUeUe(0);
      int nOverlapsTot(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {
        if (isSelected[i] == 0)
          continue;

        for (int j = (i + 1); j < static_cast<int>(jet.size()); j++) {
          if (isSelected[j] == 0)
            continue;
          if (overlap(jet[i], jet[j], rJet))
            nOverlapsJetJet++;
          if (overlap(jet[i], ue1[j], rJet) || overlap(jet[i], ue2[j], rJet))
            nOverlapsJetUe++;
          if (overlap(ue1[i], ue1[j], rJet) || overlap(ue1[i], ue2[j], rJet) || overlap(ue2[i], ue2[j], rJet))
            nOverlapsUeUe++;
        }
      }
      nOverlapsTot = nOverlapsJetJet + nOverlapsJetUe + nOverlapsUeUe;
      if (nJetsSelected > nJetsPerEventMax)
        continue;
      if (requireNoOverlap && nOverlapsTot > 0)
        continue;

      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        if (isSelected[i] == 0)
          continue;

        for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]
          if (!passedTrackSelection(track))
            continue;
          if (requirePvContributor && !(track.isPVContributor()))
            continue;
          if (track.sign() > 0)
            continue;
          if (std::fabs(track.dcaXY()) > maxDcaxy)
            continue;
          if (std::fabs(track.dcaZ()) > maxDcaz)
            continue;
          if (!track.has_mcParticle())
            continue;
          const auto particle = track.mcParticle();
          if (particle.pdgCode() != -2212)
            continue;

          TVector3 particleDirection(track.px(), track.py(), track.pz());
          float deltaEtaJet = particleDirection.Eta() - jet[i].Eta();
          float deltaPhiJet = getDeltaPhi(particleDirection.Phi(), jet[i].Phi());
          float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          float deltaEtaUe1 = particleDirection.Eta() - ue1[i].Eta();
          float deltaPhiUe1 = getDeltaPhi(particleDirection.Phi(), ue1[i].Phi());
          float deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          float deltaEtaUe2 = particleDirection.Eta() - ue2[i].Eta();
          float deltaPhiUe2 = getDeltaPhi(particleDirection.Phi(), ue2[i].Phi());
          float deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          if (deltaRjet < rJet) {
            registryMC.fill(HIST("antiproton_all_jet"), track.pt());
            if (particle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_jet"), track.pt());
            }
          }
          if (deltaRUe1 < rJet || deltaRUe2 < rJet) {
            registryMC.fill(HIST("antiproton_all_ue"), track.pt());
            if (particle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_ue"), track.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiInJets, processSecondaryAntiprotons, "process secondary antiprotons", false);

  void processAntiprotonReweighting(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mcCollisions) { // o2-linter: disable=[const-ref-in-for-loop]

      registryMC.fill(HIST("number_of_events_mc"), 7.5);

      // Selection on z_{vertex}
      if (std::fabs(mccollision.posZ()) > 10)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 8.5);

      // MC Particles per Collision
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // List of Tracks
      std::vector<TVector3> trk;
      std::vector<int> ntrk;

      for (auto& particle : mcParticlesPerColl) { // o2-linter: disable=[const-ref-in-for-loop]
        if (particle.isPhysicalPrimary() && particle.pdgCode() == -2212) {
          registryMC.fill(HIST("antiproton_eta_pt_pythia"), particle.pt(), particle.eta());
        }

        // Select Primary Particles
        double dx = particle.vx() - mccollision.posX();
        double dy = particle.vy() - mccollision.posY();
        double dz = particle.vz() - mccollision.posZ();
        double dcaxy = std::sqrt(dx * dx + dy * dy);
        double dcaz = std::fabs(dz);
        if (dcaxy > 0.25)
          continue;
        if (dcaz > 2.0)
          continue;
        if (std::fabs(particle.eta()) > 0.8)
          continue;
        if (particle.pt() < 0.15)
          continue;

        // PDG Selection
        int pdg = std::fabs(particle.pdgCode());
        if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
          continue;

        TVector3 momentum(particle.px(), particle.py(), particle.pz());
        trk.push_back(momentum);
        ntrk.push_back(1);
      }

      // Anti-kt Jet Finder
      int nParticlesRemoved(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;
      std::vector<int> nParticlesInjet;

      do {
        double dijMin(1e+300), diBmin(1e+300);
        int iMin(0), jMin(0), iBmin(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diBmin) {
            diBmin = diB;
            iBmin = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
            if (trk[j].Mag() == 0)
              continue;
            double dij = calculateDij(trk[i], trk[j], rJet);
            if (dij < dijMin) {
              dijMin = dij;
              iMin = i;
              jMin = j;
            }
          }
        }
        if (dijMin < diBmin) {
          trk[iMin] = trk[iMin] + trk[jMin];
          ntrk[iMin] = ntrk[iMin] + ntrk[jMin];
          trk[jMin].SetXYZ(0, 0, 0);
          ntrk[jMin] = 0;
          nParticlesRemoved++;
        }
        if (dijMin > diBmin) {
          jet.push_back(trk[iBmin]);
          nParticlesInjet.push_back(ntrk[iBmin]);
          trk[iBmin].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
      } while (nParticlesRemoved < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      int nJetsSelected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        // Initialization
        isSelected.push_back(0);

        // Jet fully contained inside acceptance
        if ((std::fabs(jet[i].Eta()) + rJet) > (maxEta - 0.5))
          continue;
        if (nParticlesInjet[i] < minNparticlesInJet)
          continue;

        // Perpendicular cones
        TVector3 ueAxis1(0, 0, 0);
        TVector3 ueAxis2(0, 0, 0);
        getPerpendicularAxis(jet[i], ueAxis1, +1);
        getPerpendicularAxis(jet[i], ueAxis2, -1);
        ue1.push_back(ueAxis1);
        ue2.push_back(ueAxis2);

        double ptJetCorr = jet[i].Pt() - averagePtUEMC;
        if (ptJetCorr < minJetPt)
          continue;

        nJetsSelected++;
        isSelected[i] = 1;
      }
      if (nJetsSelected == 0)
        continue;

      // Overlaps
      int nOverlapsJetJet(0);
      int nOverlapsJetUe(0);
      int nOverlapsUeUe(0);
      int nOverlapsTot(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) {
        if (isSelected[i] == 0)
          continue;

        for (int j = (i + 1); j < static_cast<int>(jet.size()); j++) {
          if (isSelected[j] == 0)
            continue;
          if (overlap(jet[i], jet[j], rJet))
            nOverlapsJetJet++;
          if (overlap(jet[i], ue1[j], rJet) || overlap(jet[i], ue2[j], rJet))
            nOverlapsJetUe++;
          if (overlap(ue1[i], ue1[j], rJet) || overlap(ue1[i], ue2[j], rJet) || overlap(ue2[i], ue2[j], rJet))
            nOverlapsUeUe++;
        }
      }
      nOverlapsTot = nOverlapsJetJet + nOverlapsJetUe + nOverlapsUeUe;
      if (nJetsSelected > nJetsPerEventMax)
        continue;
      if (requireNoOverlap && nOverlapsTot > 0)
        continue;

      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        if (isSelected[i] == 0)
          continue;

        // Generated Particles
        for (auto& particle : mcParticlesPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

          if (!particle.isPhysicalPrimary())
            continue;
          if (particle.pdgCode() != -2212)
            continue;

          TVector3 particleDirection(particle.px(), particle.py(), particle.pz());
          double deltaEtaJet = particleDirection.Eta() - jet[i].Eta();
          double deltaPhiJet = getDeltaPhi(particleDirection.Phi(), jet[i].Phi());
          double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = particleDirection.Eta() - ue1[i].Eta();
          double deltaPhiUe1 = getDeltaPhi(particleDirection.Phi(), ue1[i].Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = particleDirection.Eta() - ue2[i].Eta();
          double deltaPhiUe2 = getDeltaPhi(particleDirection.Phi(), ue2[i].Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          if (deltaRjet < rJet) {
            registryMC.fill(HIST("antiproton_eta_pt_jet"), particle.pt(), particle.eta());
          }
          if (deltaRUe1 < rJet || deltaRUe2 < rJet) {
            registryMC.fill(HIST("antiproton_eta_pt_ue"), particle.pt(), particle.eta());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiInJets, processAntiprotonReweighting, "Process antiproton reweighting", false);

  void processGhosts(SelectedCollisions::iterator const& collision, FullNucleiTracks const& tracks)
  {
    // Event Selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Track Selection for Jet Finding
    std::vector<TVector3> trk;
    for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]

      if (!passedTrackSelectionForJetReconstruction(track))
        continue;
      TVector3 momentum(track.px(), track.py(), track.pz());
      trk.push_back(momentum);
    }
    // Reject Empty Events
    if (static_cast<int>(trk.size()) < 1)
      return;

    // Generate Ghosts
    for (int i = 0; i < nGhosts; i++) { // o2-linter: disable=[const-ref-in-for-loop]

      double eta = gRandom->Uniform(-0.8, 0.8);
      double phi = gRandom->Uniform(0.0, TwoPI);
      double pt = 1e-100;
      TVector3 ghost;
      ghost.SetPtEtaPhi(pt, eta, phi);
      trk.push_back(ghost);
    }

    // Anti-kt Jet Finder
    int nParticlesRemoved(0);
    std::vector<TVector3> jet;
    std::vector<double> jetArea;

    do {
      double dijMin(1e+300), diBmin(1e+300);
      int iMin(0), jMin(0), iBmin(0);
      int nGhostsInJet(0);
      for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
        if (trk[i].Mag() == 0)
          continue;
        double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
        if (diB < diBmin) {
          diBmin = diB;
          iBmin = i;
        }
        for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[j].Mag() == 0)
            continue;
          double dij = calculateDij(trk[i], trk[j], rJet);
          if (dij < dijMin) {
            dijMin = dij;
            iMin = i;
            jMin = j;
          }
        }
      }
      if (dijMin < diBmin) {
        if (trk[iMin].Pt() == 1e-100)
          nGhostsInJet++;
        if (trk[jMin].Pt() == 1e-100)
          nGhostsInJet++;
        trk[iMin] = trk[iMin] + trk[jMin];
        trk[jMin].SetXYZ(0, 0, 0);
        nParticlesRemoved++;
      }
      if (dijMin > diBmin) {
        double area = (static_cast<double>(nGhostsInJet) / static_cast<double>(nGhosts)) * TwoPI * 1.6;
        double alphaJet = area / (PI * rJet * rJet);
        jetArea.push_back(alphaJet);
        jet.push_back(trk[iBmin]);
        trk[iBmin].SetXYZ(0, 0, 0);
        nParticlesRemoved++;
      }
      if (dijMin == diBmin) {
        registryQC.fill(HIST("hError"), 0.5);
        nParticlesRemoved = static_cast<int>(trk.size());
      }
    } while (nParticlesRemoved < static_cast<int>(trk.size()));

    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

      if ((std::fabs(jet[i].Eta()) + rJet) > (maxEta - 0.5))
        continue;

      double ptJetRec = jet[i].Pt() - averagePtUE;
      double ptJetCorr = getCorrectedPt(ptJetRec);
      if (ptJetCorr < minJetPt)
        continue;

      registryQC.fill(HIST("hJetArea"), jetArea[i]);
    }
  }
  PROCESS_SWITCH(NucleiInJets, processGhosts, "Process Ghosts", false);

  void processDetResponseMatrix(SimCollisions const& collisions, MCTracks const& mcTracks, aod::McCollisions const&, const aod::McParticles&)
  {
    for (const auto& collision : collisions) { // o2-linter: disable=[const-ref-in-for-loop]

      // Event Selection
      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        continue;

      // List of Tracks and Particles
      std::vector<TVector3> trk;
      std::vector<TVector3> part;
      std::vector<int> ntrk;
      auto tracksPerColl = mcTracks.sliceBy(perCollision, collision.globalIndex());

      for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        if (!track.has_mcParticle())
          continue;
        const auto particle = track.mcParticle();

        TVector3 recMomentum(track.px(), track.py(), track.pz());
        TVector3 genMomentum(particle.px(), particle.py(), particle.pz());
        trk.push_back(recMomentum);
        part.push_back(genMomentum);
        ntrk.push_back(1);
      }
      // Reject Empty Events
      if (static_cast<int>(trk.size()) < 1)
        continue;

      // Anti-kt Jet Finder
      int nParticlesRemoved(0);
      std::vector<TVector3> jetRecMomentum;
      std::vector<TVector3> jetGenMomentum;
      std::vector<int> nParticlesInjet;

      do {
        double dijMin(1e+300), diBmin(1e+300);
        int iMin(0), jMin(0), iBmin(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diBmin) {
            diBmin = diB;
            iBmin = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
            if (trk[j].Mag() == 0)
              continue;
            double dij = calculateDij(trk[i], trk[j], rJet);
            if (dij < dijMin) {
              dijMin = dij;
              iMin = i;
              jMin = j;
            }
          }
        }
        if (dijMin < diBmin) {
          trk[iMin] = trk[iMin] + trk[jMin];
          ntrk[iMin] = ntrk[iMin] + ntrk[jMin];
          part[iMin] = part[iMin] + part[jMin];
          trk[jMin].SetXYZ(0, 0, 0);
          ntrk[jMin] = 0;
          nParticlesRemoved++;
        }
        if (dijMin > diBmin) {
          jetRecMomentum.push_back(trk[iBmin]);
          jetGenMomentum.push_back(part[iBmin]);
          nParticlesInjet.push_back(ntrk[iBmin]);
          trk[iBmin].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
      } while (nParticlesRemoved < static_cast<int>(trk.size()));

      for (int i = 0; i < static_cast<int>(jetRecMomentum.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        if ((std::fabs(jetRecMomentum[i].Eta()) + rJet) > (maxEta - 0.5))
          continue;

        double ptGen = jetGenMomentum[i].Pt();
        double ptRec = jetRecMomentum[i].Pt();
        registryMC.fill(HIST("detectorResponseMatrix"), ptGen, ptRec);
      }
    }
  }
  PROCESS_SWITCH(NucleiInJets, processDetResponseMatrix, "process detector response matrix", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<NucleiInJets>(cfgc)};
}
