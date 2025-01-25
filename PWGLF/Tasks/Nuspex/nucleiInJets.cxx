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
  Configurable<bool> requireNoOverlap{"requireNoOverlap", false, "require no overlap between jets and UE cones"};

  // Track Parameters
  Configurable<double> par0{"par0", 0.004, "par 0"};
  Configurable<double> par1{"par1", 0.013, "par 1"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNclusters{"minTpcNclusters", 80, "minimum number of TPC clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 80, "minimum number of TPC crossed pad rows"};
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

  // Bethe-bloch Parametrization of ITS cluster size
  Configurable<double> bbPar0{"bbPar0", 0.00089176700, "Bethe Bloch Par 0"};
  Configurable<double> bbPar1{"bbPar1", 33.9651487037, "Bethe Bloch Par 1"};
  Configurable<double> bbPar2{"bbPar2", 0.42595677370, "Bethe Bloch Par 2"};
  Configurable<double> bbPar3{"bbPar3", 1.39638691440, "Bethe Bloch Par 3"};
  Configurable<double> bbPar4{"bbPar4", 7.97312623880, "Bethe Bloch Par 4"};
  Configurable<double> bbPar5{"bbPar5", 61.3838254956, "Bethe Bloch Par 5"};
  Configurable<double> bbPar6{"bbPar6", 2.30000000000, "Bethe Bloch Par 6"};
  Configurable<double> bbPar7{"bbPar7", 0.93827208820, "Bethe Bloch Par 7"};
  Configurable<double> bbPar8{"bbPar8", 1.0, "Bethe Bloch Par 8"};
  Configurable<double> resolClsSize{"resolClsSize", 0.214, "Resolution of cls size distribution"};
  Configurable<double> nSigmaClsSizeMax{"nSigmaClsSizeMax", 2.0, "nSigma cut on cluster size"};
  Configurable<std::string> urlToCcdb{"urlToCcdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> pathToFile{"pathToFile", "", "path to file with reweighting"};
  Configurable<std::string> histoNameWeightAntipJet{"histoNameWeightAntipJet", "", "reweighting histogram: antip in jet"};
  Configurable<std::string> histoNameWeightAntipUe{"histoNameWeightAntipUe", "", "reweighting histogram: antip in ue"};

  // Bethe-Bloch
  TF1* bbClsSize = nullptr;

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
    registryQC.add("nJets_found", "nJets_found", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("nJets_selected", "nJets_selected", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("event_selection_jets", "event_selection_jets", HistType::kTH1F, {{10, 0, 10, "counter"}});
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

    bbClsSize = new TF1("bbClsSize", betheBloch, 0.1, 10, 9);
    bbClsSize->SetParameter(0, bbPar0);
    bbClsSize->SetParameter(1, bbPar1);
    bbClsSize->SetParameter(2, bbPar2);
    bbClsSize->SetParameter(3, bbPar3);
    bbClsSize->SetParameter(4, bbPar4);
    bbClsSize->SetParameter(5, bbPar5);
    bbClsSize->SetParameter(6, bbPar6);
    bbClsSize->SetParameter(7, bbPar7);
    bbClsSize->SetParameter(8, bbPar8);
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
    if (track.tpcNClsFound() < minTpcNclusters)
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows)
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

  double trackInclination(double eta)
  {
    double lambda(0);
    double theta = 2.0 * std::atan(std::exp(-eta));
    if (theta <= o2::constants::math::PIHalf)
      lambda = o2::constants::math::PIHalf - theta;
    if (theta > o2::constants::math::PIHalf)
      lambda = theta - o2::constants::math::PIHalf;
    return lambda;
  }

  static double betheBloch(double* x, double* par)
  {
    // 5 parameters for the bethe bloch from 0 to 4
    // 1 parameter for the mip mpar[5]
    // 1 parameter for the charge exponent mpar[6]
    // 1 parameter for the mass mpar[7]
    // 1 parameter for the charge mpar[8]
    return par[5] * betheBlochAleph(x[0] / par[7], par[0], par[1], par[2], par[3], par[4]) * std::pow(par[8], par[6]);
  }

  static double betheBlochAleph(double bg, double kp1, double kp2, double kp3, double kp4, double kp5)
  {
    double beta = bg / std::sqrt(1.0 + bg * bg);
    double aa = std::pow(beta, kp4);
    double bb = std::pow(1.0 / bg, kp5);
    bb = std::log(kp3 + bb);
    return (kp2 - aa - bb) * kp1 / aa;
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
    registryQC.fill(HIST("event_selection_jets"), 0.5); // all events before jet selection

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

    // List of Tracks
    std::vector<TVector3> trk;

    for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]

      if (!passedTrackSelectionForJetReconstruction(track))
        continue;
      registryQC.fill(HIST("dcaxy_vs_pt"), track.pt(), track.dcaXY());
      registryQC.fill(HIST("dcaz_vs_pt"), track.pt(), track.dcaZ());

      TVector3 momentum(track.px(), track.py(), track.pz());
      trk.push_back(momentum);
    }

    // Anti-kt Jet Finder
    int nParticlesRemoved(0);
    std::vector<TVector3> jet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    do {
      double dijMin(1e+06), diBmin(1e+06);
      int iMin(0), jMin(0), iB_min(0);
      for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
        if (trk[i].Mag() == 0)
          continue;
        double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
        if (diB < diBmin) {
          diBmin = diB;
          iB_min = i;
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
        trk[jMin].SetXYZ(0, 0, 0);
        nParticlesRemoved++;
      }
      if (dijMin > diBmin) {
        jet.push_back(trk[iB_min]);
        trk[iB_min].SetXYZ(0, 0, 0);
        nParticlesRemoved++;
      }
    } while (nParticlesRemoved < static_cast<int>(trk.size()));

    registryQC.fill(HIST("nJets_found"), static_cast<int>(jet.size()));

    // Jet Selection
    std::vector<int> isSelected;
    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
      isSelected.push_back(0);
    }

    int nJetsSelected(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

      if ((std::fabs(jet[i].Eta()) + rJet) > maxEta)
        continue;

      // Perpendicular cones
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jet[i], ueAxis1, +1);
      getPerpendicularAxis(jet[i], ueAxis2, -1);
      ue1.push_back(ueAxis1);
      ue2.push_back(ueAxis2);

      double nPartJetPlusUE(0);
      double nPartJet(0);
      double nPartUE(0);
      double ptJetPlusUE(0);
      double ptJet(0);
      double ptUE(0);

      for (auto track : tracks) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        TVector3 selectedTrack(track.px(), track.py(), track.pz());

        double deltaEtaJet = selectedTrack.Eta() - jet[i].Eta();
        double deltaPhiJet = getDeltaPhi(selectedTrack.Phi(), jet[i].Phi());
        double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
        double deltaEtaUe1 = selectedTrack.Eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(selectedTrack.Phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = selectedTrack.Eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(selectedTrack.Phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        if (deltaRjet < rJet) {
          registryQC.fill(HIST("deltaEtadeltaPhiJet"), deltaEtaJet, deltaPhiJet);
          nPartJetPlusUE++;
          ptJetPlusUE = ptJetPlusUE + selectedTrack.Pt();
        }
        if (deltaRUe1 < rJet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEtaUe1, deltaPhiUe1);
          nPartUE++;
          ptUE = ptUE + selectedTrack.Pt();
        }
        if (deltaRUe2 < rJet) {
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEtaUe2, deltaPhiUe2);
          nPartUE++;
          ptUE = ptUE + selectedTrack.Pt();
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

      if (ptJet < minJetPt)
        continue;
      if (nPartJetPlusUE < minNparticlesInJet)
        continue;
      nJetsSelected++;
      isSelected[i] = 1;
    }
    registryQC.fill(HIST("nJets_selected"), nJetsSelected);

    if (nJetsSelected == 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);
    registryQC.fill(HIST("event_selection_jets"), 1.5); // events with pTjet>10 GeV/c selected
    //************************************************************************************************************************************

    // Overlaps
    int nOverlaps(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
      if (isSelected[i] == 0)
        continue;

      for (int j = 0; j < static_cast<int>(jet.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
        if (isSelected[j] == 0 || i == j)
          continue;
        if (overlap(jet[i], ue1[j], rJet) || overlap(jet[i], ue2[j], rJet))
          nOverlaps++;
      }
    }
    registryQC.fill(HIST("jet_ue_overlaps"), nJetsSelected, nOverlaps);

    if (nJetsSelected > nJetsPerEventMax)
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    if (requireNoOverlap && nOverlaps > 0)
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    //************************************************************************************************************************************

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

        // ITS Cluster size
        double averageItsClusterSize(0);
        int nItsCls(0);
        for (int i = 0; i < 7; i++) { // o2-linter: disable=[const-ref-in-for-loop]
          int clusterSize = track.itsClsSizeInLayer(i);
          averageItsClusterSize += static_cast<double>(clusterSize);
          if (clusterSize > 0)
            nItsCls++;
        }
        averageItsClusterSize = averageItsClusterSize / static_cast<double>(nItsCls);
        double lambda = trackInclination(track.eta());
        double avgClsCosL = averageItsClusterSize * std::cos(lambda);
        double nsigma = (avgClsCosL - bbClsSize->Eval(pt)) / (resolClsSize * bbClsSize->Eval(pt));

        bool isItsSelected = false;
        if (std::fabs(nsigma) < nSigmaClsSizeMax) {
          isItsSelected = true;
        }

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
            if (isItsSelected) {
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
            if (isItsSelected) {
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

      for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        TVector3 momentum(track.px(), track.py(), track.pz());
        trk.push_back(momentum);
      }

      // Anti-kt Jet Finder
      int nParticlesRemoved(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      do {
        double dijMin(1e+06), diBmin(1e+06);
        int iMin(0), jMin(0), iB_min(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diBmin) {
            diBmin = diB;
            iB_min = i;
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
          trk[jMin].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
        if (dijMin > diBmin) {
          jet.push_back(trk[iB_min]);
          trk[iB_min].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
      } while (nParticlesRemoved < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
        isSelected.push_back(0);
      }

      int nJetsSelected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        if ((std::fabs(jet[i].Eta()) + rJet) > maxEta)
          continue;

        // Perpendicular cones
        TVector3 ueAxis1(0, 0, 0);
        TVector3 ueAxis2(0, 0, 0);
        getPerpendicularAxis(jet[i], ueAxis1, +1);
        getPerpendicularAxis(jet[i], ueAxis2, -1);
        ue1.push_back(ueAxis1);
        ue2.push_back(ueAxis2);

        double nPartJetPlusUE(0);
        double ptJetPlusUE(0);
        double ptJet(0);
        double ptUE(0);

        for (auto track : tracksPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

          if (!passedTrackSelectionForJetReconstruction(track))
            continue;
          TVector3 selectedTrack(track.px(), track.py(), track.pz());

          double deltaEtaJet = selectedTrack.Eta() - jet[i].Eta();
          double deltaPhiJet = getDeltaPhi(selectedTrack.Phi(), jet[i].Phi());
          double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = selectedTrack.Eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(selectedTrack.Phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = selectedTrack.Eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(selectedTrack.Phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          if (deltaRjet < rJet) {
            nPartJetPlusUE++;
            ptJetPlusUE = ptJetPlusUE + selectedTrack.Pt();
          }
          if (deltaRUe1 < rJet) {
            ptUE = ptUE + selectedTrack.Pt();
          }
          if (deltaRUe2 < rJet) {
            ptUE = ptUE + selectedTrack.Pt();
          }
        }
        ptJet = ptJetPlusUE - 0.5 * ptUE;

        if (ptJet < minJetPt)
          continue;
        if (nPartJetPlusUE < minNparticlesInJet)
          continue;
        nJetsSelected++;
        isSelected[i] = 1;
      }
      if (nJetsSelected == 0)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 6.5);

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

        if (setDCAselectionPtDep) {
          if (dcaxy > (par0 + par1 / particle.pt()))
            continue;
          if (dcaz > (par0 + par1 / particle.pt()))
            continue;
        }
        if (!setDCAselectionPtDep) {
          if (dcaxy > maxDcaxy)
            continue;
          if (dcaz > maxDcaz)
            continue;
        }

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
      }

      // Anti-kt Jet Finder
      int nParticlesRemoved(0);
      std::vector<TVector3> jet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      do {
        double dijMin(1e+06), diBmin(1e+06);
        int iMin(0), jMin(0), iB_min(0);
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diBmin) {
            diBmin = diB;
            iB_min = i;
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
          trk[jMin].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
        if (dijMin > diBmin) {
          jet.push_back(trk[iB_min]);
          trk[iB_min].SetXYZ(0, 0, 0);
          nParticlesRemoved++;
        }
      } while (nParticlesRemoved < static_cast<int>(trk.size()));

      // Jet Selection
      std::vector<int> isSelected;
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
        isSelected.push_back(0);
      }

      int nJetsSelected(0);
      for (int i = 0; i < static_cast<int>(jet.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]

        if ((std::fabs(jet[i].Eta()) + rJet) > maxEta)
          continue;

        // Perpendicular cones
        TVector3 ueAxis1(0, 0, 0);
        TVector3 ueAxis2(0, 0, 0);
        getPerpendicularAxis(jet[i], ueAxis1, +1);
        getPerpendicularAxis(jet[i], ueAxis2, -1);
        ue1.push_back(ueAxis1);
        ue2.push_back(ueAxis2);

        double nPartJetPlusUE(0);
        double ptJetPlusUE(0);
        double ptJet(0);
        double ptUE(0);

        for (auto& particle : mcParticlesPerColl) { // o2-linter: disable=[const-ref-in-for-loop]

          // Select Primary Particles
          double dx = particle.vx() - mccollision.posX();
          double dy = particle.vy() - mccollision.posY();
          double dz = particle.vz() - mccollision.posZ();
          double dcaxy = std::sqrt(dx * dx + dy * dy);
          double dcaz = std::fabs(dz);

          if (setDCAselectionPtDep) {
            if (dcaxy > (par0 + par1 / particle.pt()))
              continue;
            if (dcaz > (par0 + par1 / particle.pt()))
              continue;
          }
          if (!setDCAselectionPtDep) {
            if (dcaxy > maxDcaxy)
              continue;
            if (dcaz > maxDcaz)
              continue;
          }

          if (std::fabs(particle.eta()) > 0.8)
            continue;
          if (particle.pt() < 0.15)
            continue;

          // PDG Selection
          int pdg = std::fabs(particle.pdgCode());
          if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
            continue;

          TVector3 selectedTrack(particle.px(), particle.py(), particle.pz());

          double deltaEtaJet = selectedTrack.Eta() - jet[i].Eta();
          double deltaPhiJet = getDeltaPhi(selectedTrack.Phi(), jet[i].Phi());
          double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = selectedTrack.Eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(selectedTrack.Phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = selectedTrack.Eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(selectedTrack.Phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          if (deltaRjet < rJet) {
            nPartJetPlusUE++;
            ptJetPlusUE = ptJetPlusUE + selectedTrack.Pt();
          }
          if (deltaRUe1 < rJet) {
            ptUE = ptUE + selectedTrack.Pt();
          }
          if (deltaRUe2 < rJet) {
            ptUE = ptUE + selectedTrack.Pt();
          }
        }
        ptJet = ptJetPlusUE - 0.5 * ptUE;

        if (ptJet < minJetPt)
          continue;
        if (nPartJetPlusUE < minNparticlesInJet)
          continue;
        nJetsSelected++;
        isSelected[i] = 1;
      }
      if (nJetsSelected == 0)
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<NucleiInJets>(cfgc)};
}
