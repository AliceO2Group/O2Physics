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
/// \file strangeness in jets.cxx
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since May 22, 2024

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>
#include <cmath>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using StrHadronDaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MCTracks = soa::Join<StrHadronDaughterTracks, aod::McTrackLabels>;

struct strangeness_in_jets {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global Parameters
  Configurable<int> particle_of_interest{"particle_of_interest", 0, "0=v0, 1=cascade, 2=pions"};
  Configurable<double> min_jet_pt{"min_jet_pt", 10.0, "Minimum pt of the jet"};
  Configurable<double> Rjet{"Rjet", 0.3, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<int> min_nPartInJet{"min_nPartInJet", 2, "Minimum number of particles inside jet"};
  Configurable<int> n_jets_per_event_max{"n_jets_per_event_max", 1000, "Maximum number of jets per event"};
  Configurable<bool> requireNoOverlap{"requireNoOverlap", false, "require no overlap between jets and UE cones"};
  Configurable<float> par0{"par0", 0.004f, "par 0"};
  Configurable<float> par1{"par1", 0.013f, "par 1"};
  // Axis parameters
  struct : ConfigurableGroup {
    ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};

    ConfigurableAxis binsDca{"binsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};

  } binsOptions;

  // Track Parameters
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> ptMin_V0_proton{"ptMin_V0_proton", 0.3f, "pt min of proton from V0"};
  Configurable<float> ptMax_V0_proton{"ptMax_V0_proton", 10.0f, "pt max of proton from V0"};
  Configurable<float> ptMin_V0_pion{"ptMin_V0_pion", 0.1f, "pt min of pion from V0"};
  Configurable<float> ptMax_V0_pion{"ptMax_V0_pion", 1.5f, "pt max of pion from V0"};
  Configurable<float> ptMin_K0_pion{"ptMin_K0_pion", 0.3f, "pt min of pion from K0"};
  Configurable<float> ptMax_K0_pion{"ptMax_K0_pion", 10.0f, "pt max of pion from K0"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<float> nsigmaTOFminPionMC{"nsigmaTOFminPionMC", -6.0f, "Minimum nsigma TOF for pion in MC"};
  Configurable<float> nsigmaTOFmaxPionMC{"nsigmaTOFmaxPionMC", +6.0f, "Maximum nsigma TOF for pion in MC"};
  Configurable<float> dcaxyMax{"dcaxyMax", 0.1f, "Maximum DCAxy to primary vertex"};
  Configurable<float> dcazMax{"dcazMax", 0.1f, "Maximum DCAz to primary vertex"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};

  // V0 Parameters
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 40.0f, "Maximum V0 Radius"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Cascade Parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 0.1f, "Minimum Cascade Radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 40.0f, "Maximum Cascade Radius"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum Cascade CosPA"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.1f, "Minimum DCA bachelor to PV"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.1f, "Minimum DCA V0 to PV"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Re-weighting
  Configurable<bool> applyReweighting{"applyReweighting", true, "apply reweighting"};
  Configurable<std::string> url_to_ccdb{"url_to_ccdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> path_to_file{"path_to_file", "", "path to file with reweighting"};
  Configurable<std::string> histo_name_weight_piplus_jet{"histo_name_weight_piplus_jet", "", "reweighting histogram: Pi Plus in jet"};
  Configurable<std::string> histo_name_weight_piplus_ue{"histo_name_weight_piplus_ue", "", "reweighting histogram: Pi Plus in ue"};
  Configurable<std::string> histo_name_weight_piminus_jet{"histo_name_weight_piminus_jet", "", "reweighting histogram: Pi Minus in jet"};
  Configurable<std::string> histo_name_weight_piminus_ue{"histo_name_weight_piminus_ue", "", "reweighting histogram: Pi Minus in ue"};
  Configurable<std::string> histo_name_weight_k0_jet{"histo_name_weight_k0_jet", "", "reweighting histogram: K0 in jet"};
  Configurable<std::string> histo_name_weight_k0_ue{"histo_name_weight_k0_ue", "", "reweighting histogram: K0 in ue"};
  Configurable<std::string> histo_name_weight_lambda_jet{"histo_name_weight_lambda_jet", "", "reweighting histogram: lambda in jet"};
  Configurable<std::string> histo_name_weight_lambda_ue{"histo_name_weight_lambda_ue", "", "reweighting histogram: lambda in ue"};
  Configurable<std::string> histo_name_weight_antilambda_jet{"histo_name_weight_antilambda_jet", "", "reweighting histogram: antilambda in jet"};
  Configurable<std::string> histo_name_weight_antilambda_ue{"histo_name_weight_antilambda_ue", "", "reweighting histogram: antilambda in ue"};

  // Two-dimensional weights
  TH2F* twod_weights_piplus_jet = nullptr;
  TH2F* twod_weights_piplus_ue = nullptr;
  TH2F* twod_weights_piminus_jet = nullptr;
  TH2F* twod_weights_piminus_ue = nullptr;
  TH2F* twod_weights_k0_jet;
  TH2F* twod_weights_k0_ue;
  TH2F* twod_weights_lambda_jet;
  TH2F* twod_weights_lambda_ue;
  TH2F* twod_weights_antilambda_jet;
  TH2F* twod_weights_antilambda_ue;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // List of Particles
  enum option { vzeros,
                cascades,
                pions };

  void init(InitContext const&)
  {
    ccdb->setURL(url_to_ccdb.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (applyReweighting) {
      GetReweightingHistograms(ccdb, TString(path_to_file), TString(histo_name_weight_k0_jet), TString(histo_name_weight_k0_ue), TString(histo_name_weight_lambda_jet), TString(histo_name_weight_lambda_ue), TString(histo_name_weight_antilambda_jet), TString(histo_name_weight_antilambda_ue));
    } else {
      twod_weights_k0_jet = nullptr;
      twod_weights_k0_ue = nullptr;
      twod_weights_lambda_jet = nullptr;
      twod_weights_lambda_ue = nullptr;
      twod_weights_antilambda_jet = nullptr;
      twod_weights_antilambda_ue = nullptr;
    }

    // Event Counters
    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
    registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {{101, 0, 101, "Multiplicity percentile"}});
    registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});

    // QC Histograms
    registryQC.add("deltaEtadeltaPhi_jet", "deltaEtadeltaPhi_jet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, o2::constants::math::PIHalf, "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhi_ue", "deltaEtadeltaPhi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, o2::constants::math::PIHalf, "#Delta#phi"}});
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
    registryQC.add("survivedK0", "survivedK0", HistType::kTH1F, {{20, 0, 20, "step"}});

    // Multiplicity Binning
    std::vector<double> multBinning = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec multAxis = {multBinning, "FT0C percentile"};
    const AxisSpec etaAxis{18, -0.9, 0.9, "#eta"};
    const AxisSpec ptAxis{binsOptions.binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nsigmaTOFAxis{200, -10, 10, "n#sigma_{TOF}"};
    const AxisSpec nsigmaTPCAxis{200, -10, 10, "n#sigma_{TPC}"};
    const AxisSpec dcaAxis{binsOptions.binsDca, "DCA_{xy} (cm)"};

    // Histograms for pions (data)
    registryData.add("piplus_tpc_in_jet", "piplus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxis, nsigmaTPCAxis});
    registryData.add("piplus_tof_in_jet", "piplus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxis, nsigmaTOFAxis});
    registryData.add("piplus_tpc_in_ue", "piplus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxis, nsigmaTPCAxis});
    registryData.add("piplus_tof_in_ue", "piplus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxis, nsigmaTOFAxis});
    registryData.add("piminus_tpc_in_jet", "piminus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxis, nsigmaTPCAxis});
    registryData.add("piminus_tof_in_jet", "piminus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxis, nsigmaTOFAxis});
    registryData.add("piminus_tpc_in_ue", "piminus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxis, nsigmaTPCAxis});
    registryData.add("piminus_tof_in_ue", "piminus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxis, nsigmaTOFAxis});

    registryData.add("piplus_dcaxy_in_jet", "piplus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryData.add("piplus_dcaxy_in_ue", "piplus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryData.add("piminus_dcaxy_in_jet", "piminus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryData.add("piminus_dcaxy_in_ue", "piminus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});

    // Histograms for lambda (data)
    registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});

    // Histograms for K0s (data)
    registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"}});

    // Histograms for xi (data)
    registryData.add("XiPos_in_jet", "XiPos_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiPos_in_ue", "XiPos_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiNeg_in_jet", "XiNeg_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiNeg_in_ue", "XiNeg_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});

    // Histograms for omega (data)
    registryData.add("OmegaPos_in_jet", "OmegaPos_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaPos_in_ue", "OmegaPos_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaNeg_in_jet", "OmegaNeg_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaNeg_in_ue", "OmegaNeg_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});

    // Histograms for efficiency (generated)
    registryMC.add("K0s_generated_jet", "K0s_generated_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("K0s_generated_ue", "K0s_generated_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_generated_jet", "Lambda_generated_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_generated_ue", "Lambda_generated_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_generated_jet", "AntiLambda_generated_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_generated_ue", "AntiLambda_generated_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiPos_generated", "XiPos_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiNeg_generated", "XiNeg_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaPos_generated", "OmegaPos_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaNeg_generated", "OmegaNeg_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Histograms for efficiency (reconstructed)
    registryMC.add("K0s_reconstructed_jet", "K0s_reconstructed_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("K0s_reconstructed_ue", "K0s_reconstructed_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_reconstructed_jet", "Lambda_reconstructed_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_reconstructed_ue", "Lambda_reconstructed_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_reconstructed_jet", "AntiLambda_reconstructed_jet", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_reconstructed_ue", "AntiLambda_reconstructed_ue", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiPos_reconstructed", "XiPos_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiNeg_reconstructed", "XiNeg_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaPos_reconstructed", "OmegaPos_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaNeg_reconstructed", "OmegaNeg_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Histograms for secondary hadrons
    registryMC.add("K0s_reconstructed_incl", "K0s_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_reconstructed_incl", "Lambda_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_reconstructed_incl", "AntiLambda_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Histograms for 2d reweighting (pion)
    registryMC.add("Pion_eta_pt_jet", "Pion_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
    registryMC.add("Pion_eta_pt_ue", "Pion_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
    registryMC.add("Pion_eta_pt_pythia", "Pion_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});

    // Histograms for 2d reweighting (K0s)
    registryMC.add("K0s_eta_pt_jet", "K0s_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("K0s_eta_pt_ue", "K0s_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("K0s_eta_pt_pythia", "K0s_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});

    // Histograms for 2d reweighting (Lambda)
    registryMC.add("Lambda_eta_pt_jet", "Lambda_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Lambda_eta_pt_ue", "Lambda_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Lambda_eta_pt_pythia", "Lambda_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiLambda_eta_pt_jet", "AntiLambda_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiLambda_eta_pt_ue", "AntiLambda_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiLambda_eta_pt_pythia", "AntiLambda_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});

    // Histograms for 2d reweighting (Xi)
    registryMC.add("Xi_eta_pt_jet", "Xi_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Xi_eta_pt_ue", "Xi_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Xi_eta_pt_pythia", "Xi_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiXi_eta_pt_jet", "AntiXi_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiXi_eta_pt_ue", "AntiXi_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiXi_eta_pt_pythia", "AntiXi_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});

    // Histograms for 2d reweighting (Omega)
    registryMC.add("Omega_eta_pt_jet", "Omega_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Omega_eta_pt_ue", "Omega_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("Omega_eta_pt_pythia", "Omega_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiOmega_eta_pt_jet", "AntiOmega_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiOmega_eta_pt_ue", "AntiOmega_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});
    registryMC.add("AntiOmega_eta_pt_pythia", "AntiOmega_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, etaAxis});

    // Histograms for efficiency (pions)
    registryMC.add("pi_plus_gen_in_jet", "pi_plus_gen_in_jet", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_plus_rec_in_jet_tpc", "pi_plus_rec_in_jet_tpc", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_plus_rec_in_jet_tof", "pi_plus_rec_in_jet_tof", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_plus_gen_in_ue", "pi_plus_gen_in_ue", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_plus_rec_in_ue_tpc", "pi_plus_rec_in_ue_tpc", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_plus_rec_in_ue_tof", "pi_plus_rec_in_ue_tof", HistType::kTH2F, {multBinning, ptAxis});

    registryMC.add("pi_minus_gen_in_jet", "pi_minus_gen_in_jet", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_minus_rec_in_jet_tpc", "pi_minus_rec_in_jet_tpc", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_minus_rec_in_jet_tof", "pi_minus_rec_in_jet_tof", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_minus_gen_in_ue", "pi_minus_gen_in_ue", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_minus_rec_in_ue_tpc", "pi_minus_rec_in_ue_tpc", HistType::kTH2F, {multBinning, ptAxis});
    registryMC.add("pi_minus_rec_in_ue_tof", "pi_minus_rec_in_ue_tof", HistType::kTH2F, {multBinning, ptAxis});

    // MC Templates
    registryMC.add("piplus_dcaxy_prim", "piplus_dcaxy_prim", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryMC.add("piminus_dcaxy_prim", "piminus_dcaxy_prim", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryMC.add("piplus_dcaxy_sec", "piplus_dcaxy_sec", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
    registryMC.add("piminus_dcaxy_sec", "piminus_dcaxy_sec", HistType::kTH3F, {multBinning, ptAxis, dcaAxis});
  }

  template <typename chargedTrack>
  bool passedTrackSelectionForJetReconstruction(const chargedTrack& track)
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
    if (std::fabs(track.dcaXY()) > (par0 + par1 / track.pt()))
      return false;
    if (std::fabs(track.dcaZ()) > (par0 + par1 / track.pt()))
      return false;
    return true;
  }

  template <typename pionTrack>
  bool passedTrackSelectionForPions(const pionTrack& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (std::fabs(track.dcaZ()) > dcazMax)
      return false;
    return true;
  }

  // Lambda Selections
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum of Lambda Daughters
    TVector3 proton(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (proton.Pt() < ptMin_V0_proton)
      return false;
    if (proton.Pt() > ptMax_V0_proton)
      return false;
    if (pion.Pt() < ptMin_V0_pion)
      return false;
    if (pion.Pt() > ptMax_V0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // AntiLambda Selections
  template <typename AntiLambda, typename TrackPos, typename TrackNeg>
  bool passedAntiLambdaSelection(const AntiLambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum AntiLambda Daughters
    TVector3 pion(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 proton(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (proton.Pt() < ptMin_V0_proton)
      return false;
    if (proton.Pt() > ptMax_V0_proton)
      return false;
    if (pion.Pt() < ptMin_V0_pion)
      return false;
    if (pion.Pt() > ptMax_V0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // K0s Selections
  template <typename K0short, typename TrackPos, typename TrackNeg>
  bool passedK0ShortSelection(const K0short& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum K0s Daughters
    TVector3 pion_pos(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion_neg(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (pion_pos.Pt() < ptMin_K0_pion)
      return false;
    if (pion_pos.Pt() > ptMax_K0_pion)
      return false;
    if (pion_neg.Pt() < ptMin_K0_pion)
      return false;
    if (pion_neg.Pt() > ptMax_K0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Xi Selections
  template <typename Xi, typename TrackPos, typename TrackNeg, typename TrackBac, typename Coll>
  bool passedXiSelection(const Xi& casc, const TrackPos& ptrack, const TrackNeg& ntrack, const TrackBac& btrack, const Coll& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Xi+ Selection (Xi+ -> antiL + pi+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Xi- Selection (Xi- -> L + pi-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin || btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaPi() < nsigmaTOFmin || btrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Omega Selections
  template <typename Omega, typename TrackPos, typename TrackNeg, typename TrackBac, typename Coll>
  bool passedOmegaSelection(const Omega& casc, const TrackPos& ptrack, const TrackNeg& ntrack, const TrackBac& btrack, const Coll& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Omega+ Selection (Omega+ -> antiL + K+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Omega- Selection (Omega- -> L + K-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin || btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaKa() < nsigmaTOFmin || btrack.tofNSigmaKa() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Single-Track Selection
  template <typename Track>
  bool passedSingleTrackSelection(const Track& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (requireTOF && (!track.hasTOF()))
      return false;
    return true;
  }

  // Pion Selection
  template <typename pionTrack>
  bool isHighPurityPion(const pionTrack& track)
  {
    if (track.p() < 0.6 && std::fabs(track.tpcNSigmaPi()) < 3.0)
      return true;
    if (track.p() > 0.6 && std::fabs(track.tpcNSigmaPi()) < 3.0 && std::fabs(track.tofNSigmaPi()) < 3.0)
      return true;
    return false;
  }

  float Minimum(float x1, float x2)
  {
    float x_min(x1);
    if (x1 < x2)
      x_min = x1;
    if (x1 >= x2)
      x_min = x2;
    return x_min;
  }

  double GetDeltaPhi(double a1, double a2)
  {
    double delta_phi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = std::fabs(phi1 - phi2);

    if (diff <= o2::constants::math::PI)
      delta_phi = diff;
    if (diff > o2::constants::math::PI)
      delta_phi = o2::constants::math::TwoPI - diff;

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

    // Protection agains delta<0
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
    double d = std::sqrt(dx * dx + dy * dy);
    if (d < 2.0 * R)
      return true;
    return false;
  }

  void GetReweightingHistograms(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, TString filepath, TString histname_k0_jet, TString histname_k0_ue, TString histname_lambda_jet, TString histname_lambda_ue, TString histname_antilambda_jet, TString histname_antilambda_ue)
  {
    TList* l = ccdbObj->get<TList>(filepath.Data());
    if (!l) {
      LOGP(error, "Could not open the file {}", Form("%s", filepath.Data()));
      return;
    }

    if (histo_name_weight_piplus_jet.value != "") {
      twod_weights_piplus_jet = ccdbObj->get<TH2F>(filepath.Data());
      LOG(info) << "Getting weight histogram for piplus in jet from " << histo_name_weight_piplus_jet.value;
    }
    if (histo_name_weight_piplus_ue.value != "") {
      twod_weights_piplus_ue = ccdbObj->get<TH2F>(filepath.Data());
      LOG(info) << "Getting weight histogram for piplus in ue from " << histo_name_weight_piplus_ue.value;
    }
    if (histo_name_weight_piminus_jet.value != "") {
      twod_weights_piminus_jet = ccdbObj->get<TH2F>(filepath.Data());
      LOG(info) << "Getting weight histogram for piminus in jet from " << histo_name_weight_piminus_jet.value;
    }
    if (histo_name_weight_piminus_ue.value != "") {
      twod_weights_piminus_ue = ccdbObj->get<TH2F>(filepath.Data());
      LOG(info) << "Getting weight histogram for piminus in ue from " << histo_name_weight_piminus_ue.value;
    }

    twod_weights_k0_jet = static_cast<TH2F*>(l->FindObject(Form("%s", histname_k0_jet.Data())));
    if (!twod_weights_k0_jet) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_k0_jet.Data()));
      return;
    }
    twod_weights_k0_ue = static_cast<TH2F*>(l->FindObject(Form("%s", histname_k0_ue.Data())));
    if (!twod_weights_k0_ue) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_k0_ue.Data()));
      return;
    }
    twod_weights_lambda_jet = static_cast<TH2F*>(l->FindObject(Form("%s", histname_lambda_jet.Data())));
    if (!twod_weights_lambda_jet) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_lambda_jet.Data()));
      return;
    }
    twod_weights_lambda_ue = static_cast<TH2F*>(l->FindObject(Form("%s", histname_lambda_ue.Data())));
    if (!twod_weights_lambda_ue) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_lambda_ue.Data()));
      return;
    }
    twod_weights_antilambda_jet = static_cast<TH2F*>(l->FindObject(Form("%s", histname_antilambda_jet.Data())));
    if (!twod_weights_antilambda_jet) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_antilambda_jet.Data()));
      return;
    }
    twod_weights_antilambda_ue = static_cast<TH2F*>(l->FindObject(Form("%s", histname_antilambda_ue.Data())));
    if (!twod_weights_antilambda_ue) {
      LOGP(error, "Could not open histogram {}", Form("%s", histname_antilambda_ue.Data()));
      return;
    }

    LOGP(info, "Opened histogram {}", Form("%s", histname_k0_jet.Data()));
    LOGP(info, "Opened histogram {}", Form("%s", histname_k0_ue.Data()));
    LOGP(info, "Opened histogram {}", Form("%s", histname_lambda_jet.Data()));
    LOGP(info, "Opened histogram {}", Form("%s", histname_lambda_ue.Data()));
    LOGP(info, "Opened histogram {}", Form("%s", histname_antilambda_jet.Data()));
    LOGP(info, "Opened histogram {}", Form("%s", histname_antilambda_ue.Data()));
  }

  void processData(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s, aod::CascDataExt const& Cascades, StrHadronDaughterTracks const& tracks)
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

    // List of Tracks
    std::vector<TVector3> trk;

    for (const auto& track : tracks) {

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
    registryQC.fill(HIST("nJets_found"), static_cast<int>(jet.size()));

    // Jet Selection
    std::vector<int> isSelected;
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {
      isSelected.push_back(0);
    }

    int n_jets_selected(0);
    for (int i = 0; i < static_cast<int>(jet.size()); i++) {

      if ((std::fabs(jet[i].Eta()) + Rjet) > etaMax)
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

      for (const auto& track : tracks) {

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;
        TVector3 sel_track(track.px(), track.py(), track.pz());

        double deltaEta_jet = sel_track.Eta() - jet[i].Eta();
        double deltaPhi_jet = GetDeltaPhi(sel_track.Phi(), jet[i].Phi());
        double deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        double deltaEta_ue1 = sel_track.Eta() - ue_axis1.Eta();
        double deltaPhi_ue1 = GetDeltaPhi(sel_track.Phi(), ue_axis1.Phi());
        double deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        double deltaEta_ue2 = sel_track.Eta() - ue_axis2.Eta();
        double deltaPhi_ue2 = GetDeltaPhi(sel_track.Phi(), ue_axis2.Phi());
        double deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

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
          registryQC.fill(HIST("deltaEtadeltaPhi_ue"), deltaEta_ue1, deltaPhi_ue1);
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

    // Event multiplicity
    float multiplicity = collision.centFT0M();

    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);

    for (int i = 0; i < static_cast<int>(jet.size()); i++) {

      if (isSelected[i] == 0)
        continue;

      // Vzeros
      if (particle_of_interest == option::vzeros) {
        for (const auto& v0 : fullV0s) {

          const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
          const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();
          TVector3 v0dir(v0.px(), v0.py(), v0.pz());

          float deltaEta_jet = v0dir.Eta() - jet[i].Eta();
          float deltaPhi_jet = GetDeltaPhi(v0dir.Phi(), jet[i].Phi());
          float deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          float deltaEta_ue1 = v0dir.Eta() - ue1[i].Eta();
          float deltaPhi_ue1 = GetDeltaPhi(v0dir.Phi(), ue1[i].Phi());
          float deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          float deltaEta_ue2 = v0dir.Eta() - ue2[i].Eta();
          float deltaPhi_ue2 = GetDeltaPhi(v0dir.Phi(), ue2[i].Phi());
          float deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          // K0s
          if (passedK0ShortSelection(v0, pos, neg)) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("K0s_in_jet"), multiplicity, v0.pt(), v0.mK0Short());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("K0s_in_ue"), multiplicity, v0.pt(), v0.mK0Short());
            }
          }
          // Lambda
          if (passedLambdaSelection(v0, pos, neg)) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("Lambda_in_jet"), multiplicity, v0.pt(), v0.mLambda());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("Lambda_in_ue"), multiplicity, v0.pt(), v0.mLambda());
            }
          }
          // AntiLambda
          if (passedAntiLambdaSelection(v0, pos, neg)) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("AntiLambda_in_jet"), multiplicity, v0.pt(), v0.mAntiLambda());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("AntiLambda_in_ue"), multiplicity, v0.pt(), v0.mAntiLambda());
            }
          }
        }
      }

      // Cascades
      if (particle_of_interest == option::cascades) {
        for (const auto& casc : Cascades) {

          auto bach = casc.bachelor_as<StrHadronDaughterTracks>();
          auto pos = casc.posTrack_as<StrHadronDaughterTracks>();
          auto neg = casc.negTrack_as<StrHadronDaughterTracks>();

          TVector3 cascade_dir(casc.px(), casc.py(), casc.pz());
          float deltaEta_jet = cascade_dir.Eta() - jet[i].Eta();
          float deltaPhi_jet = GetDeltaPhi(cascade_dir.Phi(), jet[i].Phi());
          float deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          float deltaEta_ue1 = cascade_dir.Eta() - ue1[i].Eta();
          float deltaPhi_ue1 = GetDeltaPhi(cascade_dir.Phi(), ue1[i].Phi());
          float deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          float deltaEta_ue2 = cascade_dir.Eta() - ue2[i].Eta();
          float deltaPhi_ue2 = GetDeltaPhi(cascade_dir.Phi(), ue2[i].Phi());
          float deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          // Xi+
          if (passedXiSelection(casc, pos, neg, bach, collision) &&
              bach.sign() > 0) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("XiPos_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("XiPos_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Xi-
          if (passedXiSelection(casc, pos, neg, bach, collision) &&
              bach.sign() < 0) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("XiNeg_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("XiNeg_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Omega+
          if (passedOmegaSelection(casc, pos, neg, bach, collision) &&
              bach.sign() > 0) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("OmegaPos_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("OmegaPos_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
          // Omega-
          if (passedOmegaSelection(casc, pos, neg, bach, collision) &&
              bach.sign() < 0) {
            if (deltaR_jet < Rjet) {
              registryData.fill(HIST("OmegaNeg_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryData.fill(HIST("OmegaNeg_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
        }
      }

      // Pions
      if (particle_of_interest == option::pions) {
        for (const auto& track : tracks) {

          if (!passedTrackSelectionForPions(track))
            continue;

          TVector3 track_dir(track.px(), track.py(), track.pz());
          float deltaEta_jet = track_dir.Eta() - jet[i].Eta();
          float deltaPhi_jet = GetDeltaPhi(track_dir.Phi(), jet[i].Phi());
          float deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          float deltaEta_ue1 = track_dir.Eta() - ue1[i].Eta();
          float deltaPhi_ue1 = GetDeltaPhi(track_dir.Phi(), ue1[i].Phi());
          float deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          float deltaEta_ue2 = track_dir.Eta() - ue2[i].Eta();
          float deltaPhi_ue2 = GetDeltaPhi(track_dir.Phi(), ue2[i].Phi());
          float deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          bool isInJet = false;
          bool isInUe = false;
          if (deltaR_jet < Rjet)
            isInJet = true;
          if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet)
            isInUe = true;

          if (isHighPurityPion(track) && track.sign() > 0) {
            if (isInJet)
              registryData.fill(HIST("piplus_dcaxy_in_jet"), multiplicity, track.pt(), track.dcaXY());
            if (isInUe)
              registryData.fill(HIST("piplus_dcaxy_in_ue"), multiplicity, track.pt(), track.dcaXY());
          }
          if (isHighPurityPion(track) && track.sign() < 0) {
            if (isInJet)
              registryData.fill(HIST("piminus_dcaxy_in_jet"), multiplicity, track.pt(), track.dcaXY());
            if (isInUe)
              registryData.fill(HIST("piminus_dcaxy_in_ue"), multiplicity, track.pt(), track.dcaXY());
          }

          // DCAxy Selection
          if (std::fabs(track.dcaXY()) > dcaxyMax)
            continue;

          // TPC
          if (isInJet && track.sign() > 0) {
            registryData.fill(HIST("piplus_tpc_in_jet"), multiplicity, track.pt(), track.tpcNSigmaPi());
          }
          if (isInUe && track.sign() > 0) {
            registryData.fill(HIST("piplus_tpc_in_ue"), multiplicity, track.pt(), track.tpcNSigmaPi());
          }
          if (isInJet && track.sign() < 0) {
            registryData.fill(HIST("piminus_tpc_in_jet"), multiplicity, track.pt(), track.tpcNSigmaPi());
          }
          if (isInUe && track.sign() < 0) {
            registryData.fill(HIST("piminus_tpc_in_ue"), multiplicity, track.pt(), track.tpcNSigmaPi());
          }
          if (track.tpcNSigmaPi() < nsigmaTPCmin || track.tpcNSigmaPi() > nsigmaTPCmax)
            continue;
          if (!track.hasTOF())
            continue;

          // TOF
          if (isInJet && track.sign() > 0) {
            registryData.fill(HIST("piplus_tof_in_jet"), multiplicity, track.pt(), track.tofNSigmaPi());
          }
          if (isInUe && track.sign() > 0) {
            registryData.fill(HIST("piplus_tof_in_ue"), multiplicity, track.pt(), track.tofNSigmaPi());
          }
          if (isInJet && track.sign() < 0) {
            registryData.fill(HIST("piminus_tof_in_jet"), multiplicity, track.pt(), track.tofNSigmaPi());
          }
          if (isInUe && track.sign() < 0) {
            registryData.fill(HIST("piminus_tof_in_ue"), multiplicity, track.pt(), track.tofNSigmaPi());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processData, "Process data", true);

  void processK0s(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const&)
  {
    registryData.fill(HIST("number_of_events_data"), 10.5);
    if (!collision.sel8())
      return;
    registryData.fill(HIST("number_of_events_data"), 11.5);
    if (std::fabs(collision.posZ()) > zVtx)
      return;
    registryData.fill(HIST("number_of_events_data"), 12.5);

    for (const auto& v0 : fullV0s) {
      const auto& ptrack = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& ntrack = v0.negTrack_as<StrHadronDaughterTracks>();

      registryQC.fill(HIST("survivedK0"), 0.5);

      // Single-Track Selections
      if (!passedSingleTrackSelection(ptrack))
        continue;
      registryQC.fill(HIST("survivedK0"), 1.5);

      if (!passedSingleTrackSelection(ntrack))
        continue;
      registryQC.fill(HIST("survivedK0"), 2.5);

      // Momentum K0s Daughters
      TVector3 pion_pos(v0.pxpos(), v0.pypos(), v0.pzpos());
      TVector3 pion_neg(v0.pxneg(), v0.pyneg(), v0.pzneg());

      if (pion_pos.Pt() < ptMin_K0_pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 3.5);

      if (pion_pos.Pt() > ptMax_K0_pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 4.5);

      if (pion_neg.Pt() < ptMin_K0_pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 5.5);

      if (pion_neg.Pt() > ptMax_K0_pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 6.5);

      // V0 Selections
      if (v0.v0cosPA() < v0cospaMin)
        continue;
      registryQC.fill(HIST("survivedK0"), 7.5);

      if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
        continue;
      registryQC.fill(HIST("survivedK0"), 8.5);

      if (v0.dcaV0daughters() > dcaV0DaughtersMax)
        continue;
      registryQC.fill(HIST("survivedK0"), 9.5);

      if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
        continue;
      registryQC.fill(HIST("survivedK0"), 10.5);

      if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
        continue;
      registryQC.fill(HIST("survivedK0"), 11.5);

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        continue;
      registryQC.fill(HIST("survivedK0"), 12.5);

      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        continue;
      registryQC.fill(HIST("survivedK0"), 13.5);

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          continue;
        registryQC.fill(HIST("survivedK0"), 14.5);

        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          continue;
        registryQC.fill(HIST("survivedK0"), 15.5);
      }
    }

    for (const auto& v0 : fullV0s) {
      const auto& ptrack = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& ntrack = v0.negTrack_as<StrHadronDaughterTracks>();
      if (!passedK0ShortSelection(v0, ptrack, ntrack))
        continue;
      registryQC.fill(HIST("survivedK0"), 16.5);
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processK0s, "Process K0s", false);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDataExt> perCollisionCasc = o2::aod::cascade::collisionId;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<MCTracks> perCollisionTrk = o2::aod::track::collisionId;

  void processMCefficiency(SimCollisions const& collisions, MCTracks const& mcTracks, aod::V0Datas const& fullV0s, aod::CascDataExt const& Cascades, const aod::McParticles& mcParticles)
  {
    for (const auto& collision : collisions) {
      registryMC.fill(HIST("number_of_events_mc"), 0.5);
      if (!collision.sel8())
        continue;

      registryMC.fill(HIST("number_of_events_mc"), 1.5);
      if (std::fabs(collision.posZ()) > 10.0)
        continue;

      registryMC.fill(HIST("number_of_events_mc"), 2.5);
      float multiplicity = collision.centFT0M();

      auto v0s_per_coll = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto casc_per_coll = Cascades.sliceBy(perCollisionCasc, collision.globalIndex());
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, collision.globalIndex());
      auto tracks_per_coll = mcTracks.sliceBy(perCollisionTrk, collision.globalIndex());

      for (const auto& v0 : v0s_per_coll) {

        const auto& pos = v0.posTrack_as<MCTracks>();
        const auto& neg = v0.negTrack_as<MCTracks>();
        if (!pos.has_mcParticle())
          continue;
        if (!neg.has_mcParticle())
          continue;

        auto posParticle = pos.mcParticle_as<aod::McParticles>();
        auto negParticle = neg.mcParticle_as<aod::McParticles>();
        if (!posParticle.has_mothers())
          continue;
        if (!negParticle.has_mothers())
          continue;

        int pdg_parent(0);
        bool isPhysPrim = false;
        for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
          for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
            if (particleMotherOfNeg == particleMotherOfPos) {
              pdg_parent = particleMotherOfNeg.pdgCode();
              isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
            }
          }
        }
        if (pdg_parent == 0)
          continue;

        if (passedK0ShortSelection(v0, pos, neg) && pdg_parent == 310) {
          registryMC.fill(HIST("K0s_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedLambdaSelection(v0, pos, neg) && pdg_parent == 3122) {
          registryMC.fill(HIST("Lambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdg_parent == -3122) {
          registryMC.fill(HIST("AntiLambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (!isPhysPrim)
          continue;

        // Momentum of V0
        TVector3 momentum_pos(posParticle.px(), posParticle.py(), posParticle.pz());
        TVector3 momentum_neg(negParticle.px(), negParticle.py(), negParticle.pz());
        TVector3 momentum_v0 = momentum_pos + momentum_neg;

        double w_k0_jet(1.0), w_k0_ue(1.0), w_lambda_jet(1.0), w_lambda_ue(1.0), w_antilambda_jet(1.0), w_antilambda_ue(1.0);
        if (applyReweighting) {
          int ix = twod_weights_k0_jet->GetXaxis()->FindBin(momentum_v0.Pt());
          int iy = twod_weights_k0_jet->GetYaxis()->FindBin(momentum_v0.Eta());
          w_k0_jet = twod_weights_k0_jet->GetBinContent(ix, iy);
          w_k0_ue = twod_weights_k0_ue->GetBinContent(ix, iy);
          w_lambda_jet = twod_weights_lambda_jet->GetBinContent(ix, iy);
          w_lambda_ue = twod_weights_lambda_ue->GetBinContent(ix, iy);
          w_antilambda_jet = twod_weights_antilambda_jet->GetBinContent(ix, iy);
          w_antilambda_ue = twod_weights_antilambda_ue->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twod_weights_k0_jet->GetNbinsX()) {
            w_k0_jet = 1.0;
            w_k0_ue = 1.0;
            w_lambda_jet = 1.0;
            w_lambda_ue = 1.0;
            w_antilambda_jet = 1.0;
            w_antilambda_ue = 1.0;
          }
          if (iy == 0 || iy > twod_weights_k0_jet->GetNbinsY()) {
            w_k0_jet = 1.0;
            w_k0_ue = 1.0;
            w_lambda_jet = 1.0;
            w_lambda_ue = 1.0;
            w_antilambda_jet = 1.0;
            w_antilambda_ue = 1.0;
          }
        }

        if (passedK0ShortSelection(v0, pos, neg) && pdg_parent == 310) {
          registryMC.fill(HIST("K0s_reconstructed_jet"), multiplicity, v0.pt(), w_k0_jet);
          registryMC.fill(HIST("K0s_reconstructed_ue"), multiplicity, v0.pt(), w_k0_ue);
        }
        if (passedLambdaSelection(v0, pos, neg) && pdg_parent == 3122) {
          registryMC.fill(HIST("Lambda_reconstructed_jet"), multiplicity, v0.pt(), w_lambda_jet);
          registryMC.fill(HIST("Lambda_reconstructed_ue"), multiplicity, v0.pt(), w_lambda_ue);
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdg_parent == -3122) {
          registryMC.fill(HIST("AntiLambda_reconstructed_jet"), multiplicity, v0.pt(), w_antilambda_jet);
          registryMC.fill(HIST("AntiLambda_reconstructed_ue"), multiplicity, v0.pt(), w_antilambda_ue);
        }
      }

      // Cascades
      for (const auto& casc : casc_per_coll) {
        auto bach = casc.template bachelor_as<MCTracks>();
        auto neg = casc.template negTrack_as<MCTracks>();
        auto pos = casc.template posTrack_as<MCTracks>();

        if (!bach.has_mcParticle())
          continue;
        if (!pos.has_mcParticle())
          continue;
        if (!neg.has_mcParticle())
          continue;

        auto posParticle = pos.mcParticle_as<aod::McParticles>();
        auto negParticle = neg.mcParticle_as<aod::McParticles>();
        auto bachParticle = bach.mcParticle_as<aod::McParticles>();
        if (!posParticle.has_mothers())
          continue;
        if (!negParticle.has_mothers())
          continue;
        if (!bachParticle.has_mothers())
          continue;

        int pdg_parent(0);
        for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
          for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
            for (const auto& particleMotherOfBach : bachParticle.mothers_as<aod::McParticles>()) {
              if (particleMotherOfNeg != particleMotherOfPos)
                continue;
              if (std::fabs(particleMotherOfNeg.pdgCode()) != 3122)
                continue;
              if (!particleMotherOfBach.isPhysicalPrimary())
                continue;

              pdg_parent = particleMotherOfBach.pdgCode();
            }
          }
        }
        if (pdg_parent == 0)
          continue;

        // Xi+
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdg_parent == -3312) {
          registryMC.fill(HIST("XiPos_reconstructed"), multiplicity, casc.pt());
        }
        // Xi-
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdg_parent == 3312) {
          registryMC.fill(HIST("XiNeg_reconstructed"), multiplicity, casc.pt());
        }
        // Omega+
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdg_parent == -3334) {
          registryMC.fill(HIST("OmegaPos_reconstructed"), multiplicity, casc.pt());
        }
        // Omega-
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdg_parent == 3334) {
          registryMC.fill(HIST("OmegaNeg_reconstructed"), multiplicity, casc.pt());
        }
      }

      // Reconstructed Tracks
      for (const auto& track : tracks_per_coll) {

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;
        // Track Selection
        if (!passedTrackSelectionForPions(track))
          continue;

        const auto particle = track.mcParticle();
        if (std::fabs(particle.pdgCode()) != 211)
          continue;

        if (particle.isPhysicalPrimary()) {
          if (track.sign() > 0)
            registryMC.fill(HIST("piplus_dcaxy_prim"), multiplicity, track.pt(), track.dcaXY());
          if (track.sign() < 0)
            registryMC.fill(HIST("piminus_dcaxy_prim"), multiplicity, track.pt(), track.dcaXY());
        }
        if (!particle.isPhysicalPrimary()) {
          if (track.sign() > 0)
            registryMC.fill(HIST("piplus_dcaxy_sec"), multiplicity, track.pt(), track.dcaXY());
          if (track.sign() < 0)
            registryMC.fill(HIST("piminus_dcaxy_sec"), multiplicity, track.pt(), track.dcaXY());
        }

        if (std::fabs(track.dcaXY()) > dcaxyMax)
          continue;

        if (track.tpcNSigmaPi() < nsigmaTPCmin || track.tpcNSigmaPi() > nsigmaTPCmax)
          continue;

        if (!particle.isPhysicalPrimary())
          continue;

        double w_piplus_jet(1.0), w_piplus_ue(1.0);
        double w_piminus_jet(1.0), w_piminus_ue(1.0);
        if (applyReweighting) {
          auto getWeight = [&](TH2F* histo) {
            if (!histo)
              return 1.0;
            const int bx = histo->GetXaxis()->FindBin(track.pt());
            const int by = histo->GetYaxis()->FindBin(track.eta());
            if (bx <= 0 || bx > histo->GetNbinsX()) {
              return 1.0;
            }
            if (by <= 0 || by > histo->GetNbinsX()) {
              return 1.0;
            }
            return histo->GetBinContent(bx, by);
          };
          w_piplus_jet = getWeight(twod_weights_piplus_jet);
          w_piplus_ue = getWeight(twod_weights_piplus_ue);
          w_piminus_jet = getWeight(twod_weights_piminus_jet);
          w_piminus_ue = getWeight(twod_weights_piminus_ue);
        }

        if (track.sign() > 0) {
          registryMC.fill(HIST("pi_plus_rec_in_jet_tpc"), multiplicity, track.pt(), w_piplus_jet);
          registryMC.fill(HIST("pi_plus_rec_in_ue_tpc"), multiplicity, track.pt(), w_piplus_ue);
        } else {
          registryMC.fill(HIST("pi_minus_rec_in_jet_tpc"), multiplicity, track.pt(), w_piminus_jet);
          registryMC.fill(HIST("pi_minus_rec_in_ue_tpc"), multiplicity, track.pt(), w_piminus_ue);
        }

        if (!track.hasTOF())
          continue;
        if (track.tofNSigmaPi() < nsigmaTOFminPionMC || track.tofNSigmaPi() > nsigmaTOFmaxPionMC)
          continue;

        if (track.sign() > 0) {
          registryMC.fill(HIST("pi_plus_rec_in_jet_tof"), multiplicity, track.pt(), w_piplus_jet);
          registryMC.fill(HIST("pi_plus_rec_in_ue_tof"), multiplicity, track.pt(), w_piplus_ue);
        } else {
          registryMC.fill(HIST("pi_minus_rec_in_jet_tof"), multiplicity, track.pt(), w_piminus_jet);
          registryMC.fill(HIST("pi_minus_rec_in_ue_tof"), multiplicity, track.pt(), w_piminus_ue);
        }
      }

      for (const auto& mcParticle : mcParticles_per_coll) {

        if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;

        double w_piplus_jet(1.0), w_piplus_ue(1.0);
        double w_piminus_jet(1.0), w_piminus_ue(1.0);
        double w_k0_jet(1.0), w_k0_ue(1.0), w_lambda_jet(1.0), w_lambda_ue(1.0), w_antilambda_jet(1.0), w_antilambda_ue(1.0);
        if (applyReweighting) {
          auto getWeight = [&](TH2F* histo) {
            if (!histo) {
              return 1.0;
            }
            const int bx = histo->GetXaxis()->FindBin(mcParticle.pt());
            const int by = histo->GetYaxis()->FindBin(mcParticle.eta());
            if (bx <= 0 || bx > histo->GetNbinsX()) {
              return 1.0;
            }
            if (by <= 0 || by > histo->GetNbinsX()) {
              return 1.0;
            }
            return histo->GetBinContent(bx, by);
          };
          w_piplus_jet = getWeight(twod_weights_piplus_jet);
          w_piplus_ue = getWeight(twod_weights_piplus_ue);
          w_piminus_jet = getWeight(twod_weights_piminus_jet);
          w_piminus_ue = getWeight(twod_weights_piminus_ue);

          int ix = twod_weights_k0_jet->GetXaxis()->FindBin(mcParticle.pt());
          int iy = twod_weights_k0_jet->GetYaxis()->FindBin(mcParticle.eta());
          w_k0_jet = twod_weights_k0_jet->GetBinContent(ix, iy);
          w_k0_ue = twod_weights_k0_ue->GetBinContent(ix, iy);
          w_lambda_jet = twod_weights_lambda_jet->GetBinContent(ix, iy);
          w_lambda_ue = twod_weights_lambda_ue->GetBinContent(ix, iy);
          w_antilambda_jet = twod_weights_antilambda_jet->GetBinContent(ix, iy);
          w_antilambda_ue = twod_weights_antilambda_ue->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twod_weights_k0_jet->GetNbinsX()) {
            w_k0_jet = 1.0;
            w_k0_ue = 1.0;
            w_lambda_jet = 1.0;
            w_lambda_ue = 1.0;
            w_antilambda_jet = 1.0;
            w_antilambda_ue = 1.0;
          }
          if (iy == 0 || iy > twod_weights_k0_jet->GetNbinsY()) {
            w_k0_jet = 1.0;
            w_k0_ue = 1.0;
            w_lambda_jet = 1.0;
            w_lambda_ue = 1.0;
            w_antilambda_jet = 1.0;
            w_antilambda_ue = 1.0;
          }
        }

        switch (mcParticle.pdgCode()) {
          case 211: // Pi+
            registryMC.fill(HIST("pi_plus_gen_in_jet"), multiplicity, mcParticle.pt(), w_piplus_jet);
            registryMC.fill(HIST("pi_plus_gen_in_ue"), multiplicity, mcParticle.pt(), w_piplus_ue);
            break;
          case -211: // Pi-
            registryMC.fill(HIST("pi_minus_gen_in_jet"), multiplicity, mcParticle.pt(), w_piminus_jet);
            registryMC.fill(HIST("pi_minus_gen_in_ue"), multiplicity, mcParticle.pt(), w_piminus_ue);
            break;
          case 310: // K0s
            registryMC.fill(HIST("K0s_generated_jet"), multiplicity, mcParticle.pt(), w_k0_jet);
            registryMC.fill(HIST("K0s_generated_ue"), multiplicity, mcParticle.pt(), w_k0_ue);
            registryMC.fill(HIST("K0s_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case 3122: // Lambda
            registryMC.fill(HIST("Lambda_generated_jet"), multiplicity, mcParticle.pt(), w_lambda_jet);
            registryMC.fill(HIST("Lambda_generated_ue"), multiplicity, mcParticle.pt(), w_lambda_ue);
            registryMC.fill(HIST("Lambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case -3122: // AntiLambda
            registryMC.fill(HIST("AntiLambda_generated_jet"), multiplicity, mcParticle.pt(), w_antilambda_jet);
            registryMC.fill(HIST("AntiLambda_generated_ue"), multiplicity, mcParticle.pt(), w_antilambda_ue);
            registryMC.fill(HIST("AntiLambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case -3312: // Xi Pos
            registryMC.fill(HIST("XiPos_generated"), multiplicity, mcParticle.pt());
            registryMC.fill(HIST("Xi_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case 3312: // Xi Neg
            registryMC.fill(HIST("XiNeg_generated"), multiplicity, mcParticle.pt());
            registryMC.fill(HIST("AntiXi_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case -3334: // Omega Pos
            registryMC.fill(HIST("OmegaPos_generated"), multiplicity, mcParticle.pt());
            registryMC.fill(HIST("Omega_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case 3334: // Omega Neg
            registryMC.fill(HIST("OmegaNeg_generated"), multiplicity, mcParticle.pt());
            registryMC.fill(HIST("AntiOmega_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processMCefficiency, "Process MC Efficiency", false);

  void processGen(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mcCollisions) {

      registryMC.fill(HIST("number_of_events_mc"), 3.5);

      // Selection on z_{vertex}
      if (std::fabs(mccollision.posZ()) > 10)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 4.5);

      // MC Particles per Collision
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // List of Tracks
      std::vector<TVector3> trk;

      for (const auto& particle : mcParticles_per_coll) {

        int pdg = std::fabs(particle.pdgCode());

        // Select Primary Particles
        double dx = particle.vx() - mccollision.posX();
        double dy = particle.vy() - mccollision.posY();
        double dz = particle.vz() - mccollision.posZ();
        double dcaxy = std::sqrt(dx * dx + dy * dy);
        double dcaz = std::fabs(dz);
        if (particle.pt() < 0.15)
          continue;
        if (dcaxy > (par0 + par1 / particle.pt()))
          continue;
        if (dcaz > (par0 + par1 / particle.pt()))
          continue;
        if (std::fabs(particle.eta()) > 0.8)
          continue;

        // PDG Selection
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
        for (int i = 0; i < static_cast<int>(trk.size()); i++) { // o2-linter: disable=[const-ref-in-for-loop]
          if (trk[i].Mag() == 0)
            continue;
          double diB = 1.0 / (trk[i].Pt() * trk[i].Pt());
          if (diB < diB_min) {
            diB_min = diB;
            iB_min = i;
          }
          for (int j = (i + 1); j < static_cast<int>(trk.size()); j++) { // o2-linter: disable=[const-ref-in-for-loop]
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

        if ((std::fabs(jet[i].Eta()) + Rjet) > etaMax)
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

        for (const auto& particle : mcParticles_per_coll) {

          // Select Primary Particles
          double dx = particle.vx() - mccollision.posX();
          double dy = particle.vy() - mccollision.posY();
          double dz = particle.vz() - mccollision.posZ();
          double dcaxy = std::sqrt(dx * dx + dy * dy);
          double dcaz = std::fabs(dz);
          if (particle.pt() < 0.15)
            continue;
          if (dcaxy > (par0 + par1 / particle.pt()))
            continue;
          if (dcaz > (par0 + par1 / particle.pt()))
            continue;
          if (std::fabs(particle.eta()) > 0.8)
            continue;

          // PDG Selection
          int pdg = std::fabs(particle.pdgCode());
          if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
            continue;

          TVector3 sel_track(particle.px(), particle.py(), particle.pz());

          double deltaEta_jet = sel_track.Eta() - jet[i].Eta();
          double deltaPhi_jet = GetDeltaPhi(sel_track.Phi(), jet[i].Phi());
          double deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          double deltaEta_ue1 = sel_track.Eta() - ue_axis1.Eta();
          double deltaPhi_ue1 = GetDeltaPhi(sel_track.Phi(), ue_axis1.Phi());
          double deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          double deltaEta_ue2 = sel_track.Eta() - ue_axis2.Eta();
          double deltaPhi_ue2 = GetDeltaPhi(sel_track.Phi(), ue_axis2.Phi());
          double deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

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
        for (const auto& particle : mcParticles_per_coll) {

          if (!particle.isPhysicalPrimary())
            continue;

          TVector3 particle_dir(particle.px(), particle.py(), particle.pz());
          double deltaEta_jet = particle_dir.Eta() - jet[i].Eta();
          double deltaPhi_jet = GetDeltaPhi(particle_dir.Phi(), jet[i].Phi());
          double deltaR_jet = std::sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
          double deltaEta_ue1 = particle_dir.Eta() - ue1[i].Eta();
          double deltaPhi_ue1 = GetDeltaPhi(particle_dir.Phi(), ue1[i].Phi());
          double deltaR_ue1 = std::sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
          double deltaEta_ue2 = particle_dir.Eta() - ue2[i].Eta();
          double deltaPhi_ue2 = GetDeltaPhi(particle_dir.Phi(), ue2[i].Phi());
          double deltaR_ue2 = std::sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

          int pdg = particle.pdgCode();

          if (pdg == 211) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("Pion_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("Pion_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == 310) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("K0s_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("K0s_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == 3122) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("Lambda_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("Lambda_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == -3122) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("AntiLambda_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("AntiLambda_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == 3312) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("Xi_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("Xi_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == -3312) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("AntiXi_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("AntiXi_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == 3334) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("Omega_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("Omega_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
          if (pdg == -3334) {
            if (deltaR_jet < Rjet) {
              registryMC.fill(HIST("AntiOmega_eta_pt_jet"), particle.pt(), particle.eta());
            }
            if (deltaR_ue1 < Rjet || deltaR_ue2 < Rjet) {
              registryMC.fill(HIST("AntiOmega_eta_pt_ue"), particle.pt(), particle.eta());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processGen, "Process generated MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangeness_in_jets>(cfgc)};
}
