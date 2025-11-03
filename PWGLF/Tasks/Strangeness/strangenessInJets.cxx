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
/// \file strangenessInJets.cxx
///
/// \brief task for analysis of strangeness in jets
/// \author Alberto Calivà (alberto.caliva@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \author Nicolò Jacazio (nicolo.jacazio@cern.ch)
/// \author Sara Pucillo (sara.pucillo@cern.ch)
///
/// \since May 22, 2024

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <cmath>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;
using std::array;

// Define convenient aliases for joined AOD tables
using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using DaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using DaughterTracksMC = soa::Join<DaughterTracks, aod::McTrackLabels>;

struct StrangenessInJets {

  // Instantiate the CCDB service and API interface
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Instantiate the Zorro processor for skimmed data and define an output object
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Define histogram registries
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global analysis parameters
  Configurable<int> particleOfInterest{"particleOfInterest", 0, "0 = K0 and Lambda, 1 = Xi and Omega, 2 = pion, 3 = kaon, 4 = proton"};
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum reconstructed pt of the jet (GeV/c)"};
  Configurable<double> rJet{"rJet", 0.3, "Jet resolution parameter (R)"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum z-vertex position"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from detector edge"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Enable processing of skimmed data"};
  Configurable<std::string> triggerName{"triggerName", "fOmega", "Software trigger name"};

  // Track analysis parameters
  Configurable<int> minITSnCls{"minITSnCls", 4, "Minimum number of ITS clusters"};
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80, "Minimum number of TPC crossed rows"};
  Configurable<double> maxChi2TPC{"maxChi2TPC", 4.0f, "Maximum chi2 per cluster TPC"};
  Configurable<double> etaMin{"etaMin", -0.8f, "Minimum eta"};
  Configurable<double> etaMax{"etaMax", +0.8f, "Maximum eta"};
  Configurable<double> ptMinV0Proton{"ptMinV0Proton", 0.3f, "Minimum pt of protons from V0"};
  Configurable<double> ptMaxV0Proton{"ptMaxV0Proton", 10.0f, "Maximum pt of protons from V0"};
  Configurable<double> ptMinV0Pion{"ptMinV0Pion", 0.1f, "Minimum pt of pions from V0"};
  Configurable<double> ptMaxV0Pion{"ptMaxV0Pion", 1.5f, "Maximum pt of pions from V0"};
  Configurable<double> ptMinK0Pion{"ptMinK0Pion", 0.3f, "Minimum pt of pions from K0"};
  Configurable<double> ptMaxK0Pion{"ptMaxK0Pion", 10.0f, "Maximum pt of pions from K0"};
  Configurable<double> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<double> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<double> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<double> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<bool> requireITS{"requireITS", false, "Require ITS hit"};
  Configurable<bool> requireTOF{"requireTOF", false, "Require TOF hit"};

  // V0 analysis parameters
  Configurable<double> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<double> maximumV0Radius{"maximumV0Radius", 40.0f, "Maximum V0 Radius"};
  Configurable<double> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA of negative track to primary vertex"};
  Configurable<double> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA of positive track to primary vertex"};
  Configurable<double> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 cosine of pointing angle"};
  Configurable<double> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA between V0 daughters"};

  // Cascade analysis parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 0.1f, "Minimum cascade radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 40.0f, "Maximum cascade radius"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum cascade cosine of pointing angle"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.1f, "Minimum DCA of bachelor to primary vertex"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.1f, "Minimum DCA of V0 to primary vertex"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA between daughters"};
  Configurable<float> deltaMassXi{"deltaMassXi", 0.02f, "Mass window for Xi rejection"};
  Configurable<float> deltaMassOmega{"deltaMassOmega", 0.02f, "Mass window for Omega rejection"};
  Configurable<float> deltaMassLambda{"deltaMassLambda", 0.02f, "Mass window for Lambda inclusion"};

  struct : ConfigurableGroup {
    ConfigurableAxis longLivedBinsNsigma{"longLivedBinsNsigma", {200, -10.f, 10.f}, "Binning of nSigma axis"};
    ConfigurableAxis longLivedBinsPt{"longLivedBinsPt", {VARIABLE_WIDTH, -5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.18, -0.16, -0.14, -0.12, -0.1, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};
    ConfigurableAxis longLivedBinsDca{"longLivedBinsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};
  } longLivedOptions;

  // List of Particles
  enum ParticleOfInterest { kV0Particles,
                            kCascades,
                            kPions,
                            kKaons,
                            kProtons };

  // Instantiate utility class for jet background subtraction
  JetBkgSubUtils backgroundSub;

  // Initialize CCDB access and histogram registry for Zorro processing
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerName.value);
      zorro.populateHistRegistry(registryData, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // Define binning and axis specifications for multiplicity, eta, pT, PID, and invariant mass histograms
    std::vector<double> multBinning = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec multAxis = {multBinning, "FT0C percentile"};
    const AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassK0sAxis{200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassLambdaAxis{200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassXiAxis{200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassOmegaAxis{200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"};
    const AxisSpec ptAxisLongLived{longLivedOptions.longLivedBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nsigmaTOFAxis{longLivedOptions.longLivedBinsNsigma, "n#sigma_{TOF}"};
    const AxisSpec nsigmaTPCAxis{longLivedOptions.longLivedBinsNsigma, "n#sigma_{TPC}"};
    const AxisSpec dcaAxis{longLivedOptions.longLivedBinsDca, "DCA_{xy} (cm)"};

    // Histograms for real data
    if (doprocessData) {

      // Event counters
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
      registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {{101, 0, 101, "Multiplicity percentile"}});

      // Histograms for analysis of strange hadrons
      switch (particleOfInterest) {
        case ParticleOfInterest::kV0Particles:
          registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          break;
        case ParticleOfInterest::kCascades:
          registryData.add("XiPos_in_jet", "XiPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiPos_in_ue", "XiPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_jet", "XiNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_ue", "XiNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("OmegaPos_in_jet", "OmegaPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaPos_in_ue", "OmegaPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_jet", "OmegaNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_ue", "OmegaNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          break;
        case ParticleOfInterest::kPions:
        case ParticleOfInterest::kKaons:
        case ParticleOfInterest::kProtons:
          registryData.add("ll_in_jet", "ll_in_jet", HistType::kTHnSparseF, {multBinning, ptAxisLongLived, nsigmaTPCAxis, nsigmaTOFAxis, dcaAxis});
          registryData.add("ll_in_ue", "ll_in_ue", HistType::kTHnSparseF, {multBinning, ptAxisLongLived, nsigmaTPCAxis, nsigmaTOFAxis, dcaAxis});
          break;
        default:
          LOG(fatal) << "Cannot interpret particle " << particleOfInterest;
          break;
      }
    }

    // Histograms for mc generated
    if (doprocessMCgenerated) {

      // Event counter
      registryMC.add("number_of_events_mc_gen", "number of gen events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});

      // Histograms for analysis
      switch (particleOfInterest) {
        case ParticleOfInterest::kV0Particles:
          registryMC.add("K0s_generated_jet", "K0s_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("K0s_generated_ue", "K0s_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_generated_jet", "Lambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_generated_ue", "Lambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_generated_jet", "AntiLambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_generated_ue", "AntiLambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          break;
        case ParticleOfInterest::kCascades:
          registryMC.add("XiPos_generated_jet", "XiPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiPos_generated_ue", "XiPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiNeg_generated_jet", "XiNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiNeg_generated_ue", "XiNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaPos_generated_jet", "OmegaPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaPos_generated_ue", "OmegaPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaNeg_generated_jet", "OmegaNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaNeg_generated_ue", "OmegaNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
          break;
        case ParticleOfInterest::kPions:
        case ParticleOfInterest::kKaons:
        case ParticleOfInterest::kProtons:
          registryMC.add("ll_generated_in_jet", "ll_generated_in_jet", HistType::kTH2F, {multBinning, ptAxisLongLived});
          registryMC.add("ll_generated_in_ue", "ll_generated_in_ue", HistType::kTH2F, {multBinning, ptAxisLongLived});
          break;
        default:
          LOG(fatal) << "Cannot interpret particle " << particleOfInterest;
      }
    }

    // Histograms for mc reconstructed
    if (doprocessMCreconstructed) {

      // Event counter
      registryMC.add("number_of_events_mc_rec", "number of rec events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});

      // Histograms for analysis
      switch (particleOfInterest) {
        case ParticleOfInterest::kV0Particles:
          registryMC.add("K0s_reconstructed_jet", "K0s_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("K0s_reconstructed_ue", "K0s_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_reconstructed_jet", "Lambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_reconstructed_ue", "Lambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_reconstructed_jet", "AntiLambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_reconstructed_ue", "AntiLambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          // Histograms for secondary hadrons
          registryMC.add("K0s_reconstructed_jet_incl", "K0s_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("K0s_reconstructed_ue_incl", "K0s_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_reconstructed_jet_incl", "Lambda_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("Lambda_reconstructed_ue_incl", "Lambda_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_reconstructed_jet_incl", "AntiLambda_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("AntiLambda_reconstructed_ue_incl", "AntiLambda_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
          break;
        case ParticleOfInterest::kCascades:
          registryMC.add("XiPos_reconstructed_jet", "XiPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiPos_reconstructed_ue", "XiPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiNeg_reconstructed_jet", "XiNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("XiNeg_reconstructed_ue", "XiNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaPos_reconstructed_jet", "OmegaPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaPos_reconstructed_ue", "OmegaPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaNeg_reconstructed_jet", "OmegaNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
          registryMC.add("OmegaNeg_reconstructed_ue", "OmegaNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
          break;
        case ParticleOfInterest::kPions:
        case ParticleOfInterest::kKaons:
        case ParticleOfInterest::kProtons:
          registryMC.add("ll_reconstructed_in_jet", "ll_reconstructed_in_jet", HistType::kTH2F, {multBinning, ptAxisLongLived});
          registryMC.add("ll_reconstructed_in_ue", "ll_reconstructed_in_ue", HistType::kTH2F, {multBinning, ptAxisLongLived});
          break;
        default:
          LOG(fatal) << "Cannot interpret particle " << particleOfInterest;
      }
    }
  }

  /*
  // Calculation of perpendicular axes
  void getPerpendicularAxis(TVector3 p, TVector3& u, double sign)
  {
    // initialization
    double ux(0), uy(0), uz(0);

    // components of vector p
    const double px = p.X();
    const double py = p.Y();
    const double pz = p.Z();

    // protection 1
    if (px == 0 && py != 0) {
      uy = -(pz * pz) / py;
      ux = sign * std::sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // protection 2
    if (py == 0 && px != 0) {
      ux = -(pz * pz) / px;
      uy = sign * std::sqrt(px * px - (pz * pz * pz * pz) / (px * px));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // equation parameters
    const double a = px * px + py * py;
    const double b = 2.0 * px * pz * pz;
    const double c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    const double delta = b * b - 4.0 * a * c;

    // protection agains delta<0
    if (delta < 0) {
      return;
    }

    // solutions
    ux = (-b + sign * std::sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }
  */

  // Delta phi calculation
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

  // Check if particle is a physical primary or a decay product of a heavy-flavor hadron
  bool isPhysicalPrimaryOrFromHF(aod::McParticle const& particle, aod::McParticles const& mcParticles)
  {
    // Keep only pi, K, p, e, mu
    int pdg = std::abs(particle.pdgCode());
    if (!(pdg == PDG_t::kPiPlus || pdg == PDG_t::kKPlus || pdg == PDG_t::kProton || pdg == PDG_t::kElectron || pdg == PDG_t::kMuonMinus))
      return false;

    // Constants for identifying heavy-flavor (charm and bottom) content from PDG codes
    static constexpr int kCharmQuark = 4;
    static constexpr int kBottomQuark = 5;
    static constexpr int hundreds = 100;
    static constexpr int thousands = 1000;

    // Check if particle is from heavy-flavor decay
    bool fromHF = false;
    if (particle.has_mothers()) {
      auto mother = mcParticles.iteratorAt(particle.mothersIds()[0]);
      int motherPdg = std::abs(mother.pdgCode());
      fromHF = (motherPdg / hundreds == kCharmQuark || motherPdg / hundreds == kBottomQuark || motherPdg / thousands == kCharmQuark || motherPdg / thousands == kBottomQuark);
    }

    // Select only physical primary particles or from heavy-flavor
    return (particle.isPhysicalPrimary() || fromHF);
  }

  // Compute two transverse directions orthogonal to vector p
  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    // Get momentum components
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    // Precompute squared terms
    double px2 = px * px;
    double py2 = py * py;
    double pz2 = pz * pz;
    double pz4 = pz2 * pz2;

    // Case 1: vector along z-axis -> undefined perpendiculars
    if (px == 0 && py == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    // Case 2: px = 0 -> avoid division by zero
    if (px == 0 && py != 0) {
      double ux = std::sqrt(py2 - pz4 / py2);
      double uy = -pz2 / py;
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(-ux, uy, pz);
      return;
    }

    // Case 3: py = 0 -> avoid division by zero
    if (py == 0 && px != 0) {
      double ux = -pz2 / px;
      double uy = std::sqrt(px2 - pz4 / px2);
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(ux, -uy, pz);
      return;
    }

    // General case: solve quadratic for perpendicular vectors
    double a = px2 + py2;
    double b = 2.0 * px * pz2;
    double c = pz4 - py2 * py2 - px2 * py2;
    double delta = b * b - 4.0 * a * c;

    // Invalid or degenerate solutions
    if (delta < 0 || a == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    // Solution 1
    double u1x = (-b + std::sqrt(delta)) / (2.0 * a);
    double u1y = (-pz2 - px * u1x) / py;
    u1.SetXYZ(u1x, u1y, pz);

    // Solution 2
    double u2x = (-b - std::sqrt(delta)) / (2.0 * a);
    double u2y = (-pz2 - px * u2x) / py;
    u2.SetXYZ(u2x, u2y, pz);
  }

  // Find ITS hit
  template <typename TrackIts>
  bool hasITSHitOnLayer(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-track selection for particles inside jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    const int minTpcCr = 70;
    const double maxChi2Tpc = 4.0;
    const double maxChi2Its = 36.0;
    const double maxPseudorapidity = 0.8;
    const double minPtTrack = 0.1;
    const double dcaxyMaxTrackPar0 = 0.0105;
    const double dcaxyMaxTrackPar1 = 0.035;
    const double dcaxyMaxTrackPar2 = 1.1;
    const double dcazMaxTrack = 2.0;

    if (!track.hasITS())
      return false;
    if ((!hasITSHitOnLayer(track, 1)) && (!hasITSHitOnLayer(track, 2)) && (!hasITSHitOnLayer(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcCr)
      return false;
    if (track.tpcChi2NCl() > maxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > maxChi2Its)
      return false;
    if (std::fabs(track.eta()) > maxPseudorapidity)
      return false;
    if (track.pt() < minPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (dcaxyMaxTrackPar0 + dcaxyMaxTrackPar1 / std::pow(track.pt(), dcaxyMaxTrackPar2)))
      return false;
    if (std::fabs(track.dcaZ()) > dcazMaxTrack)
      return false;
    return true;
  }

  // Lambda selections
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-track selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum of lambda daughters
    TVector3 proton(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selection on pt of Lambda daughters
    if (proton.Pt() < ptMinV0Proton || proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion || pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC): positive track = proton, negative track = pion
    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF): positive track = proton, negative track = pion
    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // AntiLambda selections
  template <typename AntiLambda, typename TrackPos, typename TrackNeg>
  bool passedAntiLambdaSelection(const AntiLambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-track selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum AntiLambda daughters
    TVector3 pion(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 proton(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selections on pt of Antilambda daughters
    if (proton.Pt() < ptMinV0Proton || proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion || pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC): negative track = proton, positive track = pion
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    // PID selections (TOF): negative track = proton, positive track = pion
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // K0s selections
  template <typename K0short, typename TrackPos, typename TrackNeg>
  bool passedK0ShortSelection(const K0short& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum of K0s daughters
    TVector3 pionPos(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pionNeg(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selections on pt of K0s daughters
    if (pionPos.Pt() < ptMinK0Pion || pionPos.Pt() > ptMaxK0Pion)
      return false;
    if (pionNeg.Pt() < ptMinK0Pion || pionNeg.Pt() > ptMaxK0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
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
    // Single-track selections on cascade daughters
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Xi+ selection (Xi+ -> antiL + pi+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMinV0Proton || ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion || ptrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), o2::constants::physics::MassProton);
      pPion.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), o2::constants::physics::MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - o2::constants::physics::MassLambda0) > deltaMassLambda)
        return false;
    }

    // Xi- selection (Xi- -> L + pi-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMinV0Proton || ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion || ntrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), o2::constants::physics::MassProton);
      pPion.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), o2::constants::physics::MassPionCharged);
      const double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - o2::constants::physics::MassLambda0) > deltaMassLambda)
        return false;
    }

    // V0 selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // PID selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin || btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaPi() < nsigmaTOFmin || btrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Reject candidates compatible with Omega
    if (std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < deltaMassOmega)
      return false;
    return true;
  }

  // Omega selections
  template <typename Omega, typename TrackPos, typename TrackNeg, typename TrackBac, typename Coll>
  bool passedOmegaSelection(const Omega& casc, const TrackPos& ptrack, const TrackNeg& ntrack, const TrackBac& btrack, const Coll& coll)
  {
    // Single-track selections on cascade daughters
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Omega+ selection (Omega+ -> antiL + K+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMinV0Proton || ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion || ptrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), o2::constants::physics::MassProton);
      pPion.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), o2::constants::physics::MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - o2::constants::physics::MassLambda0) > deltaMassLambda)
        return false;
    }

    // Omega- selection (Omega- -> L + K-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMinV0Proton || ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion || ntrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), o2::constants::physics::MassProton);
      pPion.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), o2::constants::physics::MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - o2::constants::physics::MassLambda0) > deltaMassLambda)
        return false;
    }

    // V0 selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // PID selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin || btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaKa() < nsigmaTOFmin || btrack.tofNSigmaKa() > nsigmaTOFmax)
        return false;
    }

    // Reject candidates compatible with Xi
    if (std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < deltaMassXi)
      return false;
    return true;
  }

  // Single-track selection
  template <typename Track>
  bool passedSingleTrackSelection(const Track& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
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

  // Process data
  void processData(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s,
                   aod::CascDataExt const& Cascades, DaughterTracks const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    // Fill event counter before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Get the bunch crossing (BC) information associated with the collision
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

    // Initialize CCDB objects using the BC info
    initCCDB(bc);

    // If skimmed processing is enabled, skip this event unless it passes Zorro selection
    if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) {
      return;
    }

    // Fill event counter after zorro selection
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Fill event counter after event selection
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // Loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {

      // Require that tracks pass selection criteria
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }

    // Reject empty events
    if (fjParticles.size() < 1)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

    // Jet selection
    bool isAtLeastOneJetSelected = false;
    std::vector<TVector3> selectedJet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    // Loop over reconstructed jets
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
        continue;

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      if (jetMinusBkg.pt() < minJetPt)
        continue;
      isAtLeastOneJetSelected = true;

      // Calculation of perpendicular cones
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
        continue;
      }

      // Store jet and UE axes
      selectedJet.emplace_back(jetAxis);
      ue1.emplace_back(ueAxis1);
      ue2.emplace_back(ueAxis2);
    }
    if (!isAtLeastOneJetSelected)
      return;

    // Fill event counter with events with at least one jet
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Event multiplicity
    const float multiplicity = collision.centFT0M();

    // Fill event multiplicity
    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);

    // Loop over selected jets
    for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {
      switch (particleOfInterest) {
        case ParticleOfInterest::kV0Particles: // V0s
          for (const auto& v0 : fullV0s) {

            // Get V0 daughters
            const auto& pos = v0.posTrack_as<DaughterTracks>();
            const auto& neg = v0.negTrack_as<DaughterTracks>();
            TVector3 v0dir(v0.px(), v0.py(), v0.pz());

            // Calculate distance from jet and UE axes
            const float deltaEtaJet = v0dir.Eta() - selectedJet[i].Eta();
            const float deltaPhiJet = getDeltaPhi(v0dir.Phi(), selectedJet[i].Phi());
            const float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            const float deltaEtaUe1 = v0dir.Eta() - ue1[i].Eta();
            const float deltaPhiUe1 = getDeltaPhi(v0dir.Phi(), ue1[i].Phi());
            const float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            const float deltaEtaUe2 = v0dir.Eta() - ue2[i].Eta();
            const float deltaPhiUe2 = getDeltaPhi(v0dir.Phi(), ue2[i].Phi());
            const float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // K0s
            if (passedK0ShortSelection(v0, pos, neg)) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("K0s_in_jet"), multiplicity, v0.pt(), v0.mK0Short());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("K0s_in_ue"), multiplicity, v0.pt(), v0.mK0Short());
              }
            }
            // Lambda
            if (passedLambdaSelection(v0, pos, neg)) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("Lambda_in_jet"), multiplicity, v0.pt(), v0.mLambda());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("Lambda_in_ue"), multiplicity, v0.pt(), v0.mLambda());
              }
            }
            // AntiLambda
            if (passedAntiLambdaSelection(v0, pos, neg)) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("AntiLambda_in_jet"), multiplicity, v0.pt(), v0.mAntiLambda());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("AntiLambda_in_ue"), multiplicity, v0.pt(), v0.mAntiLambda());
              }
            }
          }
          break;
        case ParticleOfInterest::kCascades: // Cascades
          for (const auto& casc : Cascades) {
            // Get cascade daughters
            const auto& bach = casc.bachelor_as<DaughterTracks>();
            const auto& pos = casc.posTrack_as<DaughterTracks>();
            const auto& neg = casc.negTrack_as<DaughterTracks>();
            TVector3 cascadeDir(casc.px(), casc.py(), casc.pz());

            // Calculate distance from jet and UE axes
            const double deltaEtaJet = cascadeDir.Eta() - selectedJet[i].Eta();
            const double deltaPhiJet = getDeltaPhi(cascadeDir.Phi(), selectedJet[i].Phi());
            const double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            const double deltaEtaUe1 = cascadeDir.Eta() - ue1[i].Eta();
            const double deltaPhiUe1 = getDeltaPhi(cascadeDir.Phi(), ue1[i].Phi());
            const double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            const double deltaEtaUe2 = cascadeDir.Eta() - ue2[i].Eta();
            const double deltaPhiUe2 = getDeltaPhi(cascadeDir.Phi(), ue2[i].Phi());
            const double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // Xi+
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("XiPos_in_jet"), multiplicity, casc.pt(), casc.mXi());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("XiPos_in_ue"), multiplicity, casc.pt(), casc.mXi());
              }
            }
            // Xi-
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("XiNeg_in_jet"), multiplicity, casc.pt(), casc.mXi());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("XiNeg_in_ue"), multiplicity, casc.pt(), casc.mXi());
              }
            }
            // Omega+
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("OmegaPos_in_jet"), multiplicity, casc.pt(), casc.mOmega());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("OmegaPos_in_ue"), multiplicity, casc.pt(), casc.mOmega());
              }
            }
            // Omega-
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
              if (deltaRjet < rJet) {
                registryData.fill(HIST("OmegaNeg_in_jet"), multiplicity, casc.pt(), casc.mOmega());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryData.fill(HIST("OmegaNeg_in_ue"), multiplicity, casc.pt(), casc.mOmega());
              }
            }
          }
          break;
        case ParticleOfInterest::kPions:
        case ParticleOfInterest::kKaons:
        case ParticleOfInterest::kProtons:
          for (const auto& trk : tracks) {

            if (!passedSingleTrackSelection(trk)) {
              continue;
            }

            const double deltaEtaJet = trk.eta() - selectedJet[i].Eta();
            const double deltaPhiJet = getDeltaPhi(trk.phi(), selectedJet[i].Phi());
            const double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            const double deltaEtaUe1 = trk.eta() - ue1[i].Eta();
            const double deltaPhiUe1 = getDeltaPhi(trk.phi(), ue1[i].Phi());
            const double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            const double deltaEtaUe2 = trk.eta() - ue2[i].Eta();
            const double deltaPhiUe2 = getDeltaPhi(trk.phi(), ue2[i].Phi());
            const double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            float nsigmaTPC = 0.f;
            float nsigmaTOF = 0.f;
            switch (particleOfInterest) {
              case ParticleOfInterest::kPions:
                nsigmaTPC = trk.tpcNSigmaPi();
                nsigmaTOF = trk.tofNSigmaPi();
                break;
              case ParticleOfInterest::kKaons:
                nsigmaTPC = trk.tpcNSigmaKa();
                nsigmaTOF = trk.tofNSigmaKa();
                break;
              case ParticleOfInterest::kProtons:
                nsigmaTPC = trk.tpcNSigmaPr();
                nsigmaTOF = trk.tofNSigmaPr();
                break;
            }

            if (deltaRjet < rJet) {
              registryData.fill(HIST("ll_in_jet"), multiplicity, trk.pt() * trk.sign(), nsigmaTPC, nsigmaTOF, trk.dcaXY());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("ll_in_ue"), multiplicity, trk.pt() * trk.sign(), nsigmaTPC, nsigmaTOF, trk.dcaXY());
            }
          }
          break;
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processData, "Process data", true);

  // Define per-collision preslices for V0s, cascades, MC particles, and daughter tracks
  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDataExt> perCollisionCasc = o2::aod::cascade::collisionId;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<DaughterTracksMC> perCollisionTrk = o2::aod::track::collisionId;

  // Generated MC events
  void processMCgenerated(soa::Join<aod::McCollisions, aod::McCentFT0Ms> const& collisions, aod::McParticles const& mcParticles)
  {
    // Define per-event particle containers
    std::vector<fastjet::PseudoJet> fjParticles;
    std::vector<TVector3> strHadronMomentum;
    std::vector<int> pdg;

    // Jet and area definitions
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));

    // Loop over all simulated collision events
    for (const auto& collision : collisions) {

      // Clear containers at the start of the event loop
      fjParticles.clear();
      strHadronMomentum.clear();
      pdg.clear();

      // Fill event counter before any selection
      registryMC.fill(HIST("number_of_events_mc_gen"), 0.5);

      // Need to apply event selection to simulated events
      registryMC.fill(HIST("number_of_events_mc_gen"), 1.5);

      // Require vertex position within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Fill event counter after selection on z-vertex
      registryMC.fill(HIST("number_of_events_mc_gen"), 2.5);

      // Multiplicity of generated event
      double genMultiplicity = collision.centFT0M();

      // MC particles per collision
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, collision.globalIndex());

      // Loop over all MC particles and select physical primaries within acceptance
      for (const auto& particle : mcParticlesPerColl) {

        // Select physical primary particles or HF decay products
        if (!isPhysicalPrimaryOrFromHF(particle, mcParticles))
          continue;

        double minPtParticle = 0.1;
        if (particle.eta() < etaMin || particle.eta() > etaMax || particle.pt() < minPtParticle)
          continue;

        // Build 4-momentum assuming charged pion mass
        static constexpr float kMassPionChargedSquared = o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged;
        const double energy = std::sqrt(particle.p() * particle.p() + kMassPionChargedSquared);
        fastjet::PseudoJet fourMomentum(particle.px(), particle.py(), particle.pz(), energy);
        fourMomentum.set_user_index(particle.pdgCode());
        fjParticles.emplace_back(fourMomentum);

        // Store properties of strange hadrons
        int pdgAbs = std::abs(particle.pdgCode());
        if (pdgAbs == kK0Short || pdgAbs == kLambda0 || pdgAbs == kXiMinus || pdgAbs == kOmegaMinus) {
          pdg.emplace_back(particle.pdgCode());
          strHadronMomentum.emplace_back(particle.px(), particle.py(), particle.pz());
        }
      }

      // Skip events with no particles
      if (fjParticles.size() < 1)
        continue;
      registryMC.fill(HIST("number_of_events_mc_gen"), 3.5);

      // Cluster MC particles into jets using anti-kt algorithm
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

      // Estimate background energy density (rho) in perpendicular cone
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

      // Loop over clustered jets
      for (const auto& jet : jets) {

        // Jet must be fully contained in acceptance
        if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
          continue;

        // Subtract background energy from jet
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);

        // Apply jet pT threshold
        if (jetMinusBkg.pt() < minJetPt)
          continue;
        registryMC.fill(HIST("number_of_events_mc_gen"), 4.5);

        // Set up two perpendicular cone axes for underlying event estimation
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        double coneRadius = std::sqrt(jet.area() / PI);
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
        if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
          continue;
        }

        // Loop over strange hadrons
        int index = -1;
        for (const auto& hadron : strHadronMomentum) {

          // Particle index
          index++;

          // Compute distance of particles from jet and UE axes
          double deltaEtaJet = hadron.Eta() - jetAxis.Eta();
          double deltaPhiJet = getDeltaPhi(hadron.Phi(), jetAxis.Phi());
          double deltaRJet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = hadron.Eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(hadron.Phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = hadron.Eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(hadron.Phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Select particles inside jet
          if (deltaRJet < coneRadius) {
            switch (pdg[index]) {
              case kK0Short:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("K0s_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kLambda0:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("Lambda_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kLambda0Bar:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("AntiLambda_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kXiMinus:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("XiNeg_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kXiPlusBar:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("XiPos_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kOmegaMinus:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("OmegaNeg_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kOmegaPlusBar:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("OmegaPos_generated_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kPiPlus:
                if (particleOfInterest == ParticleOfInterest::kPions) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kKPlus:
                if (particleOfInterest == ParticleOfInterest::kKaons) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kProton:
                if (particleOfInterest == ParticleOfInterest::kProtons) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt());
                }
                break;
              case kPiMinus:
                if (particleOfInterest == ParticleOfInterest::kPions) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              case kKMinus:
                if (particleOfInterest == ParticleOfInterest::kKaons) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              case kProtonBar:
                if (particleOfInterest == ParticleOfInterest::kProtons) {
                  registryMC.fill(HIST("ll_generated_in_jet"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              default:
                break;
            }
          }

          // Select particles inside UE cones
          if (deltaRUe1 < coneRadius || deltaRUe2 < coneRadius) {
            switch (pdg[index]) {
              case kK0Short:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("K0s_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kLambda0:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("Lambda_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kLambda0Bar:
                if (particleOfInterest == ParticleOfInterest::kV0Particles) {
                  registryMC.fill(HIST("AntiLambda_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kXiMinus:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("XiNeg_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kXiPlusBar:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("XiPos_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kOmegaMinus:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("OmegaNeg_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kOmegaPlusBar:
                if (particleOfInterest == ParticleOfInterest::kCascades) {
                  registryMC.fill(HIST("OmegaPos_generated_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kPiPlus:
                if (particleOfInterest == ParticleOfInterest::kPions) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kKPlus:
                if (particleOfInterest == ParticleOfInterest::kKaons) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kProton:
                if (particleOfInterest == ParticleOfInterest::kProtons) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt());
                }
                break;
              case kPiMinus:
                if (particleOfInterest == ParticleOfInterest::kPions) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              case kKMinus:
                if (particleOfInterest == ParticleOfInterest::kKaons) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              case kProtonBar:
                if (particleOfInterest == ParticleOfInterest::kProtons) {
                  registryMC.fill(HIST("ll_generated_in_ue"), genMultiplicity, hadron.Pt() * -1.f);
                }
                break;
              default:
                break;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processMCgenerated, "process generated events", false);

  // Reconstructed MC events
  void processMCreconstructed(SimCollisions const& collisions, soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&,
                              DaughterTracksMC const& mcTracks, aod::V0Datas const& fullV0s,
                              aod::CascDataExt const& Cascades, const aod::McParticles&)
  {
    // Define per-event containers
    std::vector<fastjet::PseudoJet> fjParticles;
    std::vector<TVector3> selectedJet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    // Jet and area definitions
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));

    // Loop over reconstructed collisions
    for (const auto& collision : collisions) {

      if (!collision.has_mcCollision()) {
        continue;
      }

      const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();

      // Clear containers at the start of the event loop
      fjParticles.clear();
      selectedJet.clear();
      ue1.clear();
      ue2.clear();

      // Fill event counter before any selection
      registryMC.fill(HIST("number_of_events_mc_rec"), 0.5);
      if (!collision.sel8())
        continue;

      // Fill event counter after event selection
      registryMC.fill(HIST("number_of_events_mc_rec"), 1.5);
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Fill event counter after selection on z-vertex
      registryMC.fill(HIST("number_of_events_mc_rec"), 2.5);

      // Event multiplicity
      const float multiplicity = mcCollision.centFT0M();

      // Number of V0 and cascades per collision
      auto v0sPerColl = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto cascPerColl = Cascades.sliceBy(perCollisionCasc, collision.globalIndex());
      auto tracksPerColl = mcTracks.sliceBy(perCollisionTrk, collision.globalIndex());

      // Loop over reconstructed tracks
      for (auto const& track : tracksPerColl) {
        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        // 4-momentum representation of a particle
        fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
        fjParticles.emplace_back(fourMomentum);
      }

      // Reject empty events
      if (fjParticles.size() < 1)
        continue;
      registryMC.fill(HIST("number_of_events_mc_rec"), 3.5);

      // Cluster particles using the anti-kt algorithm
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

      // Jet selection
      bool isAtLeastOneJetSelected = false;

      // Loop over clustered jets
      for (const auto& jet : jets) {

        // jet must be fully contained in the acceptance
        if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
          continue;

        // jet pt must be larger than threshold
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
        if (jetMinusBkg.pt() < minJetPt)
          continue;
        isAtLeastOneJetSelected = true;

        // Perpendicular cones
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
        if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
          continue;
        }

        // Store selected jet and UE cone axes
        selectedJet.emplace_back(jetAxis);
        ue1.emplace_back(ueAxis1);
        ue2.emplace_back(ueAxis2);
      }
      if (!isAtLeastOneJetSelected)
        continue;

      // Fill event counter for events with at least one selected jet
      registryMC.fill(HIST("number_of_events_mc_rec"), 4.5);

      // Loop over selected jets
      for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {

        // V0 particles
        if (particleOfInterest == ParticleOfInterest::kV0Particles) {
          for (const auto& v0 : v0sPerColl) {
            const auto& pos = v0.posTrack_as<DaughterTracksMC>();
            const auto& neg = v0.negTrack_as<DaughterTracksMC>();
            TVector3 v0dir(v0.px(), v0.py(), v0.pz());

            // Get MC particles
            if (!pos.has_mcParticle() || !neg.has_mcParticle())
              continue;
            auto posParticle = pos.mcParticle_as<aod::McParticles>();
            auto negParticle = neg.mcParticle_as<aod::McParticles>();
            if (!posParticle.has_mothers() || !negParticle.has_mothers())
              continue;

            // Select particles originating from the same parent
            int pdgParent(0);
            bool isPhysPrim = false;
            for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
              for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos) {
                  pdgParent = particleMotherOfNeg.pdgCode();
                  isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
                }
              }
            }
            if (pdgParent == 0)
              continue;

            // Compute distance from jet and UE axes
            double deltaEtaJet = v0dir.Eta() - selectedJet[i].Eta();
            double deltaPhiJet = getDeltaPhi(v0dir.Phi(), selectedJet[i].Phi());
            double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            double deltaEtaUe1 = v0dir.Eta() - ue1[i].Eta();
            double deltaPhiUe1 = getDeltaPhi(v0dir.Phi(), ue1[i].Phi());
            double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            double deltaEtaUe2 = v0dir.Eta() - ue2[i].Eta();
            double deltaPhiUe2 = getDeltaPhi(v0dir.Phi(), ue2[i].Phi());
            double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // K0s
            if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_ue"), multiplicity, v0.pt());
              }
            }
            // Lambda
            if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0 && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_ue"), multiplicity, v0.pt());
              }
            }
            // AntiLambda
            if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_ue"), multiplicity, v0.pt());
              }
            }

            // Fill inclusive spectra
            // K0s
            if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
            // Lambda
            if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
            // AntiLambda
            if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
          }
        }

        // Cascades
        if (particleOfInterest == ParticleOfInterest::kCascades) {
          for (const auto& casc : cascPerColl) {
            auto bach = casc.bachelor_as<DaughterTracksMC>();
            auto pos = casc.posTrack_as<DaughterTracksMC>();
            auto neg = casc.negTrack_as<DaughterTracksMC>();

            // Get MC particles
            if (!bach.has_mcParticle() || !pos.has_mcParticle() || !neg.has_mcParticle())
              continue;
            auto posParticle = pos.mcParticle_as<aod::McParticles>();
            auto negParticle = neg.mcParticle_as<aod::McParticles>();
            auto bachParticle = bach.mcParticle_as<aod::McParticles>();
            if (!posParticle.has_mothers() || !negParticle.has_mothers() || !bachParticle.has_mothers())
              continue;

            // Select particles originating from the same parent
            int pdgParent(0);
            bool isPhysPrim = false;
            for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
              for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
                for (const auto& particleMotherOfBach : bachParticle.mothers_as<aod::McParticles>()) {
                  if (particleMotherOfNeg != particleMotherOfPos)
                    continue;
                  if (std::abs(particleMotherOfNeg.pdgCode()) != kLambda0)
                    continue;
                  isPhysPrim = particleMotherOfBach.isPhysicalPrimary();
                  pdgParent = particleMotherOfBach.pdgCode();
                }
              }
            }
            if (pdgParent == 0)
              continue;
            if (!isPhysPrim)
              continue;

            // Compute distances from jet and UE axes
            TVector3 cascadeDir(casc.px(), casc.py(), casc.pz());
            double deltaEtaJet = cascadeDir.Eta() - selectedJet[i].Eta();
            double deltaPhiJet = getDeltaPhi(cascadeDir.Phi(), selectedJet[i].Phi());
            double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            double deltaEtaUe1 = cascadeDir.Eta() - ue1[i].Eta();
            double deltaPhiUe1 = getDeltaPhi(cascadeDir.Phi(), ue1[i].Phi());
            double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            double deltaEtaUe2 = cascadeDir.Eta() - ue2[i].Eta();
            double deltaPhiUe2 = getDeltaPhi(cascadeDir.Phi(), ue2[i].Phi());
            double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // Xi+
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() > 0 && pdgParent == kXiPlusBar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("XiPos_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("XiPos_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Xi-
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() < 0 && pdgParent == kXiMinus) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("XiNeg_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("XiNeg_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Omega+
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() > 0 && pdgParent == kOmegaPlusBar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("OmegaPos_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("OmegaPos_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Omega-
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() < 0 && pdgParent == kOmegaMinus) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("OmegaNeg_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("OmegaNeg_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processMCreconstructed, "process reconstructed events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StrangenessInJets>(cfgc)};
}
