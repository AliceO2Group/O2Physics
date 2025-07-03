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
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch), Nicolò Jacazio (nicolo.jacazio@cern.ch), Sara Pucillo (sara.pucillo@cern.ch)
/// \since May 22, 2024

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using StrHadronDaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                                          aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                          aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MCTracks = soa::Join<StrHadronDaughterTracks, aod::McTrackLabels>;

struct StrangenessInJets {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global Parameters
  Configurable<int> particleOfInterest{"particleOfInterest", 0, "0=v0, 1=cascade, 2=pions, 3=k+-, 4=proton"};
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet"};
  Configurable<double> rJet{"rJet", 0.3, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};

  // Axis parameters
  struct : ConfigurableGroup {
    ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};

    ConfigurableAxis binsDca{"binsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};

  } binsOptions;

  // Track Parameters
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> minTpcNcrossedRowsOverFindable{"minTpcNcrossedRowsOverFindable", 0.8, "crossed rows/findable"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> ptMinV0Proton{"ptMinV0Proton", 0.3f, "pt min of proton from V0"};
  Configurable<float> ptMaxV0Proton{"ptMaxV0Proton", 10.0f, "pt max of proton from V0"};
  Configurable<float> ptMinV0Pion{"ptMinV0Pion", 0.1f, "pt min of pion from V0"};
  Configurable<float> ptMaxV0Pion{"ptMaxV0Pion", 1.5f, "pt max of pion from V0"};
  Configurable<float> ptMinK0Pion{"ptMinK0Pion", 0.3f, "pt min of pion from K0"};
  Configurable<float> ptMaxK0Pion{"ptMaxK0Pion", 10.0f, "pt max of pion from K0"};
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
  Configurable<std::string> urlToCcdb{"urlToCcdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> pathToFile{"pathToFile", "", "path to file with reweighting"};
  Configurable<std::string> histoNameWeightPiplusJet{"histoNameWeightPiplusJet", "", "reweighting histogram: Pi Plus in jet"};
  Configurable<std::string> histoNameWeightPiplusUe{"histoNameWeightPiplusUe", "", "reweighting histogram: Pi Plus in ue"};
  Configurable<std::string> histoNameWeightPiminusJet{"histoNameWeightPiminusJet", "", "reweighting histogram: Pi Minus in jet"};
  Configurable<std::string> histoNameWeightPiminusUe{"histoNameWeightPiminusUe", "", "reweighting histogram: Pi Minus in ue"};
  Configurable<std::string> histoNameWeightK0Jet{"histoNameWeightK0Jet", "", "reweighting histogram: K0 in jet"};
  Configurable<std::string> histoNameWeightK0Ue{"histoNameWeightK0Ue", "", "reweighting histogram: K0 in ue"};
  Configurable<std::string> histoNameWeightLambdaJet{"histoNameWeightLambdaJet", "", "reweighting histogram: lambda in jet"};
  Configurable<std::string> histoNameWeightLambdaUe{"histoNameWeightLambdaUe", "", "reweighting histogram: lambda in ue"};
  Configurable<std::string> histoNameWeightAntilambdaJet{"histoNameWeightAntilambdaJet", "", "reweighting histogram: antilambda in jet"};
  Configurable<std::string> histoNameWeightAntilambdaUe{"histoNameWeightAntilambdaUe", "", "reweighting histogram: antilambda in ue"};
  Configurable<std::string> histoNameWeightsXiInJet{"histoNameWeightsXiInJet", "", "reweighting histogram: xi in jet"};
  Configurable<std::string> histoNameWeightsXiUe{"histoNameWeightsXiUe", "", "reweighting histogram: xi in ue"};
  Configurable<std::string> histoNameWeightsAntiXiInJet{"histoNameWeightsAntiXiInJet", "", "reweighting histogram: antixi in jet"};
  Configurable<std::string> histoNameWeightsAntiXiUe{"histoNameWeightsAntiXiUe", "", "reweighting histogram: antixi in ue"};
  Configurable<std::string> histoNameWeightsOmegaInJet{"histoNameWeightsOmegaInJet", "", "reweighting histogram: omega in jet"};
  Configurable<std::string> histoNameWeightsOmegaUe{"histoNameWeightsOmegaUe", "", "reweighting histogram: omega in ue"};
  Configurable<std::string> histoNameWeightsAntiOmegaInJet{"histoNameWeightsAntiOmegaInJet", "", "reweighting histogram: antiomega in jet"};
  Configurable<std::string> histoNameWeightsAntiOmegaUe{"histoNameWeightsAntiOmegaUe", "", "reweighting histogram: antiomega in ue"};

  // Two-dimensional weights
  TH2F* twodWeightsPiplusJet = nullptr;
  TH2F* twodWeightsPiplusUe = nullptr;
  TH2F* twodWeightsPiminusJet = nullptr;
  TH2F* twodWeightsPiminusUe = nullptr;
  TH2F* twodWeightsK0Jet;
  TH2F* twodWeightsK0Ue;
  TH2F* twodWeightsLambdaJet;
  TH2F* twodWeightsLambdaUe;
  TH2F* twodWeightsAntilambdaJet;
  TH2F* twodWeightsAntilambdaUe;
  TH2F* twodWeightsXiInJet;
  TH2F* twodWeightsXiUe;
  TH2F* twodWeightsAntiXiInJet;
  TH2F* twodWeightsAntiXiUe;
  TH2F* twodWeightsOmegaInJet;
  TH2F* twodWeightsOmegaUe;
  TH2F* twodWeightsAntiOmegaInJet;
  TH2F* twodWeightsAntiOmegaUe;

  // List of Particles
  enum Option { KZeroLambda,
                CascadePart,
                ChargedPions,
                ChargedKaon,
                ProtonAntiproton };

  // Jet background subtraction
  JetBkgSubUtils backgroundSub;

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fOmega");
      zorro.populateHistRegistry(registryData, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    ccdb->setURL(urlToCcdb.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (applyReweighting) {
      getReweightingHistograms(ccdb);
    } else {
      twodWeightsK0Jet = nullptr;
      twodWeightsK0Ue = nullptr;
      twodWeightsLambdaJet = nullptr;
      twodWeightsLambdaUe = nullptr;
      twodWeightsAntilambdaJet = nullptr;
      twodWeightsAntilambdaUe = nullptr;
      twodWeightsXiInJet = nullptr;
      twodWeightsXiUe = nullptr;
      twodWeightsAntiXiInJet = nullptr;
      twodWeightsAntiXiUe = nullptr;
      twodWeightsOmegaInJet = nullptr;
      twodWeightsOmegaUe = nullptr;
      twodWeightsAntiOmegaInJet = nullptr;
      twodWeightsAntiOmegaUe = nullptr;
    }

    // Event Counters
    if (doprocessData) {
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
      registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {{101, 0, 101, "Multiplicity percentile"}});
    }

    if (doprocessMCefficiency) {
      registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});
    }

    // QC Histograms
    registryQC.add("deltaEtadeltaPhiJet", "deltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, o2::constants::math::PIHalf, "#Delta#phi"}});
    registryQC.add("deltaEtadeltaPhi_ue", "deltaEtadeltaPhi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, o2::constants::math::PIHalf, "#Delta#phi"}});
    registryQC.add("NchJetPlusUE", "NchJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchJet", "NchJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("NchUE", "NchUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQC.add("sumPtJetPlusUE", "sumPtJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtJet", "sumPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("sumPtUE", "sumPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQC.add("nJets_found", "nJets_found", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    registryQC.add("nJets_selected", "nJets_selected", HistType::kTH1F, {{10, 0, 10, "#it{n}_{Jet}"}});
    // registryQC.add("dcaxy_vs_pt", "dcaxy_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    // registryQC.add("dcaz_vs_pt", "dcaz_vs_pt", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    registryQC.add("jet_jet_overlaps", "jet_jet_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("jet_ue_overlaps", "jet_ue_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("ue_ue_overlaps", "ue_ue_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("tot_overlaps", "tot_overlaps", HistType::kTH2F, {{20, 0.0, 20.0, "#it{n}_{jet}"}, {200, 0.0, 200.0, "#it{n}_{overlaps}"}});
    registryQC.add("survivedK0", "survivedK0", HistType::kTH1F, {{20, 0, 20, "step"}});

    // Multiplicity Binning
    std::vector<double> multBinning = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec multAxis = {multBinning, "FT0C percentile"};
    const AxisSpec etaAxis{18, -0.9, 0.9, "#eta"};
    const AxisSpec ptAxisPi{binsOptions.binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nsigmaTOFAxis{200, -10, 10, "n#sigma_{TOF}"};
    const AxisSpec nsigmaTPCAxis{200, -10, 10, "n#sigma_{TPC}"};
    const AxisSpec dcaAxis{binsOptions.binsDca, "DCA_{xy} (cm)"};
    const AxisSpec invMassK0sAxis{200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassLambdaAxis{200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassXiAxis{200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassOmegaAxis{200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"};

    // Histograms for pions (data)
    if (doprocessData || doprocessK0s) {
      switch (particleOfInterest) {
        case KZeroLambda:
          // Histograms for lambda (data)
          registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          // Histograms for K0s (data)
          registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          break;
        case CascadePart:
          // Histograms for xi (data)
          registryData.add("XiPos_in_jet", "XiPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiPos_in_ue", "XiPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_jet", "XiNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_ue", "XiNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          // Histograms for omega (data)
          registryData.add("OmegaPos_in_jet", "OmegaPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaPos_in_ue", "OmegaPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_jet", "OmegaNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_ue", "OmegaNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          break;
        case ChargedPions:
          registryData.add("piplus_tpc_in_jet", "piplus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_jet", "piplus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_tpc_in_ue", "piplus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_ue", "piplus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_jet", "piminus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_jet", "piminus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_ue", "piminus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_ue", "piminus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_dcaxy_in_jet", "piplus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piplus_dcaxy_in_ue", "piplus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_jet", "piminus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_ue", "piminus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          break;
        case ChargedKaon: // Here we keep the name but just change the title to avoid rewriting everything
          registryData.add("piplus_tpc_in_jet", "kaplus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_jet", "kaplus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_tpc_in_ue", "kaplus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_ue", "kaplus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_jet", "kaminus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_jet", "kaminus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_ue", "kaminus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_ue", "kaminus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_dcaxy_in_jet", "kaplus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piplus_dcaxy_in_ue", "kaplus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_jet", "kaminus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_ue", "kaminus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          break;
        case ProtonAntiproton: // Here we keep the name but just change the title to avoid rewriting everything
          registryData.add("piplus_tpc_in_jet", "prplus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_jet", "prplus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_tpc_in_ue", "prplus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piplus_tof_in_ue", "prplus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_jet", "prminus_tpc_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_jet", "prminus_tof_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piminus_tpc_in_ue", "prminus_tpc_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTPCAxis});
          registryData.add("piminus_tof_in_ue", "prminus_tof_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, nsigmaTOFAxis});
          registryData.add("piplus_dcaxy_in_jet", "prplus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piplus_dcaxy_in_ue", "prplus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_jet", "prminus_dcaxy_in_jet", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          registryData.add("piminus_dcaxy_in_ue", "prminus_dcaxy_in_ue", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
          break;
        default:
          LOG(fatal) << "Cannot interpret particle " << particleOfInterest;
          break;
      }
    }

    // Histograms for efficiency (generated)
    if (doprocessMCefficiency) {
      registryMC.add("K0s_generated_jet", "K0s_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("K0s_generated_ue", "K0s_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_generated_jet", "Lambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_generated_ue", "Lambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_generated_jet", "AntiLambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_generated_ue", "AntiLambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_generated_jet", "XiPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_generated_ue", "XiPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_generated_jet", "XiNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_generated_ue", "XiNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_generated_jet", "OmegaPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_generated_ue", "OmegaPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_generated_jet", "OmegaNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_generated_ue", "OmegaNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});

      // Histograms for efficiency (reconstructed)
      registryMC.add("K0s_reconstructed_jet", "K0s_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("K0s_reconstructed_ue", "K0s_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_jet", "Lambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_ue", "Lambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_jet", "AntiLambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_ue", "AntiLambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_reconstructed_jet", "XiPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_reconstructed_ue", "XiPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_reconstructed_jet", "XiNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_reconstructed_ue", "XiNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_reconstructed_jet", "OmegaPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_reconstructed_ue", "OmegaPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_reconstructed_jet", "OmegaNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_reconstructed_ue", "OmegaNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});

      // Histograms for secondary hadrons
      registryMC.add("K0s_reconstructed_incl", "K0s_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_incl", "Lambda_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_incl", "AntiLambda_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});

      registryMC.add("XiPos_reconstructed_incl", "XiPos_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_reconstructed_incl", "XiNeg_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_reconstructed_incl", "OmegaPos_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_reconstructed_incl", "OmegaNeg_reconstructed_incl", HistType::kTH2F, {multBinning, ptAxis});

      // Histograms for secondary lambda in jet and UE
      // registryMC.add("Secondary_Lambda_InJet", "Secondary_Lambda_InJet", HistType::kTH1F, {ptAxis});
      // registryMC.add("Secondary_Lambda_Ue", "Secondary_Lambda_Ue", HistType::kTH1F, {ptAxis});
      // registryMC.add("Secondary_AntiLambda_InJet", "Secondary_AntiLambda_InJet", HistType::kTH1F, {ptAxis});
      // registryMC.add("Secondary_AntiLambda_Ue", "Secondary_AntiLambda_Ue", HistType::kTH1F, {ptAxis});

      // Histograms for 2d reweighting (pion)
      registryMC.add("mc_pi_plus_eta_pt/jet", "", HistType::kTH2F, {ptAxisPi, etaAxis});
      registryMC.addClone("mc_pi_plus_eta_pt/jet", "mc_pi_plus_eta_pt/ue");
      registryMC.addClone("mc_pi_plus_eta_pt/jet", "mc_pi_plus_eta_pt/pythia");

      registryMC.addClone("mc_pi_plus_eta_pt/", "mc_pi_minus_eta_pt/");

      registryMC.addClone("mc_pi_plus_eta_pt/", "mc_ka_plus_eta_pt/");
      registryMC.addClone("mc_pi_minus_eta_pt/", "mc_ka_minus_eta_pt/");

      registryMC.addClone("mc_pi_plus_eta_pt/", "mc_pr_plus_eta_pt/");
      registryMC.addClone("mc_pi_minus_eta_pt/", "mc_pr_minus_eta_pt/");

      // Histograms for 2d reweighting (K0s)
      registryMC.add("K0s_eta_pt_jet", "K0s_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("K0s_eta_pt_ue", "K0s_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("K0s_eta_pt_pythia", "K0s_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});

      // Histograms for 2d reweighting (Lambda)
      registryMC.add("Lambda_eta_pt_jet", "Lambda_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("Lambda_eta_pt_ue", "Lambda_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("Lambda_eta_pt_pythia", "Lambda_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("AntiLambda_eta_pt_jet", "AntiLambda_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("AntiLambda_eta_pt_ue", "AntiLambda_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("AntiLambda_eta_pt_pythia", "AntiLambda_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});

      // Histograms for 2d reweighting (Xi)
      registryMC.add("XiNeg_eta_pt_jet", "Xi_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("XiNeg_eta_pt_ue", "Xi_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("XiNeg_eta_pt_pythia", "Xi_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("XiPos_eta_pt_jet", "AntiXi_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("XiPos_eta_pt_ue", "AntiXi_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("XiPos_eta_pt_pythia", "AntiXi_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});

      // Histograms for 2d reweighting (Omega)
      registryMC.add("OmegaNeg_eta_pt_jet", "Omega_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("OmegaNeg_eta_pt_ue", "Omega_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("OmegaNeg_eta_pt_pythia", "Omega_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("OmegaPos_eta_pt_jet", "AntiOmega_eta_pt_jet", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("OmegaPos_eta_pt_ue", "AntiOmega_eta_pt_ue", HistType::kTH2F, {ptAxis, etaAxis});
      registryMC.add("OmegaPos_eta_pt_pythia", "AntiOmega_eta_pt_pythia", HistType::kTH2F, {ptAxis, etaAxis});

      // Histograms for efficiency (pions)
      registryMC.add("mc_pi_plus/in_jet/gen", "", HistType::kTH2F, {multBinning, ptAxisPi});
      registryMC.addClone("mc_pi_plus/in_jet/gen", "mc_pi_plus/in_jet/rec_tpc");
      registryMC.addClone("mc_pi_plus/in_jet/gen", "mc_pi_plus/in_jet/rec_tof");
      registryMC.addClone("mc_pi_plus/in_jet/", "mc_pi_plus/in_ue/");

      registryMC.addClone("mc_pi_plus/", "mc_pi_minus/");

      // Histograms for efficiency (kaons)
      registryMC.addClone("mc_pi_plus/", "mc_ka_plus/");
      registryMC.addClone("mc_pi_minus/", "mc_ka_minus/");

      // Histograms for efficiency (protons)
      registryMC.addClone("mc_pi_plus/", "mc_pr_plus/");
      registryMC.addClone("mc_pi_minus/", "mc_pr_minus/");

      // MC Templates
      registryMC.add("pi_plus_dcaxy/prm", "", HistType::kTH3F, {multBinning, ptAxisPi, dcaAxis});
      registryMC.addClone("pi_plus_dcaxy/prm", "pi_plus_dcaxy/sec");

      registryMC.addClone("pi_plus_dcaxy/", "pi_minus_dcaxy/");

      registryMC.addClone("pi_plus_dcaxy/", "kaplus_dcaxy/");
      registryMC.addClone("pi_minus_dcaxy/", "kaminus_dcaxy/");

      registryMC.addClone("pi_plus_dcaxy/", "prplus_dcaxy/");
      registryMC.addClone("pi_minus_dcaxy/", "prminus_dcaxy/");
    }
  }

  void getPerpendicularAxis(TVector3 p, TVector3& u, double sign)
  {
    // initialization
    double ux(0), uy(0), uz(0);

    // components of vector p
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

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
    double a = px * px + py * py;
    double b = 2.0 * px * pz * pz;
    double c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    double delta = b * b - 4.0 * a * c;

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

  // ITS hit
  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-Track Selection for Particles inside Jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    const int minTpcCr = 70;
    const double minCrFindable = 0.8;
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
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcCr)
      return false;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < minCrFindable)
      return false;
    if (track.tpcChi2NCl() > maxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > maxChi2Its)
      return false;
    if (track.eta() < -maxPseudorapidity || track.eta() > maxPseudorapidity)
      return false;
    if (track.pt() < minPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (dcaxyMaxTrackPar0 + dcaxyMaxTrackPar1 / std::pow(track.pt(), dcaxyMaxTrackPar2)))
      return false;
    if (std::fabs(track.dcaZ()) > dcazMaxTrack)
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
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < minTpcNcrossedRowsOverFindable)
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

    if (proton.Pt() < ptMinV0Proton)
      return false;
    if (proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion)
      return false;
    if (pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 Selections
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

    if (proton.Pt() < ptMinV0Proton)
      return false;
    if (proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion)
      return false;
    if (pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 Selections
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
    TVector3 pionPos(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pionNeg(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (pionPos.Pt() < ptMinK0Pion)
      return false;
    if (pionPos.Pt() > ptMaxK0Pion)
      return false;
    if (pionNeg.Pt() < ptMinK0Pion)
      return false;
    if (pionNeg.Pt() > ptMaxK0Pion)
      return false;

    // V0 Selections
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
      if (ntrack.pt() < ptMinV0Proton)
        return false;
      if (ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion)
        return false;
      if (ptrack.pt() > ptMaxV0Pion)
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
      if (ptrack.pt() < ptMinV0Proton)
        return false;
      if (ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion)
        return false;
      if (ntrack.pt() > ptMaxV0Pion)
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
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade Selections
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
      if (ntrack.pt() < ptMinV0Proton)
        return false;
      if (ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion)
        return false;
      if (ptrack.pt() > ptMaxV0Pion)
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
      if (ptrack.pt() < ptMinV0Proton)
        return false;
      if (ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion)
        return false;
      if (ntrack.pt() > ptMaxV0Pion)
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
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade Selections
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
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < minTpcNcrossedRowsOverFindable)
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
  bool isHighPurityPion(const pionTrack& track, const float nsigmaTPC, const float nsigmaTOF)
  {
    double nsigmaPi = 2.0;
    double pThreshold = 0.6;
    if (track.p() < pThreshold && std::fabs(nsigmaTPC) < nsigmaPi)
      return true;
    if (track.p() > pThreshold && std::fabs(nsigmaTPC) < nsigmaPi && std::fabs(nsigmaTOF) < nsigmaPi)
      return true;
    return false;
  }

  double getCorrectedPt(double ptRec)
  {
    // to be developed
    return ptRec;
  }

  void getReweightingHistograms(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj)
  {
    auto getWeightHistoObj = [&](Configurable<std::string> name, TH2F*& histo) {
      if (name.value == "") {
        LOG(info) << "Getting weight histogram for " << name.name << " from " << name.value;
        histo = ccdbObj->get<TH2F>(name);
      }
    };

    getWeightHistoObj(histoNameWeightPiplusJet, twodWeightsPiplusJet);
    getWeightHistoObj(histoNameWeightPiplusUe, twodWeightsPiplusUe);
    getWeightHistoObj(histoNameWeightPiminusJet, twodWeightsPiminusJet);
    getWeightHistoObj(histoNameWeightPiminusUe, twodWeightsPiminusUe);

    TList* l = ccdbObj->get<TList>(pathToFile.value.c_str());
    if (!l) {
      LOG(error) << "Could not open the file " << pathToFile.value;
      return;
    }
    l->ls();

    auto get2DWeightHisto = [&](Configurable<std::string> name, TH2F*& histo) {
      LOG(info) << "Looking for 2D weight histogram '" << name.value << "' for " << name.name;
      if (name.value == "") {
        LOG(info) << " -> Skipping";
        return;
      }
      histo = static_cast<TH2F*>(l->FindObject(name.value.c_str()));
      if (!histo) {
        LOG(error) << "Could not open histogram '" << name.value << "'";
        return;
      }
      LOG(info) << "Opened histogram " << histo->ClassName() << " " << histo->GetName();
    };

    get2DWeightHisto(histoNameWeightK0Jet, twodWeightsK0Jet);
    get2DWeightHisto(histoNameWeightK0Ue, twodWeightsK0Ue);
    get2DWeightHisto(histoNameWeightLambdaJet, twodWeightsLambdaJet);
    get2DWeightHisto(histoNameWeightLambdaUe, twodWeightsLambdaUe);
    get2DWeightHisto(histoNameWeightAntilambdaJet, twodWeightsAntilambdaJet);
    get2DWeightHisto(histoNameWeightAntilambdaUe, twodWeightsAntilambdaUe);

    // Xi
    get2DWeightHisto(histoNameWeightsXiInJet, twodWeightsXiInJet);
    get2DWeightHisto(histoNameWeightsXiUe, twodWeightsXiUe);
    get2DWeightHisto(histoNameWeightsAntiXiInJet, twodWeightsAntiXiInJet);
    get2DWeightHisto(histoNameWeightsAntiXiUe, twodWeightsAntiXiUe);

    // Omega
    get2DWeightHisto(histoNameWeightsOmegaInJet, twodWeightsOmegaInJet);
    get2DWeightHisto(histoNameWeightsOmegaUe, twodWeightsOmegaUe);
    get2DWeightHisto(histoNameWeightsAntiOmegaInJet, twodWeightsAntiOmegaInJet);
    get2DWeightHisto(histoNameWeightsAntiOmegaUe, twodWeightsAntiOmegaUe);
  }

  void processData(SelCollisions::iterator const& collision,
                   aod::V0Datas const& fullV0s,
                   aod::CascDataExt const& Cascades,
                   StrHadronDaughterTracks const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    // event counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) {
      return;
    }

    // event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // event counter: after event selection
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {

      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }

    // reject empty events
    if (fjParticles.size() < 1)
      return;
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

    // jet selection
    bool isAtLeastOneJetSelected = false;
    std::vector<TVector3> selectedJet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    for (const auto& jet : jets) { // o2-linter: disable=const-ref-in-for-loop (required by backgroundSub)

      // jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
        continue;

      // jet pt must be larger than threshold
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
      if (getCorrectedPt(jetMinusBkg.pt()) < minJetPt)
        continue;
      isAtLeastOneJetSelected = true;

      // perpendicular cone
      // double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jetAxis, ueAxis1, +1);
      getPerpendicularAxis(jetAxis, ueAxis2, -1);

      selectedJet.emplace_back(jetAxis);
      ue1.emplace_back(ueAxis1);
      ue2.emplace_back(ueAxis2);
    }
    if (!isAtLeastOneJetSelected)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Event multiplicity
    const float multiplicity = collision.centFT0M();

    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);

    for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {

      // KZeroLambda
      if (particleOfInterest == Option::KZeroLambda) {
        for (const auto& v0 : fullV0s) {

          const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
          const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();
          TVector3 v0dir(v0.px(), v0.py(), v0.pz());

          float deltaEtaJet = v0dir.Eta() - selectedJet[i].Eta();
          float deltaPhiJet = getDeltaPhi(v0dir.Phi(), selectedJet[i].Phi());
          float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          float deltaEtaUe1 = v0dir.Eta() - ue1[i].Eta();
          float deltaPhiUe1 = getDeltaPhi(v0dir.Phi(), ue1[i].Phi());
          float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          float deltaEtaUe2 = v0dir.Eta() - ue2[i].Eta();
          float deltaPhiUe2 = getDeltaPhi(v0dir.Phi(), ue2[i].Phi());
          float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

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
      }

      // Cascades
      if (particleOfInterest == Option::CascadePart) {
        for (const auto& casc : Cascades) {

          auto bach = casc.bachelor_as<StrHadronDaughterTracks>();
          auto pos = casc.posTrack_as<StrHadronDaughterTracks>();
          auto neg = casc.negTrack_as<StrHadronDaughterTracks>();

          TVector3 cascadeDir(casc.px(), casc.py(), casc.pz());
          float deltaEtaJet = cascadeDir.Eta() - selectedJet[i].Eta();
          float deltaPhiJet = getDeltaPhi(cascadeDir.Phi(), selectedJet[i].Phi());
          float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          float deltaEtaUe1 = cascadeDir.Eta() - ue1[i].Eta();
          float deltaPhiUe1 = getDeltaPhi(cascadeDir.Phi(), ue1[i].Phi());
          float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          float deltaEtaUe2 = cascadeDir.Eta() - ue2[i].Eta();
          float deltaPhiUe2 = getDeltaPhi(cascadeDir.Phi(), ue2[i].Phi());
          float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Xi+
          if (passedXiSelection(casc, pos, neg, bach, collision) &&
              bach.sign() > 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("XiPos_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("XiPos_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Xi-
          if (passedXiSelection(casc, pos, neg, bach, collision) &&
              bach.sign() < 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("XiNeg_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("XiNeg_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Omega+
          if (passedOmegaSelection(casc, pos, neg, bach, collision) &&
              bach.sign() > 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("OmegaPos_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("OmegaPos_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
          // Omega-
          if (passedOmegaSelection(casc, pos, neg, bach, collision) &&
              bach.sign() < 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("OmegaNeg_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("OmegaNeg_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
        }
      }

      // Pions
      if (particleOfInterest == Option::ChargedPions) {
        for (const auto& track : tracks) {

          if (!passedTrackSelectionForPions(track))
            continue;

          TVector3 trackDir(track.px(), track.py(), track.pz());
          const float deltaEtaJet = trackDir.Eta() - selectedJet[i].Eta();
          const float deltaPhiJet = getDeltaPhi(trackDir.Phi(), selectedJet[i].Phi());
          const float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          const float deltaEtaUe1 = trackDir.Eta() - ue1[i].Eta();
          const float deltaPhiUe1 = getDeltaPhi(trackDir.Phi(), ue1[i].Phi());
          const float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          const float deltaEtaUe2 = trackDir.Eta() - ue2[i].Eta();
          const float deltaPhiUe2 = getDeltaPhi(trackDir.Phi(), ue2[i].Phi());
          const float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          float nsigmaTPC = 999.f;
          float nsigmaTOF = 999.f;
          switch (particleOfInterest) {
            case Option::ChargedPions:
              nsigmaTPC = track.tpcNSigmaPi();
              nsigmaTOF = track.tofNSigmaPi();
              break;
            case Option::ChargedKaon:
              nsigmaTPC = track.tpcNSigmaKa();
              nsigmaTOF = track.tofNSigmaKa();
              break;
            case Option::ProtonAntiproton:
              nsigmaTPC = track.tpcNSigmaPr();
              nsigmaTOF = track.tofNSigmaPr();
              break;
          }

          bool isInJet = false;
          bool isUe = false;
          if (deltaRjet < rJet)
            isInJet = true;
          if (deltaRue1 < rJet || deltaRue2 < rJet)
            isUe = true;

          if (isHighPurityPion(track, nsigmaTPC, nsigmaTOF)) {
            if (track.sign() > 0) {
              if (isInJet)
                registryData.fill(HIST("piplus_dcaxy_in_jet"), multiplicity, track.pt(), track.dcaXY());
              if (isUe)
                registryData.fill(HIST("piplus_dcaxy_in_ue"), multiplicity, track.pt(), track.dcaXY());
            } else {
              if (isInJet)
                registryData.fill(HIST("piminus_dcaxy_in_jet"), multiplicity, track.pt(), track.dcaXY());
              if (isUe)
                registryData.fill(HIST("piminus_dcaxy_in_ue"), multiplicity, track.pt(), track.dcaXY());
            }
          }

          // DCAxy Selection
          if (std::fabs(track.dcaXY()) > dcaxyMax)
            continue;

          // TPC
          switch (track.sign()) {
            case 1:
              if (isInJet) {
                registryData.fill(HIST("piplus_tpc_in_jet"), multiplicity, track.pt(), nsigmaTPC);
              }
              if (isUe) {
                registryData.fill(HIST("piplus_tpc_in_ue"), multiplicity, track.pt(), nsigmaTPC);
              }
              break;
            case -1:
              if (isInJet) {
                registryData.fill(HIST("piminus_tpc_in_jet"), multiplicity, track.pt(), nsigmaTPC);
              }
              if (isUe) {
                registryData.fill(HIST("piminus_tpc_in_ue"), multiplicity, track.pt(), nsigmaTPC);
              }
              break;
            default:
              LOG(fatal) << "Error in the charge";
          }
          if (nsigmaTPC < nsigmaTPCmin || nsigmaTPC > nsigmaTPCmax)
            continue;
          if (!track.hasTOF())
            continue;

          // TOF
          switch (track.sign()) {
            case 1:
              if (isInJet) {
                registryData.fill(HIST("piplus_tof_in_jet"), multiplicity, track.pt(), nsigmaTOF);
              }
              if (isUe) {
                registryData.fill(HIST("piplus_tof_in_ue"), multiplicity, track.pt(), nsigmaTOF);
              }
              break;
            case -1:
              if (isInJet) {
                registryData.fill(HIST("piminus_tof_in_jet"), multiplicity, track.pt(), nsigmaTOF);
              }
              if (isUe) {
                registryData.fill(HIST("piminus_tof_in_ue"), multiplicity, track.pt(), nsigmaTOF);
              }
              break;
            default:
              LOG(fatal) << "Error in the charge";
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processData, "Process data", true);

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
      TVector3 pionPos(v0.pxpos(), v0.pypos(), v0.pzpos());
      TVector3 pionNeg(v0.pxneg(), v0.pyneg(), v0.pzneg());

      if (pionPos.Pt() < ptMinK0Pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 3.5);

      if (pionPos.Pt() > ptMaxK0Pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 4.5);

      if (pionNeg.Pt() < ptMinK0Pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 5.5);

      if (pionNeg.Pt() > ptMaxK0Pion)
        continue;
      registryQC.fill(HIST("survivedK0"), 6.5);

      // V0 Selections
      if (v0.v0cosPA() < v0cospaMin)
        continue;
      registryQC.fill(HIST("survivedK0"), 7.5);

      if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
        continue;
      registryQC.fill(HIST("survivedK0"), 8.5);

      if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
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
  PROCESS_SWITCH(StrangenessInJets, processK0s, "Process K0s", false);

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
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      registryMC.fill(HIST("number_of_events_mc"), 2.5);
      float multiplicity = collision.centFT0M();

      auto v0sPerColl = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto cascPerColl = Cascades.sliceBy(perCollisionCasc, collision.globalIndex());
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, collision.globalIndex());
      auto tracksPerColl = mcTracks.sliceBy(perCollisionTrk, collision.globalIndex());

      for (const auto& v0 : v0sPerColl) {

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

        // Generated Momentum of V0
        TVector3 momentumPos(posParticle.px(), posParticle.py(), posParticle.pz());
        TVector3 momentumNeg(negParticle.px(), negParticle.py(), negParticle.pz());
        TVector3 momentumV0 = momentumPos + momentumNeg;

        // Feed-down for lambda
        // if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0) {
        //  if (!isPhysPrim) {
        //    double wSecLambdaInJet(1.0);
        //    double wSecLambdaUe(1.0);
        //    int idMother = posParticle.mothersIds()[0];
        //    const auto& mother = mcParticles.iteratorAt(idMother);
        //    int idGrandMother = mother.mothersIds()[0];
        //    const auto& grandMother = mcParticles.iteratorAt(idGrandMother);
        //    switch (grandMother.pdgCode()) {
        //      case kXiMinus:
        //      case kXiPlusBar:
        //        // case 3322:
        //        // case -3322:
        //        if (weightsXiInJet) {
        //          int ibinXiInJet = weightsXiInJet->GetXaxis()->FindBin(grandMother.pt());
        //          wSecLambdaInJet = weightsXiInJet->GetBinContent(ibinXiInJet);
        //        }
        //        if (weightsXiUe) {
        //          int ibinXiUe = weightsXiUe->GetXaxis()->FindBin(grandMother.pt());
        //          wSecLambdaUe = weightsXiUe->GetBinContent(ibinXiUe);
        //        }
        //        break;
        //      default:
        //        break;
        //    }
        //    registryMC.fill(HIST("Secondary_Lambda_InJet"), v0.pt(), wSecLambdaInJet);
        //    registryMC.fill(HIST("Secondary_Lambda_Ue"), v0.pt(), wSecLambdaUe);
        //  }
        //}

        // Feed-down for antilambda
        // if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar) {
        //  if (!isPhysPrim) {
        //    double wSecAntiLambdaInJet(1.0);
        //    double wSecAntiLambdaUe(1.0);
        //    int idMother = posParticle.mothersIds()[0];
        //    const auto& mother = mcParticles.iteratorAt(idMother);
        //    int idGrandMother = mother.mothersIds()[0];
        //    const auto& grandMother = mcParticles.iteratorAt(idGrandMother);
        //    switch (grandMother.pdgCode()) {
        //      case kXiMinus:
        //      case kXiPlusBar:
        //        // case 3322:
        //        // case -3322:
        //        if (weightsAntiXiInJet) {
        //          int ibinAntiXiInJet = weightsAntiXiInJet->GetXaxis()->FindBin(grandMother.pt());
        //          wSecAntiLambdaInJet = weightsAntiXiInJet->GetBinContent(ibinAntiXiInJet);
        //        }
        //        if (weightsAntiXiUe) {
        //          int ibinAntiXiUe = weightsAntiXiUe->GetXaxis()->FindBin(grandMother.pt());
        //          wSecAntiLambdaUe = weightsAntiXiUe->GetBinContent(ibinAntiXiUe);
        //        }
        //        break;
        //      default:
        //        break;
        //    }
        //    registryMC.fill(HIST("Secondary_AntiLambda_InJet"), v0.pt(), wSecAntiLambdaInJet);
        //    registryMC.fill(HIST("Secondary_AntiLambda_Ue"), v0.pt(), wSecAntiLambdaUe);
        //  }
        //}

        if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short) {
          registryMC.fill(HIST("K0s_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0) {
          registryMC.fill(HIST("Lambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar) {
          registryMC.fill(HIST("AntiLambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (!isPhysPrim)
          continue;

        double wK0jet(1.0), wK0Ue(1.0), wLambdaJet(1.0), wLambdaUe(1.0), wAntilambdaJet(1.0), wAntilambdaUe(1.0);
        if (applyReweighting) {
          int ix = twodWeightsK0Jet->GetXaxis()->FindBin(momentumV0.Pt());
          int iy = twodWeightsK0Jet->GetYaxis()->FindBin(momentumV0.Eta());
          wK0jet = twodWeightsK0Jet->GetBinContent(ix, iy);
          wK0Ue = twodWeightsK0Ue->GetBinContent(ix, iy);
          wLambdaJet = twodWeightsLambdaJet->GetBinContent(ix, iy);
          wLambdaUe = twodWeightsLambdaUe->GetBinContent(ix, iy);
          wAntilambdaJet = twodWeightsAntilambdaJet->GetBinContent(ix, iy);
          wAntilambdaUe = twodWeightsAntilambdaUe->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twodWeightsK0Jet->GetNbinsX()) {
            wK0jet = 1.0;
            wK0Ue = 1.0;
            wLambdaJet = 1.0;
            wLambdaUe = 1.0;
            wAntilambdaJet = 1.0;
            wAntilambdaUe = 1.0;
          }
          if (iy == 0 || iy > twodWeightsK0Jet->GetNbinsY()) {
            wK0jet = 1.0;
            wK0Ue = 1.0;
            wLambdaJet = 1.0;
            wLambdaUe = 1.0;
            wAntilambdaJet = 1.0;
            wAntilambdaUe = 1.0;
          }
        }

        if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short) {
          registryMC.fill(HIST("K0s_reconstructed_jet"), multiplicity, v0.pt(), wK0jet);
          registryMC.fill(HIST("K0s_reconstructed_ue"), multiplicity, v0.pt(), wK0Ue);
        }
        if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0) {
          registryMC.fill(HIST("Lambda_reconstructed_jet"), multiplicity, v0.pt(), wLambdaJet);
          registryMC.fill(HIST("Lambda_reconstructed_ue"), multiplicity, v0.pt(), wLambdaUe);
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar) {
          registryMC.fill(HIST("AntiLambda_reconstructed_jet"), multiplicity, v0.pt(), wAntilambdaJet);
          registryMC.fill(HIST("AntiLambda_reconstructed_ue"), multiplicity, v0.pt(), wAntilambdaUe);
        }
      }

      // Cascades
      for (const auto& casc : cascPerColl) {
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

        int pdgParent(0);
        bool isPhysPrim = false;
        for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
          for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
            for (const auto& particleMotherOfBach : bachParticle.mothers_as<aod::McParticles>()) {
              if (particleMotherOfNeg != particleMotherOfPos)
                continue;
              if (std::abs(particleMotherOfNeg.pdgCode()) != kLambda0)
                continue;
              if (!particleMotherOfBach.isPhysicalPrimary())
                continue;

              pdgParent = particleMotherOfBach.pdgCode();
              isPhysPrim = particleMotherOfBach.isPhysicalPrimary();
            }
          }
        }
        if (pdgParent == 0)
          continue;

        // Generated Momentum of casc
        TVector3 momentumPos(posParticle.px(), posParticle.py(), posParticle.pz());
        TVector3 momentumNeg(negParticle.px(), negParticle.py(), negParticle.pz());
        TVector3 momentumBach(bachParticle.px(), bachParticle.py(), bachParticle.pz());
        TVector3 momentumV0 = momentumPos + momentumNeg;
        TVector3 momentumCasc = momentumV0 + momentumBach;

        // Xi+
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdgParent == kXiPlusBar) {
          registryMC.fill(HIST("XiPos_reconstructed_incl"), multiplicity, casc.pt());
        }
        // Xi-
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdgParent == kXiMinus) {
          registryMC.fill(HIST("XiNeg_reconstructed_incl"), multiplicity, casc.pt());
        }
        // Omega+
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdgParent == kOmegaPlusBar) {
          registryMC.fill(HIST("OmegaPos_reconstructed_incl"), multiplicity, casc.pt());
        }
        // Omega-
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdgParent == kOmegaMinus) {
          registryMC.fill(HIST("OmegaNeg_reconstructed_incl"), multiplicity, casc.pt());
        }

        if (!isPhysPrim)
          continue;

        double wXiPosJet(1.0), wXiPosUe(1.0), wXiNegJet(1.0), wXiNegUe(1.0), wOmegaPosJet(1.0), wOmegaPosUe(1.0), wOmegaNegJet(1.0), wOmegaNegUe(1.0);
        if (applyReweighting) {
          int ix = twodWeightsXiInJet->GetXaxis()->FindBin(momentumCasc.Pt());
          int iy = twodWeightsXiInJet->GetYaxis()->FindBin(momentumCasc.Eta());
          wXiPosJet = twodWeightsAntiXiInJet->GetBinContent(ix, iy);
          wXiPosUe = twodWeightsAntiXiUe->GetBinContent(ix, iy);
          wXiNegJet = twodWeightsXiInJet->GetBinContent(ix, iy);
          wXiNegUe = twodWeightsXiUe->GetBinContent(ix, iy);
          wOmegaPosJet = twodWeightsAntiOmegaInJet->GetBinContent(ix, iy);
          wOmegaPosUe = twodWeightsAntiOmegaUe->GetBinContent(ix, iy);
          wOmegaNegJet = twodWeightsOmegaInJet->GetBinContent(ix, iy);
          wOmegaNegUe = twodWeightsOmegaUe->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twodWeightsXiInJet->GetNbinsX()) {
            wXiPosJet = 1.0;
            wXiPosUe = 1.0;
            wXiNegJet = 1.0;
            wXiNegUe = 1.0;
            wOmegaPosJet = 1.0;
            wOmegaPosUe = 1.0;
            wOmegaNegJet = 1.0;
            wOmegaNegUe = 1.0;
          }
          if (iy == 0 || iy > twodWeightsXiInJet->GetNbinsY()) {
            wXiPosJet = 1.0;
            wXiPosUe = 1.0;
            wXiNegJet = 1.0;
            wXiNegUe = 1.0;
            wOmegaPosJet = 1.0;
            wOmegaPosUe = 1.0;
            wOmegaNegJet = 1.0;
            wOmegaNegUe = 1.0;
          }
        }

        // Xi+
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdgParent == kXiPlusBar) {
          registryMC.fill(HIST("XiPos_reconstructed_jet"), multiplicity, casc.pt(), wXiPosJet);
          registryMC.fill(HIST("XiPos_reconstructed_ue"), multiplicity, casc.pt(), wXiPosUe);
        }
        // Xi-
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdgParent == kXiMinus) {
          registryMC.fill(HIST("XiNeg_reconstructed_jet"), multiplicity, casc.pt(), wXiNegJet);
          registryMC.fill(HIST("XiNeg_reconstructed_ue"), multiplicity, casc.pt(), wXiNegUe);
        }
        // Omega+
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdgParent == kOmegaPlusBar) {
          registryMC.fill(HIST("OmegaPos_reconstructed_jet"), multiplicity, casc.pt(), wOmegaPosJet);
          registryMC.fill(HIST("OmegaPos_reconstructed_ue"), multiplicity, casc.pt(), wOmegaPosUe);
        }
        // Omega-
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdgParent == kOmegaMinus) {
          registryMC.fill(HIST("OmegaNeg_reconstructed_jet"), multiplicity, casc.pt(), wOmegaNegJet);
          registryMC.fill(HIST("OmegaNeg_reconstructed_ue"), multiplicity, casc.pt(), wOmegaNegUe);
        }
      }

      /*
      // Reconstructed Tracks
      for (const auto& track : tracksPerColl) {

        // Get MC Particle
        if (!track.has_mcParticle())
          continue;
        // Track Selection
        if (!passedTrackSelectionForPions(track))
          continue;

        const auto particle = track.mcParticle();
        switch (std::abs(particle.pdgCode())) {
          case kPiPlus:
            if (particle.isPhysicalPrimary()) {
              if (track.sign() > 0)
                registryMC.fill(HIST("pi_plus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("pi_minus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
            } else {
              if (track.sign() > 0)
                registryMC.fill(HIST("pi_plus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("pi_minus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
            }
            break;
          case kKPlus:
            if (particle.isPhysicalPrimary()) {
              if (track.sign() > 0)
                registryMC.fill(HIST("ka_plus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("ka_minus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
            } else {
              if (track.sign() > 0)
                registryMC.fill(HIST("ka_plus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("ka_minus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
            }
            break;
          case kProton:
            if (particle.isPhysicalPrimary()) {
              if (track.sign() > 0)
                registryMC.fill(HIST("pr_plus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("pr_minus_dcaxy/prm"), multiplicity, track.pt(), track.dcaXY());
            } else {
              if (track.sign() > 0)
                registryMC.fill(HIST("pr_plus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
              else
                registryMC.fill(HIST("pr_minus_dcaxy/sec"), multiplicity, track.pt(), track.dcaXY());
            }
            break;
          default:
            continue;
        }

        if (std::fabs(track.dcaXY()) > dcaxyMax)
          continue;

        if (track.tpcNSigmaPi() < nsigmaTPCmin || track.tpcNSigmaPi() > nsigmaTPCmax)
          continue;

        if (!particle.isPhysicalPrimary())
          continue;

        double wPiplusJet(1.0), wPiplusUe(1.0);
        double wPiminusJet(1.0), wPiminusUe(1.0);
        if (applyReweighting) {
          auto getWeight = [&](TH2F* histo) {
            if (!histo)
              return 1.0;
            const int bx = histo->GetXaxis()->FindBin(track.pt());
            if (bx <= 0 || bx > histo->GetNbinsX()) {
              return 1.0;
            }
            const int by = histo->GetYaxis()->FindBin(track.eta());
            if (by <= 0 || by > histo->GetNbinsX()) {
              return 1.0;
            }
            return histo->GetBinContent(bx, by);
          };
          wPiplusJet = getWeight(twodWeightsPiplusJet);
          wPiplusUe = getWeight(twodWeightsPiplusUe);
          wPiminusJet = getWeight(twodWeightsPiminusJet);
          wPiminusUe = getWeight(twodWeightsPiminusUe);
        }

        if (track.sign() > 0) {
          registryMC.fill(HIST("mc_pi_plus/in_jet/rec_tpc"), multiplicity, track.pt(), wPiplusJet);
          registryMC.fill(HIST("mc_pi_plus/in_ue/rec_tpc"), multiplicity, track.pt(), wPiplusUe);
        } else {
          registryMC.fill(HIST("pi_minus_rec_in_jet_tpc"), multiplicity, track.pt(), wPiminusJet);
          registryMC.fill(HIST("pi_minus_rec_in_ue_tpc"), multiplicity, track.pt(), wPiminusUe);
        }

        if (!track.hasTOF())
          continue;
        if (track.tofNSigmaPi() < nsigmaTOFminPionMC || track.tofNSigmaPi() > nsigmaTOFmaxPionMC)
          continue;

        if (track.sign() > 0) {
          registryMC.fill(HIST("mc_pi_plus/in_jet/rec_tof"), multiplicity, track.pt(), wPiplusJet);
          registryMC.fill(HIST("mc_pi_plus/in_ue/rec_tof"), multiplicity, track.pt(), wPiplusUe);
        } else {
          registryMC.fill(HIST("pi_minus_rec_in_jet_tof"), multiplicity, track.pt(), wPiminusJet);
          registryMC.fill(HIST("pi_minus_rec_in_ue_tof"), multiplicity, track.pt(), wPiminusUe);
        }
      }
      */
      for (const auto& mcParticle : mcParticlesPerColl) {

        if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;

        /*
        double wPiplusJet(1.0), wPiplusUe(1.0);
        double wPiminusJet(1.0), wPiminusUe(1.0);
        double wKaplusJet(1.0), wKaplusUe(1.0);
        double wKaminusJet(1.0), wKaminusUe(1.0);
        double wPrplusJet(1.0), wPrplusUe(1.0);
        double wPrminusJet(1.0), wPrminusUe(1.0);
        */
        double wK0jet(1.0), wK0Ue(1.0), wLambdaJet(1.0), wLambdaUe(1.0), wAntilambdaJet(1.0), wAntilambdaUe(1.0);
        if (applyReweighting) {
          /*
          auto getWeight = [&](TH2F* histo) {
            if (!histo) {
              return 1.0;
            }
            const int bx = histo->GetXaxis()->FindBin(mcParticle.pt());
            if (bx <= 0 || bx > histo->GetNbinsX()) {
              return 1.0;
            }
            const int by = histo->GetYaxis()->FindBin(mcParticle.eta());
            if (by <= 0 || by > histo->GetNbinsX()) {
              return 1.0;
            }
            return histo->GetBinContent(bx, by);
          };
          wPiplusJet = getWeight(twodWeightsPiplusJet);
          wPiplusUe = getWeight(twodWeightsPiplusUe);
          wPiminusJet = getWeight(twodWeightsPiminusJet);
          wPiminusUe = getWeight(twodWeightsPiminusUe);
          */
          int ix = twodWeightsK0Jet->GetXaxis()->FindBin(mcParticle.pt());
          int iy = twodWeightsK0Jet->GetYaxis()->FindBin(mcParticle.eta());
          wK0jet = twodWeightsK0Jet->GetBinContent(ix, iy);
          wK0Ue = twodWeightsK0Ue->GetBinContent(ix, iy);
          wLambdaJet = twodWeightsLambdaJet->GetBinContent(ix, iy);
          wLambdaUe = twodWeightsLambdaUe->GetBinContent(ix, iy);
          wAntilambdaJet = twodWeightsAntilambdaJet->GetBinContent(ix, iy);
          wAntilambdaUe = twodWeightsAntilambdaUe->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twodWeightsK0Jet->GetNbinsX()) {
            wK0jet = 1.0;
            wK0Ue = 1.0;
            wLambdaJet = 1.0;
            wLambdaUe = 1.0;
            wAntilambdaJet = 1.0;
            wAntilambdaUe = 1.0;
          }
          if (iy == 0 || iy > twodWeightsK0Jet->GetNbinsY()) {
            wK0jet = 1.0;
            wK0Ue = 1.0;
            wLambdaJet = 1.0;
            wLambdaUe = 1.0;
            wAntilambdaJet = 1.0;
            wAntilambdaUe = 1.0;
          }
        }

        double wXiPosJet(1.0), wXiPosUe(1.0), wXiNegJet(1.0), wXiNegUe(1.0), wOmegaPosJet(1.0), wOmegaPosUe(1.0), wOmegaNegJet(1.0), wOmegaNegUe(1.0);
        if (applyReweighting) {
          int ix = twodWeightsXiInJet->GetXaxis()->FindBin(mcParticle.pt());
          int iy = twodWeightsXiInJet->GetYaxis()->FindBin(mcParticle.eta());
          wXiPosJet = twodWeightsAntiXiInJet->GetBinContent(ix, iy);
          wXiPosUe = twodWeightsAntiXiUe->GetBinContent(ix, iy);
          wXiNegJet = twodWeightsXiInJet->GetBinContent(ix, iy);
          wXiNegUe = twodWeightsXiUe->GetBinContent(ix, iy);
          wOmegaPosJet = twodWeightsAntiOmegaInJet->GetBinContent(ix, iy);
          wOmegaPosUe = twodWeightsAntiOmegaUe->GetBinContent(ix, iy);
          wOmegaNegJet = twodWeightsOmegaInJet->GetBinContent(ix, iy);
          wOmegaNegUe = twodWeightsOmegaUe->GetBinContent(ix, iy);

          // protections
          if (ix == 0 || ix > twodWeightsXiInJet->GetNbinsX()) {
            wXiPosJet = 1.0;
            wXiPosUe = 1.0;
            wXiNegJet = 1.0;
            wXiNegUe = 1.0;
            wOmegaPosJet = 1.0;
            wOmegaPosUe = 1.0;
            wOmegaNegJet = 1.0;
            wOmegaNegUe = 1.0;
          }
          if (iy == 0 || iy > twodWeightsXiInJet->GetNbinsY()) {
            wXiPosJet = 1.0;
            wXiPosUe = 1.0;
            wXiNegJet = 1.0;
            wXiNegUe = 1.0;
            wOmegaPosJet = 1.0;
            wOmegaPosUe = 1.0;
            wOmegaNegJet = 1.0;
            wOmegaNegUe = 1.0;
          }
        }

        switch (mcParticle.pdgCode()) {
          /*
          case kPiPlus: // Pi+
            registryMC.fill(HIST("mc_pi_plus/in_jet/gen"), multiplicity, mcParticle.pt(), wPiplusJet);
            registryMC.fill(HIST("mc_pi_plus/in_ue/gen"), multiplicity, mcParticle.pt(), wPiplusUe);
            registryMC.fill(HIST("pi_plus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kPiMinus: // Pi-
            registryMC.fill(HIST("mc_pi_minus/in_jet/gen"), multiplicity, mcParticle.pt(), wPiminusJet);
            registryMC.fill(HIST("mc_pi_minus/in_ue/gen"), multiplicity, mcParticle.pt(), wPiminusUe);
            registryMC.fill(HIST("pi_minus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kKPlus: // Ka+
            registryMC.fill(HIST("mc_ka_plus/in_jet/gen"), multiplicity, mcParticle.pt(), wKaplusJet);
            registryMC.fill(HIST("mc_ka_plus/in_ue/gen"), multiplicity, mcParticle.pt(), wKaplusUe);
            registryMC.fill(HIST("ka_plus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kKMinus: // Ka-
            registryMC.fill(HIST("mc_ka_minus/in_jet/gen"), multiplicity, mcParticle.pt(), wKaminusJet);
            registryMC.fill(HIST("mc_ka_minus/in_ue/gen"), multiplicity, mcParticle.pt(), wKaminusUe);
            registryMC.fill(HIST("ka_minus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kProton: // Pr+
            registryMC.fill(HIST("mc_pr_plus/in_jet/gen"), multiplicity, mcParticle.pt(), wPrplusJet);
            registryMC.fill(HIST("mc_pr_plus/in_ue/gen"), multiplicity, mcParticle.pt(), wPrplusUe);
            registryMC.fill(HIST("pr_plus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kProtonBar: // Pr-
            registryMC.fill(HIST("mc_pr_minus/in_jet/gen"), multiplicity, mcParticle.pt(), wPrminusJet);
            registryMC.fill(HIST("mc_pr_minus/in_ue/gen"), multiplicity, mcParticle.pt(), wPrminusUe);
            registryMC.fill(HIST("pr_minus_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          */
          case kK0Short: // K0s
            registryMC.fill(HIST("K0s_generated_jet"), multiplicity, mcParticle.pt(), wK0jet);
            registryMC.fill(HIST("K0s_generated_ue"), multiplicity, mcParticle.pt(), wK0Ue);
            registryMC.fill(HIST("K0s_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kLambda0: // Lambda
            registryMC.fill(HIST("Lambda_generated_jet"), multiplicity, mcParticle.pt(), wLambdaJet);
            registryMC.fill(HIST("Lambda_generated_ue"), multiplicity, mcParticle.pt(), wLambdaUe);
            registryMC.fill(HIST("Lambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kLambda0Bar: // AntiLambda
            registryMC.fill(HIST("AntiLambda_generated_jet"), multiplicity, mcParticle.pt(), wAntilambdaJet);
            registryMC.fill(HIST("AntiLambda_generated_ue"), multiplicity, mcParticle.pt(), wAntilambdaUe);
            registryMC.fill(HIST("AntiLambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kXiPlusBar: // Xi Pos
            registryMC.fill(HIST("XiPos_generated_jet"), multiplicity, mcParticle.pt(), wXiPosJet);
            registryMC.fill(HIST("XiPos_generated_ue"), multiplicity, mcParticle.pt(), wXiPosUe);
            registryMC.fill(HIST("XiPos_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kXiMinus: // Xi Neg
            registryMC.fill(HIST("XiNeg_generated_jet"), multiplicity, mcParticle.pt(), wXiNegJet);
            registryMC.fill(HIST("XiNeg_generated_ue"), multiplicity, mcParticle.pt(), wXiNegUe);
            registryMC.fill(HIST("XiNeg_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kOmegaPlusBar: // Omega Pos
            registryMC.fill(HIST("OmegaPos_generated_jet"), multiplicity, mcParticle.pt(), wOmegaPosJet);
            registryMC.fill(HIST("OmegaPos_generated_ue"), multiplicity, mcParticle.pt(), wOmegaPosUe);
            registryMC.fill(HIST("OmegaPos_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
          case kOmegaMinus: // Omega Neg
            registryMC.fill(HIST("OmegaNeg_generated_jet"), multiplicity, mcParticle.pt(), wOmegaNegJet);
            registryMC.fill(HIST("OmegaNeg_generated_ue"), multiplicity, mcParticle.pt(), wOmegaNegUe);
            registryMC.fill(HIST("OmegaNeg_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processMCefficiency, "Process MC Efficiency", false);

  void processGen(o2::aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mcCollisions) {

      // Selection on z_{vertex}
      if (std::fabs(mccollision.posZ()) > zVtx)
        continue;

      // MC Particles per Collision
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      // loop over generated MC particles
      std::vector<fastjet::PseudoJet> fjParticles;
      for (const auto& particle : mcParticlesPerColl) {

        if (!particle.isPhysicalPrimary())
          continue;
        if (std::fabs(particle.eta()) > etaMax)
          continue;
        double ptMinPart = 0.1;
        if (particle.pt() < ptMinPart)
          continue;

        // 4-momentum representation of a particle
        double energy = std::sqrt(particle.p() * particle.p() + MassPionCharged * MassPionCharged);
        fastjet::PseudoJet fourMomentum(particle.px(), particle.py(), particle.pz(), energy);
        fjParticles.emplace_back(fourMomentum);
      }
      // reject empty events
      if (fjParticles.size() < 1)
        continue;

      // cluster particles using the anti-kt algorithm
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

      // jet selection
      bool isAtLeastOneJetSelected = false;
      std::vector<TVector3> selectedJet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      for (const auto& jet : jets) { // o2-linter: disable=const-ref-in-for-loop (required by backgroundSub)

        // jet must be fully contained in the acceptance
        if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
          continue;

        // jet pt must be larger than threshold
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
        if (getCorrectedPt(jetMinusBkg.pt()) < minJetPt)
          continue;
        isAtLeastOneJetSelected = true;

        // perpendicular cone
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        TVector3 ueAxis1(0, 0, 0);
        TVector3 ueAxis2(0, 0, 0);
        getPerpendicularAxis(jetAxis, ueAxis1, +1);
        getPerpendicularAxis(jetAxis, ueAxis2, -1);
        selectedJet.emplace_back(jetAxis);
        ue1.emplace_back(ueAxis1);
        ue2.emplace_back(ueAxis2);
      }
      if (!isAtLeastOneJetSelected)
        continue;

      // loop over selected jets / UE
      for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {
        for (const auto& particle : mcParticlesPerColl) {

          if (!particle.isPhysicalPrimary())
            continue;
          if (std::fabs(particle.eta()) > etaMax)
            continue;
          double ptMinPart = 0.1;
          if (particle.pt() < ptMinPart)
            continue;

          TVector3 particleDir(particle.px(), particle.py(), particle.pz());
          float deltaEtaJet = particleDir.Eta() - selectedJet[i].Eta();
          float deltaPhiJet = getDeltaPhi(particleDir.Phi(), selectedJet[i].Phi());
          float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          float deltaEtaUe1 = particleDir.Eta() - ue1[i].Eta();
          float deltaPhiUe1 = getDeltaPhi(particleDir.Phi(), ue1[i].Phi());
          float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          float deltaEtaUe2 = particleDir.Eta() - ue2[i].Eta();
          float deltaPhiUe2 = getDeltaPhi(particleDir.Phi(), ue2[i].Phi());
          float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // In jet
          if (deltaRjet < rJet) {
            switch (particle.pdgCode()) {
              /*
              case kPiPlus:
                registryMC.fill(HIST("pi_plus_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kPiMinus:
                registryMC.fill(HIST("pi_minus_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              */
              case kK0Short:
                registryMC.fill(HIST("K0s_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kLambda0:
                registryMC.fill(HIST("Lambda_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kLambda0Bar:
                registryMC.fill(HIST("AntiLambda_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kXiMinus:
                registryMC.fill(HIST("XiNeg_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kXiPlusBar:
                registryMC.fill(HIST("XiPos_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kOmegaMinus:
                registryMC.fill(HIST("OmegaNeg_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              case kOmegaPlusBar:
                registryMC.fill(HIST("OmegaPos_eta_pt_jet"), particle.pt(), particle.eta());
                break;
              default:
                continue;
            }
          }

          if (deltaRue1 < rJet || deltaRue2 < rJet) {
            switch (particle.pdgCode()) {
              /*
              case kPiPlus:
                registryMC.fill(HIST("pi_plus_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kPiMinus:
                registryMC.fill(HIST("pi_minus_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              */
              case kK0Short:
                registryMC.fill(HIST("K0s_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kLambda0:
                registryMC.fill(HIST("Lambda_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kLambda0Bar:
                registryMC.fill(HIST("AntiLambda_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kXiMinus:
                registryMC.fill(HIST("XiNeg_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kXiPlusBar:
                registryMC.fill(HIST("XiPos_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kOmegaMinus:
                registryMC.fill(HIST("OmegaNeg_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              case kOmegaPlusBar:
                registryMC.fill(HIST("OmegaPos_eta_pt_ue"), particle.pt(), particle.eta());
                break;
              default:
                continue;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processGen, "Process generated MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StrangenessInJets>(cfgc)};
}
