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

/// \file netchargeFluctuations.cxx
/// \brief Calculate net-charge fluctuations using nu_dyn observable
///        For charged particles
///        For RUN-3
///
/// \author Nida Malik <nida.malik@cern.ch>
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using namespace o2::constants::physics;

enum RunType {
  kRun3 = 0,
  kRun2
};

// Structure to handle net charge fluctuation analysis
struct NetchargeFluctuations {

  // Macro to define configurable parameters with default values and help text

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

  // Services for PDG and CCDB (Calibration and Condition Database)
  Service<o2::framework::O2DatabasePDG> pdgService; // Particle data group service
  Service<o2::ccdb::BasicCCDBManager> ccdb;         // CCDB manager service

  // Random number generator for statistical fluctuations, initialized with seed 0
  TRandom3* fRndm = new TRandom3(0);

  // Registry for histograms used in analysis
  HistogramRegistry histogramRegistry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // -------------------
  // Configurable parameters
  // -------------------
  // CCDB related configurations
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://alice-ccdb.cern.ch", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/n/nimalik/My/Object/pn", "Path for ccdb-object"};
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};
  Configurable<bool> cfgEffNue{"cfgEffNue", false, "efficiency correction to nu_dyn"};

  // Track and event selection cuts
  Configurable<float> vertexZcut{"vertexZcut", 10.f, "Vertex Z"};
  Configurable<float> etaCut{"etaCut", 0.8f, "Eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.2, "Pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 5.0, "Pt max cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.2, "DCA XY cut"};
  Configurable<float> dcaZCut{"dcaZCut", 2.0, "DCA Z cut"};
  Configurable<float> tpcCrossCut{"tpcCrossCut", 70., "TPC crossrows cut"};
  Configurable<float> itsChiCut{"itsChiCut", 36., "ITS chi2 cluster cut"};
  Configurable<float> tpcChiCut{"tpcChiCut", 4., "TPC chi2 cluster cut"};
  Configurable<float> centMin{"centMin", 0.0f, "cenrality min for delta eta"};
  Configurable<float> centMax{"centMax", 10.0f, "cenrality max for delta eta"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 30, "Number of subsamples for Error"};
  Configurable<int> deltaEta{"deltaEta", 8, "Delta eta bin count"};
  Configurable<double> threshold{"threshold", 1e-6, "Delta eta bin count"};

  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"}; // sel8
  Configurable<bool> cInt7Trig{"cInt7Trig", true, "kINT7 MB Trigger"};                // kINT7
  Configurable<bool> cSel7Trig{"cSel7Trig", true, "Sel7 (V0A + V0C) Selection Run2"}; // sel7
  Configurable<bool> cDcaXy{"cDcaXy", true, "Dca XY cut"};
  Configurable<bool> cDcaZ{"cDcaZ", true, "Dca Z cut"};
  Configurable<bool> cTpcCr{"cTpcCr", true, "tpc crossrows"};
  Configurable<bool> cItsChi{"cItsChi", true, "ITS chi"};
  Configurable<bool> cTpcChi{"cTpcChi", true, "TPC chi"};
  Configurable<bool> cFT0C{"cFT0C", true, "cent FT0C"};
  Configurable<bool> cFT0M{"cFT0M", false, "cent FT0M"};

  // Centrality binning configuration
  ConfigurableAxis centBining{"centBining", {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality/Multiplicity percentile bining"};
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};                          // pileup
  Configurable<bool> cfgUseGoodItsLayerAllCut{"cfgUseGoodItsLayerAllCut", false, "Good ITS Layers All"}; // pileup
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};                        // pileup
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};                     // pileup
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};                        // pileup
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};                 // pileup
  Configurable<bool> cPVcont{"cPVcont", false, "primary vertex contributor"};

  // Configurable to enable multiplicity correlation cuts
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, false, "Multiplicity correlation cut")

  // Struct grouping multiplicity vs centrality/vertex cuts and related parameters
  struct : ConfigurableGroup {

    // Flags to enable specific multiplicity correlation cuts
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, true, "Enable PV multiplicity vs T0C centrality cut")
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalFT0CCutEnabled, bool, true, "Enable globalTracks vs FT0C multiplicity cut")
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, true, "Enable globalTracks vs PV multiplicity cut")

    // Parameter values for PV multiplicity vs FT0C centrality cut (polynomial coefficients, etc.)
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars",
                                                          std::vector<double>{30.434, -0.917137, 0.0185032, -0.000198425, 7.94381e-07, 13.7406, -0.282656, 0.00556147, -6.32766e-05, 2.51648e-07},
                                                          "PV multiplicity vs T0C centrality cut parameter values"};

    // Parameter values for globalTracks vs FT0C multiplicity cut
    Configurable<std::vector<double>> cfgMultGlobalFT0CCutPars{"cfgMultGlobalFT0CCutPars",
                                                               std::vector<double>{18.9628, -0.576466, 0.0117324, -0.000126086, 5.05365e-07, 8.99921, -0.188022, 0.0037089, -4.20275e-05, 1.68234e-07},
                                                               "globalTracks vs FT0C cut parameter values"};

    // Parameter values for globalTracks vs PV multiplicity cut
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars",
                                                             std::vector<double>{0.148031, 0.616699, 0.603083, 0.112751, -0.0013846, 8.38211e-06},
                                                             "globalTracks vs PV cut parameter values"};

    // Local vectors to store the above parameters
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalFT0CPars;
    std::vector<double> multGlobalPVCutPars;

    // TF1 objects to represent low/high cut functions for the above correlations
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultGlobalFT0CCutLow = nullptr;
    TF1* fMultGlobalFT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;

  } cfgFunCoeff;

  // Histogram pointer for CCDB efficiency
  // TH1D* efficiency = nullptr;

  TH1D* efficiencyPos = nullptr;
  TH1D* efficiencyNeg = nullptr;

  // Filters for selecting collisions and tracks
  Filter collisionFilter = nabs(aod::collision::posZ) <= vertexZcut;
  Filter trackFilter = (nabs(aod::track::eta) < etaCut) && (aod::track::pt > ptMinCut) && (aod::track::pt < ptMaxCut) && (requireGlobalTrackInFilter());

  using MyCollisionsRun2 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>>;
  using MyCollisionRun2 = MyCollisionsRun2::iterator;

  using MyCollisionsRun3 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults>>;
  using MyCollisionRun3 = MyCollisionsRun3::iterator;

  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using MyTrack = MyTracks::iterator;

  using MyMCCollisionsRun2 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults, aod::McCollisionLabels>>;
  using MyMCCollisionRun2 = MyMCCollisionsRun2::iterator;

  using MyMCCollisionsRun3 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>>;
  using MyMCCollisionRun3 = MyMCCollisionsRun3::iterator;

  using MyMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::StoredTracks, aod::TrackSelection, aod::McTrackLabels>>;
  using MyMCTrack = MyMCTracks::iterator;

  void init(o2::framework::InitContext&)
  {
    // -------------------------------
    // Define histogram axes specifications
    // -------------------------------
    const AxisSpec vtxzAxis = {800, -20, 20, "V_{Z} (cm)"};
    const AxisSpec dcaAxis = {1000, -0.5, 0.5, "DCA_{xy} (cm)"};
    const AxisSpec dcazAxis = {600, -3, 3, "DCA_{z} (cm)"};
    const AxisSpec phiAxis = {70, 0, 7, "#phi "};
    const AxisSpec ptAxis = {70, 0.0, 7.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis = {20, -1., 1., "#eta"};
    const AxisSpec deltaEtaAxis = {9, 0, 1.8, "#eta"};
    const AxisSpec centAxis = {100, 0., 100., "centrality"};
    const AxisSpec multAxis = {100000, 0., 100000., "FT0M Amplitude"};
    const AxisSpec tpcChiAxis = {700, 0., 7., "Chi2"};
    const AxisSpec itsChiAxis = {400, 0., 40., "Chi2"};
    const AxisSpec crossedRowAxis = {1600, 0., 160., "TPC Crossed rows"};
    const AxisSpec eventsAxis = {10, 0, 10, ""};
    const AxisSpec signAxis = {20, -10, 10, ""};
    const AxisSpec nchAxis = {5000, 0, 5000, "Nch"};
    const AxisSpec nch1Axis = {1500, 0, 1500, "Nch"};
    const AxisSpec nchpAxis = {50000, 0, 50000, "Nch"};
    const AxisSpec cent1Axis{centBining, "Multiplicity percentile from FT0M (%)"};

    // Subsample axis (used for error estimation from subsamples)
    auto noSubsample = static_cast<int>(cfgNSubsample);
    float maxSubsample = 1.0 * noSubsample;
    AxisSpec subsampleAxis = {noSubsample, 0.0, maxSubsample, "subsample no."};

    // Add QA histograms
    histogramRegistry.add("QA/hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("QA/hDcaXY_before", "", kTH1F, {dcaAxis});
    histogramRegistry.add("QA/hphi", "", kTH1F, {phiAxis});
    histogramRegistry.add("QA/hDcaZ_before", "", kTH1F, {dcazAxis});
    histogramRegistry.add("QA/hTPCchi2perCluster_before", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("QA/hITSchi2perCluster_before", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("QA/hTPCCrossedrows_before", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("QA/hPtDcaXY_before", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("QA/hPtDcaZ_before", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("QA/hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("QA/hDcaXY_after", "", kTH1F, {dcaAxis});
    histogramRegistry.add("QA/hDcaZ_after", "", kTH1F, {dcazAxis});
    histogramRegistry.add("QA/hTPCchi2perCluster_after", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("QA/hITSchi2perCluster_after", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("QA/hTPCCrossedrows_after", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("QA/hPtDcaXY_after", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("QA/hPtDcaZ_after", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("QA/hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("QA/cent_hEta", "", kTH2F, {cent1Axis, etaAxis});
    histogramRegistry.add("QA/hPt", "", kTH1F, {ptAxis});
    histogramRegistry.add("QA/cent_hPt", "", kTH2F, {cent1Axis, ptAxis});
    histogramRegistry.add("QA/hPt_eta", "", kTH2F, {ptAxis, etaAxis});
    histogramRegistry.add("QA/hCentrality", "", kTH1F, {centAxis});
    histogramRegistry.add("QA/hMultiplicity", "", kTH1F, {multAxis});

    histogramRegistry.add("gen/hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen/hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen/hPt", "", kTH1F, {ptAxis});
    histogramRegistry.add("gen/cent_hPt", "", kTH2F, {centAxis, ptAxis});
    histogramRegistry.add("gen/hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("gen/cent_hEta", "", kTH2F, {centAxis, etaAxis});
    histogramRegistry.add("gen/hSign", "", kTH1F, {signAxis});
    histogramRegistry.add("gen/hPt_eta", "", kTH2F, {ptAxis, etaAxis});
    histogramRegistry.add("gen/cent_pos", "cent vs fpos", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_neg", "cent vs fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_termp", "cent vs termp", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_termn", "cent vs termn", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_pos_sq", "cent vs sqfpos", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_neg_sq", "cent vs sqfneg", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_posneg", "cent vs fpos*fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/cent_nch", "cent vs nch", kTProfile, {cent1Axis});
    histogramRegistry.add("gen/nch", "", kTH1F, {nchAxis});
    histogramRegistry.add("gen/delta_eta_eta", "delta_eta ", kTH1F, {etaAxis});
    histogramRegistry.add("gen/delta_eta_pos", "delta_eta vs fpos ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_neg", "delta_eta vs fneg ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_termp", "delta_eta vs termp ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_termn", "delta_eta vs termn ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_pos_sq", "delta_eta vs pos_sq ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_neg_sq", "delta_eta vs neg_sq ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_posneg", "delta_eta vs posneg ", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("gen/delta_eta_nch", "delta_eta vs nchGen ", kTProfile, {deltaEtaAxis});

    histogramRegistry.add("data/nch", "", kTH1D, {nchAxis});
    histogramRegistry.add("data/cent_nch", "", kTProfile, {cent1Axis});
    histogramRegistry.add("data/nch_pos", "", kTH1D, {nchAxis});
    histogramRegistry.add("data/cent_nch_pos", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("data/nch_neg", "", kTH1D, {nchAxis});
    histogramRegistry.add("data/cent_nch_neg", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("data/nch_negpos", "", kTH1D, {nchpAxis});
    histogramRegistry.add("data/cent_nch_negpos", "", kTH2D, {centAxis, nchpAxis});
    histogramRegistry.add("data/cent_pos", "cent vs fpos", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_neg", "cent vs fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_termp", "cent vs termp", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_termn", "cent vs termn", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_pos_sq", "cent vs sqfpos", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_neg_sq", "cent vs sqfneg", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_posneg", "cent vs fpos*fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("data/hPt_cor", "", kTH1F, {ptAxis});
    histogramRegistry.add("data/hEta_cor", "", kTH1F, {etaAxis});
    histogramRegistry.add("data/nch_nchCor", "", kTProfile, {nchAxis});
    histogramRegistry.add("data/nchCor", "", kTH1F, {nchAxis});
    histogramRegistry.add("data/cent_nchCor", "", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_pos_cor", "", kTProfile, {cent1Axis});
    histogramRegistry.add("data/cent_neg_cor", "", kTProfile, {cent1Axis});
    histogramRegistry.add("data/delta_eta_cent", "Centrality", kTH1F, {cent1Axis});
    histogramRegistry.add("data/delta_eta_eta", "eta", kTH1F, {etaAxis});
    histogramRegistry.add("data/delta_eta_nchTotal", "delta_eta vs nchTotal", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_nch", "delta_eta vs nch", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_nchCor", "delta_eta vs nchCor", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_pos", "delta_eta vs fpos", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_neg", "delta_eta vs fneg", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_termp", "delta_eta vs termp", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_termn", "delta_eta vs termn", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_pos_sq", "delta_eta vs sqfpos", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_neg_sq", "delta_eta vs sqfneg", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_posneg", "delta_eta vs fpos*fneg", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_pos_cor", "delta_eta vs fpos_cor", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("data/delta_eta_neg_cor", "delta_eta vs fneg_cor", kTProfile, {deltaEtaAxis});

    histogramRegistry.add("subsample/pos", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/neg", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/termp", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/termn", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/pos_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/neg_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/posneg", "", kTProfile2D, {cent1Axis, subsampleAxis});

    histogramRegistry.add("subsample/gen/pos", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/neg", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/termp", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/termn", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/pos_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/neg_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/gen/posneg", "", kTProfile2D, {cent1Axis, subsampleAxis});

    histogramRegistry.add("subsample/delta_eta/pos", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/neg", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/termp", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/termn", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/pos_sq", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/neg_sq", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/posneg", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});

    histogramRegistry.add("subsample/delta_eta/gen/pos", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/neg", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/termp", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/termn", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/pos_sq", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/neg_sq", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});
    histogramRegistry.add("subsample/delta_eta/gen/posneg", "", kTProfile2D, {deltaEtaAxis, subsampleAxis});

    histogramRegistry.add("QA/hCentFT0C", "", kTH1F, {centAxis});
    histogramRegistry.add("QA/hNchGlobal", "", kTH1F, {nchAxis});
    histogramRegistry.add("QA/hNchPV", "", kTH1F, {nchAxis});

    histogramRegistry.add("eff/hPt_np_gen", "", kTH1F, {ptAxis});
    histogramRegistry.add("eff/hPt_nm_gen", "", kTH1F, {ptAxis});
    histogramRegistry.add("eff/hPt_np", "", kTH1F, {ptAxis});
    histogramRegistry.add("eff/hPt_nm", "", kTH1F, {ptAxis});

    // QA histograms for multiplicity correlations
    histogramRegistry.add("MultCorrelationPlots/globalTracks_PV_bef", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histogramRegistry.add("MultCorrelationPlots/globalTracks_FT0C_bef", "", {HistType::kTH2D, {centAxis, nchAxis}});
    histogramRegistry.add("MultCorrelationPlots/PV_FT0C_bef", "", {HistType::kTH2D, {centAxis, nchAxis}});

    histogramRegistry.add("MultCorrelationPlots/globalTracks_PV_aft", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histogramRegistry.add("MultCorrelationPlots/globalTracks_FT0C_aft", "", {HistType::kTH2D, {centAxis, nchAxis}});
    histogramRegistry.add("MultCorrelationPlots/PV_FT0C_aft", "", {HistType::kTH2D, {centAxis, nchAxis}});

    cfgFunCoeff.multPVT0CCutPars = cfgFunCoeff.cfgMultPVT0CCutPars;
    cfgFunCoeff.multGlobalFT0CPars = cfgFunCoeff.cfgMultGlobalFT0CCutPars;
    cfgFunCoeff.multGlobalPVCutPars = cfgFunCoeff.cfgMultGlobalPVCutPars;

    // --- Initialize PV vs FT0C multiplicity cut functions ---
    // Lower cut function: 4th-order polynomial minus 3.5 sigma deviation

    cfgFunCoeff.fMultPVT0CCutLow =
      new TF1("fMultPVT0CCutLow",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultPVT0CCutLow->SetParameters(&(cfgFunCoeff.multPVT0CCutPars[0]));

    // Upper cut function: 4th-order polynomial plus 3.5 sigma deviation
    cfgFunCoeff.fMultPVT0CCutHigh =
      new TF1("fMultPVT0CCutHigh",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultPVT0CCutHigh->SetParameters(&(cfgFunCoeff.multPVT0CCutPars[0]));

    // --- Initialize globalTracks vs FT0C multiplicity cut functions ---
    // Lower cut function
    cfgFunCoeff.fMultGlobalFT0CCutLow =
      new TF1("fMultGlobalFT0CCutLow",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalFT0CCutLow->SetParameters(&(cfgFunCoeff.multGlobalFT0CPars[0]));

    // Upper cut function
    cfgFunCoeff.fMultGlobalFT0CCutHigh =
      new TF1("fMultGlobalFT0CCutHigh",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalFT0CCutHigh->SetParameters(&(cfgFunCoeff.multGlobalFT0CPars[0]));

    // --- Initialize globalTracks vs PV multiplicity cut functions ---
    // Lower cut: linear + cubic term minus 3.5 sigma

    cfgFunCoeff.fMultGlobalPVCutLow =
      new TF1("fMultGlobalPVCutLow",
              "[0]+[1]*x - 3.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalPVCutLow->SetParameters(&(cfgFunCoeff.multGlobalPVCutPars[0]));

    // Upper cut: linear + cubic term plus 3.5 sigma
    cfgFunCoeff.fMultGlobalPVCutHigh =
      new TF1("fMultGlobalPVCutHigh",
              "[0]+[1]*x + 3.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalPVCutHigh->SetParameters(&(cfgFunCoeff.multGlobalPVCutPars[0]));

    // --- Load efficiency histogram from CCDB
    if (cfgLoadEff) {
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      TList* list = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      efficiencyPos = reinterpret_cast<TH1D*>(list->FindObject("efficiency_Pos_Run3"));
      efficiencyNeg = reinterpret_cast<TH1D*>(list->FindObject("efficiency_Neg_Run3"));
      // Log fatal error if efficiency histogram is not found
      if (!efficiencyPos || !efficiencyNeg) {
        LOGF(info, "FATAL!! Could not find required histograms in CCDB");
      }
    }
  }

  bool eventSelected(const float& globalNch, const float& pvTrack, const float& centrality)
  {
    if (cfgFunCoeff.cfgMultPVT0CCutEnabled) {

      if (pvTrack < cfgFunCoeff.fMultPVT0CCutLow->Eval(centrality))
        return false;
      if (pvTrack > cfgFunCoeff.fMultPVT0CCutHigh->Eval(centrality))
        return false;
    }

    if (cfgFunCoeff.cfgMultGlobalFT0CCutEnabled) {

      if (globalNch < cfgFunCoeff.fMultGlobalFT0CCutLow->Eval(centrality))
        return false;
      if (globalNch > cfgFunCoeff.fMultGlobalFT0CCutHigh->Eval(centrality))
        return false;
    }

    if (cfgFunCoeff.cfgMultGlobalPVCutEnabled) {

      if (globalNch < cfgFunCoeff.fMultGlobalPVCutLow->Eval(pvTrack))
        return false;
      if (globalNch > cfgFunCoeff.fMultGlobalPVCutHigh->Eval(pvTrack))
        return false;
    }

    return true;
  }

  template <RunType run, typename C>
  bool selCollision(C const& coll, float& cent, float& mult)
  {

    if (std::abs(coll.posZ()) >= vertexZcut)
      return false;
    if constexpr (run == kRun3) {
      if (cSel8Trig && !coll.sel8()) {
        return false;
      }
      if (cFT0M) {
        cent = coll.centFT0M(); // centrality for run3 using FT0M
        mult = coll.multFT0M();
      } else if (cFT0C) {
        cent = coll.centFT0C(); // centrality for run3 using FT0C
        mult = coll.multFT0C();
      }

    } else if constexpr (run == kRun2) {
      if (cInt7Trig && !coll.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !coll.sel7()) {
        return false;
      }
      cent = coll.centRun2V0M(); // centrality for run2
      mult = coll.multFV0M();    // multiplicity for run2
    }

    if (cNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (cTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (cPileupReject && !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (cZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (cItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (cfgUseGoodItsLayerAllCut && !(coll.selection_bit(aod::evsel::kIsGoodITSLayersAll)))
      return false;

    return true;
  }

  template <typename T>
  void fillBeforeQA(T const& track)
  {
    histogramRegistry.fill(HIST("QA/hTPCchi2perCluster_before"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("QA/hITSchi2perCluster_before"), track.itsChi2NCl());
    histogramRegistry.fill(HIST("QA/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("QA/hDcaXY_before"), track.dcaXY());
    histogramRegistry.fill(HIST("QA/hDcaZ_before"), track.dcaZ());
    histogramRegistry.fill(HIST("QA/hPtDcaXY_before"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("QA/hPtDcaZ_before"), track.pt(), track.dcaZ());
  }

  template <typename T>
  void fillAfterQA(T const& track)
  {
    histogramRegistry.fill(HIST("QA/hphi"), track.phi());
    histogramRegistry.fill(HIST("QA/hDcaXY_after"), track.dcaXY());
    histogramRegistry.fill(HIST("QA/hDcaZ_after"), track.dcaZ());
    histogramRegistry.fill(HIST("QA/hPt"), track.pt());
    histogramRegistry.fill(HIST("QA/hEta"), track.eta());
    histogramRegistry.fill(HIST("QA/hPt_eta"), track.pt(), track.eta());
    histogramRegistry.fill(HIST("QA/hPtDcaXY_after"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("QA/hPtDcaZ_after"), track.pt(), track.dcaZ());
    histogramRegistry.fill(HIST("QA/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("QA/hTPCchi2perCluster_after"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("QA/hITSchi2perCluster_after"), track.itsChi2NCl());
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack())
      return false;
    if (cPVcont && !track.isPVContributor())
      return false;
    if (std::fabs(track.eta()) >= etaCut)
      return false;
    if (track.pt() <= ptMinCut || track.pt() >= ptMaxCut)
      return false;
    if (track.sign() == 0)
      return false;
    if (cDcaXy && std::fabs(track.dcaXY()) >= dcaXYCut)
      return false;
    if (cDcaZ && std::fabs(track.dcaZ()) >= dcaZCut)
      return false;
    if (cTpcCr && track.tpcNClsCrossedRows() <= tpcCrossCut)
      return false;
    if (cItsChi && track.itsChi2NCl() >= itsChiCut)
      return false;
    if (cTpcChi && track.tpcChi2NCl() >= tpcChiCut)
      return false;

    return true;
  }

  double getEfficiency(float pt, int sign)
  {
    TH1D* hEff = nullptr;

    if (sign > 0) {
      hEff = efficiencyPos;
    } else if (sign < 0) {
      hEff = efficiencyNeg;
    }

    if (!hEff) {
      return 1e-6;
    }
    int bin = hEff->GetXaxis()->FindBin(pt);
    if (bin < 1 || bin > hEff->GetNbinsX()) {
      return 1e-6;
    }
    double eff = hEff->GetBinContent(bin);
    return eff;
  }

  void fillHistograms(float nch, float cent, float fpos, float fneg, float posneg, float termp, float termn)
  {
    histogramRegistry.fill(HIST("data/nch"), nch);
    histogramRegistry.fill(HIST("data/cent_nch"), cent, nch);
    histogramRegistry.fill(HIST("data/nch_pos"), fpos);
    histogramRegistry.fill(HIST("data/cent_nch_pos"), cent, fpos);
    histogramRegistry.fill(HIST("data/nch_neg"), fneg);
    histogramRegistry.fill(HIST("data/cent_nch_neg"), cent, fneg);
    histogramRegistry.fill(HIST("data/nch_negpos"), posneg);
    histogramRegistry.fill(HIST("data/cent_nch_negpos"), cent, posneg);

    histogramRegistry.fill(HIST("data/cent_pos"), cent, fpos);
    histogramRegistry.fill(HIST("data/cent_neg"), cent, fneg);
    histogramRegistry.fill(HIST("data/cent_termp"), cent, termp);
    histogramRegistry.fill(HIST("data/cent_termn"), cent, termn);
    histogramRegistry.fill(HIST("data/cent_pos_sq"), cent, fpos * fpos);
    histogramRegistry.fill(HIST("data/cent_neg_sq"), cent, fneg * fneg);
    histogramRegistry.fill(HIST("data/cent_posneg"), cent, posneg);

    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    histogramRegistry.fill(HIST("subsample/pos"), cent, sampleIndex, fpos);
    histogramRegistry.fill(HIST("subsample/neg"), cent, sampleIndex, fneg);
    histogramRegistry.fill(HIST("subsample/termp"), cent, sampleIndex, termp);
    histogramRegistry.fill(HIST("subsample/termn"), cent, sampleIndex, termn);
    histogramRegistry.fill(HIST("subsample/pos_sq"), cent, sampleIndex, fpos * fpos);
    histogramRegistry.fill(HIST("subsample/neg_sq"), cent, sampleIndex, fneg * fneg);
    histogramRegistry.fill(HIST("subsample/posneg"), cent, sampleIndex, posneg);
  }

  template <RunType run, typename C, typename T>
  void calculationData(C const& coll, T const& tracks)
  {
    float cent = -1, mult = -1;
    histogramRegistry.fill(HIST("QA/hVtxZ_before"), coll.posZ());
    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }

    float globalNch = tracks.size();
    float pvTrack = coll.multNTracksPV();

    histogramRegistry.fill(HIST("QA/hCentFT0C"), cent);
    histogramRegistry.fill(HIST("QA/hNchGlobal"), globalNch);
    histogramRegistry.fill(HIST("QA/hNchPV"), pvTrack);

    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_PV_bef"), pvTrack, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_FT0C_bef"), cent, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/PV_FT0C_bef"), cent, pvTrack);

    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent)) {
      return;
    }
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_PV_aft"), pvTrack, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_FT0C_aft"), cent, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/PV_FT0C_aft"), cent, pvTrack);
    histogramRegistry.fill(HIST("QA/hVtxZ_after"), coll.posZ());
    histogramRegistry.fill(HIST("QA/hCentrality"), cent);
    histogramRegistry.fill(HIST("QA/hMultiplicity"), mult);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    double fposW = 0, fnegW = 0, fposW2 = 0.0, fnegW2 = 0.0;
    double posnegW = 0.0, termnW = 0.0, termpW = 0.0;
    int nch = 0;
    double nchCor = 0;
    for (const auto& track : tracks) {

      fillBeforeQA(track);
      if (!selTrack(track))
        continue;

      double eff = getEfficiency(track.pt(), track.sign());
      if (eff < threshold) {
        continue;
      }
      double weight = 1.0 / eff;

      nch += 1;
      fillAfterQA(track);
      histogramRegistry.fill(HIST("QA/cent_hEta"), cent, track.eta());
      histogramRegistry.fill(HIST("QA/cent_hPt"), cent, track.pt());
      histogramRegistry.fill(HIST("data/hPt_cor"), track.pt(), weight);
      histogramRegistry.fill(HIST("data/hEta_cor"), track.eta(), weight);

      if (track.sign() == 1) {
        histogramRegistry.fill(HIST("eff/hPt_np"), track.pt());
      } else if (track.sign() == -1) {
        histogramRegistry.fill(HIST("eff/hPt_nm"), track.pt());
      }

      nchCor += weight;
      if (track.sign() == 1) {
        fpos += 1;
        fposW += weight;
        fposW2 += weight * weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        fnegW += weight;
        fnegW2 += weight * weight;
      }
    } // track
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    termpW = fposW * fposW - fposW2;
    termnW = fnegW * fnegW - fnegW2;
    posnegW = fposW * fnegW;

    histogramRegistry.fill(HIST("data/nch_nchCor"), nch, nchCor);
    histogramRegistry.fill(HIST("data/nchCor"), nchCor);
    histogramRegistry.fill(HIST("data/cent_nchCor"), cent, nchCor);
    histogramRegistry.fill(HIST("data/cent_nch"), cent, nch);
    histogramRegistry.fill(HIST("data/cent_pos_cor"), cent, fposW);
    histogramRegistry.fill(HIST("data/cent_neg_cor"), cent, fnegW);

    if (cfgEffNue) {
      fillHistograms(nchCor, cent, fposW, fnegW, posnegW, termpW, termnW);
    } else {
      fillHistograms(nch, cent, fpos, fneg, posneg, termp, termn);
    }
  }

  template <RunType run, typename C, typename T, typename M, typename P>
  void calculationMc(C const& coll, T const& inputTracks, M const& mcCollisions, P const& mcParticles)
  {
    (void)mcCollisions;
    if (!coll.has_mcCollision()) {
      return;
    }
    histogramRegistry.fill(HIST("gen/hVtxZ_before"), coll.mcCollision().posZ());
    float cent = -1, mult = -1;
    histogramRegistry.fill(HIST("QA/hVtxZ_before"), coll.posZ());
    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }

    int globalNch = inputTracks.size();
    int pvTrack = coll.multNTracksPV();

    histogramRegistry.fill(HIST("QA/hCentFT0C"), cent);
    histogramRegistry.fill(HIST("QA/hNchGlobal"), globalNch);
    histogramRegistry.fill(HIST("QA/hNchPV"), pvTrack);

    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_PV_bef"), pvTrack, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_FT0C_bef"), cent, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/PV_FT0C_bef"), cent, pvTrack);

    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent)) {
      return;
    }
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_PV_aft"), pvTrack, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/globalTracks_FT0C_aft"), cent, globalNch);
    histogramRegistry.fill(HIST("MultCorrelationPlots/PV_FT0C_aft"), cent, pvTrack);

    histogramRegistry.fill(HIST("QA/hVtxZ_after"), coll.posZ());
    histogramRegistry.fill(HIST("QA/hCentrality"), cent);
    histogramRegistry.fill(HIST("QA/hMultiplicity"), mult);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    double fposW = 0.0, fnegW = 0.0, fposW2 = 0.0, fnegW2 = 0.0;
    double posnegW = 0.0, termnW = 0.0, termpW = 0.0;
    int nch = 0;
    double nchCor = 0.0;

    for (const auto& track : inputTracks) {
      fillBeforeQA(track);
      if (!selTrack(track))
        continue;
      nch += 1;
      fillAfterQA(track);

      if (track.sign() == 1) {
        histogramRegistry.fill(HIST("eff/hPt_np"), track.pt());
      } else if (track.sign() == -1) {
        histogramRegistry.fill(HIST("eff/hPt_nm"), track.pt());
      }

      histogramRegistry.fill(HIST("QA/cent_hEta"), cent, track.eta());
      histogramRegistry.fill(HIST("QA/cent_hPt"), cent, track.pt());

      double eff = getEfficiency(track.pt(), track.sign());
      if (eff < threshold)
        continue;
      double weight = 1.0 / eff;
      histogramRegistry.fill(HIST("data/hPt_cor"), track.pt(), weight);
      histogramRegistry.fill(HIST("data/hEta_cor"), track.eta(), weight);

      if (track.sign() == 1) {
        fpos += 1;
        fposW += weight;
        fposW2 += weight * weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        fnegW += weight;
        fnegW2 += weight * weight;
      }
      nchCor += weight;
    } // track
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    termpW = fposW * fposW - fposW2;
    termnW = fnegW * fnegW - fnegW2;
    posnegW = fposW * fnegW;

    histogramRegistry.fill(HIST("data/nch_nchCor"), nch, nchCor);
    histogramRegistry.fill(HIST("data/nchCor"), nchCor);
    histogramRegistry.fill(HIST("data/cent_nchCor"), cent, nchCor);
    histogramRegistry.fill(HIST("data/cent_nch"), cent, nch);
    histogramRegistry.fill(HIST("data/cent_pos_cor"), cent, fposW);
    histogramRegistry.fill(HIST("data/cent_neg_cor"), cent, fnegW);

    if (cfgEffNue) {
      fillHistograms(nchCor, cent, fposW, fnegW, posnegW, termpW, termnW);
    } else {
      fillHistograms(nch, cent, fpos, fneg, posneg, termp, termn);
    }

    int posGen = 0, negGen = 0, posNegGen = 0, termNGen = 0, termPGen = 0, nchGen = 0;

    const auto& mccolgen = coll.template mcCollision_as<aod::McCollisions>();
    if (std::abs(mccolgen.posZ()) >= vertexZcut)
      return;
    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);
    histogramRegistry.fill(HIST("gen/hVtxZ_after"), mccolgen.posZ());
    for (const auto& mcpart : mcpartgen) {
      if (std::fabs(mcpart.eta()) >= etaCut)
        continue;
      if (!mcpart.isPhysicalPrimary())
        continue;
      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdgService->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }
      if (sign == 0)
        continue;
      if (std::abs(pid) != kElectron &&
          std::abs(pid) != kMuonMinus &&
          std::abs(pid) != kPiPlus &&
          std::abs(pid) != kKPlus &&
          std::abs(pid) != kProton)
        continue;
      if (std::fabs(mcpart.eta()) >= etaCut)
        continue;
      if ((mcpart.pt() <= ptMinCut) || (mcpart.pt() >= ptMaxCut)) {
        continue;
      }
      if (sign == 1) {
        histogramRegistry.fill(HIST("eff/hPt_np_gen"), mcpart.pt());
      } else if (sign == -1) {
        histogramRegistry.fill(HIST("eff/hPt_nm_gen"), mcpart.pt());
      }

      histogramRegistry.fill(HIST("gen/hPt"), mcpart.pt());
      histogramRegistry.fill(HIST("gen/cent_hPt"), cent, mcpart.pt());
      histogramRegistry.fill(HIST("gen/hEta"), mcpart.eta());
      histogramRegistry.fill(HIST("gen/cent_hEta"), cent, mcpart.eta());
      histogramRegistry.fill(HIST("gen/hSign"), sign);
      histogramRegistry.fill(HIST("gen/hPt_eta"), mcpart.pt(), mcpart.eta());
      nchGen += 1;
      if (sign == 1) {
        posGen += 1;
      }
      if (sign == -1) {
        negGen += 1;
      }
    }
    termPGen = posGen * (posGen - 1);
    termNGen = negGen * (negGen - 1);
    posNegGen = posGen * negGen;
    histogramRegistry.fill(HIST("gen/cent_pos"), cent, posGen);
    histogramRegistry.fill(HIST("gen/cent_neg"), cent, negGen);
    histogramRegistry.fill(HIST("gen/cent_termp"), cent, termPGen);
    histogramRegistry.fill(HIST("gen/cent_termn"), cent, termNGen);
    histogramRegistry.fill(HIST("gen/cent_pos_sq"), cent, posGen * posGen);
    histogramRegistry.fill(HIST("gen/cent_neg_sq"), cent, negGen * negGen);
    histogramRegistry.fill(HIST("gen/cent_posneg"), cent, posNegGen);
    histogramRegistry.fill(HIST("gen/cent_nch"), cent, nchGen);
    histogramRegistry.fill(HIST("gen/nch"), nchGen);

    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    histogramRegistry.fill(HIST("subsample/gen/pos"), cent, sampleIndex, posGen);
    histogramRegistry.fill(HIST("subsample/gen/neg"), cent, sampleIndex, negGen);
    histogramRegistry.fill(HIST("subsample/gen/termp"), cent, sampleIndex, termPGen);
    histogramRegistry.fill(HIST("subsample/gen/termn"), cent, sampleIndex, termNGen);
    histogramRegistry.fill(HIST("subsample/gen/pos_sq"), cent, sampleIndex, posGen * posGen);
    histogramRegistry.fill(HIST("subsample/gen/neg_sq"), cent, sampleIndex, negGen * negGen);
    histogramRegistry.fill(HIST("subsample/gen/posneg"), cent, sampleIndex, posNegGen);

  } // void

  template <RunType run, typename C, typename T>
  void calculationDeltaEta(C const& coll, T const& tracks, float deta1, float deta2)
  {
    float cent = -1, mult = -1;
    if (!selCollision<run>(coll, cent, mult))
      return;

    int globalNch = tracks.size();
    int pvTrack = coll.multNTracksPV();
    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent))
      return;

    if (!(cent >= centMin && cent < centMax))
      return;
    histogramRegistry.fill(HIST("data/delta_eta_cent"), cent);

    int fpos = 0, fneg = 0, termp = 0, termn = 0, posneg = 0;
    double fposW = 0.0, fnegW = 0.0;
    double fposW2 = 0.0, fnegW2 = 0.0;
    double termpW = 0.0, termnW = 0.0, posnegW = 0.0;

    int nch = 0;
    double nchCor = 0.0;

    for (const auto& track : tracks) {

      if (!selTrack(track))
        continue;

      double eta = track.eta();
      if (eta < deta1 || eta > deta2)
        continue;

      double eff = getEfficiency(track.pt(), track.sign());
      if (eff < threshold)
        continue;

      double weight = 1.0 / eff;

      nch += 1;
      nchCor += weight;

      histogramRegistry.fill(HIST("data/delta_eta_eta"), eta);

      if (track.sign() == 1) {
        fpos += 1;
        fposW += weight;
        fposW2 += weight * weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        fnegW += weight;
        fnegW2 += weight * weight;
      }
    }

    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    termpW = fposW * fposW - fposW2;
    termnW = fnegW * fnegW - fnegW2;
    posnegW = fposW * fnegW;

    float deltaEtaWidth = deta2 - deta1 + 1e-5f;

    histogramRegistry.fill(HIST("data/delta_eta_nch"), deltaEtaWidth, nch);
    histogramRegistry.fill(HIST("data/delta_eta_nchCor"), deltaEtaWidth, nchCor);
    histogramRegistry.fill(HIST("data/delta_eta_pos_cor"), deltaEtaWidth, fposW);
    histogramRegistry.fill(HIST("data/delta_eta_neg_cor"), deltaEtaWidth, fnegW);

    if (cfgEffNue) {
      histogramRegistry.fill(HIST("data/delta_eta_pos"), deltaEtaWidth, fposW);
      histogramRegistry.fill(HIST("data/delta_eta_neg"), deltaEtaWidth, fnegW);
      histogramRegistry.fill(HIST("data/delta_eta_termp"), deltaEtaWidth, termpW);
      histogramRegistry.fill(HIST("data/delta_eta_termn"), deltaEtaWidth, termnW);
      histogramRegistry.fill(HIST("data/delta_eta_pos_sq"), deltaEtaWidth, fposW * fposW);
      histogramRegistry.fill(HIST("data/delta_eta_neg_sq"), deltaEtaWidth, fnegW * fnegW);
      histogramRegistry.fill(HIST("data/delta_eta_posneg"), deltaEtaWidth, posnegW);
    } else {

      histogramRegistry.fill(HIST("data/delta_eta_pos"), deltaEtaWidth, fpos);
      histogramRegistry.fill(HIST("data/delta_eta_neg"), deltaEtaWidth, fneg);
      histogramRegistry.fill(HIST("data/delta_eta_termp"), deltaEtaWidth, termp);
      histogramRegistry.fill(HIST("data/delta_eta_termn"), deltaEtaWidth, termn);
      histogramRegistry.fill(HIST("data/delta_eta_pos_sq"), deltaEtaWidth, fpos * fpos);
      histogramRegistry.fill(HIST("data/delta_eta_neg_sq"), deltaEtaWidth, fneg * fneg);
      histogramRegistry.fill(HIST("data/delta_eta_posneg"), deltaEtaWidth, posneg);
    }
    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    if (cfgEffNue) {
      histogramRegistry.fill(HIST("subsample/delta_eta/pos"), deltaEtaWidth, sampleIndex, fposW);
      histogramRegistry.fill(HIST("subsample/delta_eta/neg"), deltaEtaWidth, sampleIndex, fnegW);
      histogramRegistry.fill(HIST("subsample/delta_eta/termp"), deltaEtaWidth, sampleIndex, termpW);
      histogramRegistry.fill(HIST("subsample/delta_eta/termn"), deltaEtaWidth, sampleIndex, termnW);
      histogramRegistry.fill(HIST("subsample/delta_eta/pos_sq"), deltaEtaWidth, sampleIndex, fposW * fposW);
      histogramRegistry.fill(HIST("subsample/delta_eta/neg_sq"), deltaEtaWidth, sampleIndex, fnegW * fnegW);
      histogramRegistry.fill(HIST("subsample/delta_eta/posneg"), deltaEtaWidth, sampleIndex, posnegW);
    } else {
      histogramRegistry.fill(HIST("subsample/delta_eta/pos"), deltaEtaWidth, sampleIndex, fpos);
      histogramRegistry.fill(HIST("subsample/delta_eta/neg"), deltaEtaWidth, sampleIndex, fneg);
      histogramRegistry.fill(HIST("subsample/delta_eta/termp"), deltaEtaWidth, sampleIndex, termp);
      histogramRegistry.fill(HIST("subsample/delta_eta/termn"), deltaEtaWidth, sampleIndex, termn);
      histogramRegistry.fill(HIST("subsample/delta_eta/pos_sq"), deltaEtaWidth, sampleIndex, fpos * fpos);
      histogramRegistry.fill(HIST("subsample/delta_eta/neg_sq"), deltaEtaWidth, sampleIndex, fneg * fneg);
      histogramRegistry.fill(HIST("subsample/delta_eta/posneg"), deltaEtaWidth, sampleIndex, posneg);
    }
  }

  template <RunType run, typename C, typename T, typename M, typename P>
  void calculationMcDeltaEta(C const& coll, T const& inputTracks, M const& mcCollisions, P const& mcParticles, float deta1, float deta2)
  {
    (void)mcCollisions;

    if (!coll.has_mcCollision())
      return;

    float cent = -1, mult = -1;
    if (!selCollision<run>(coll, cent, mult))
      return;

    int globalNch = inputTracks.size();
    int pvTrack = coll.multNTracksPV();
    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent))
      return;

    if (!(cent >= centMin && cent < centMax))
      return;
    histogramRegistry.fill(HIST("data/delta_eta_cent"), cent);

    float deltaEtaWidth = deta2 - deta1 + 1e-5f;

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    int nch = 0;
    double fposW = 0.0, fnegW = 0.0;
    double fposW2 = 0.0, fnegW2 = 0.0;
    double nchCor = 0.0, termpW = 0.0, termnW = 0.0, posnegW = 0.0;

    for (const auto& track : inputTracks) {
      if (!selTrack(track))
        continue;
      double eta = track.eta();
      if (eta < deta1 || eta > deta2)
        continue;

      histogramRegistry.fill(HIST("data/delta_eta_eta"), eta);
      double eff = getEfficiency(track.pt(), track.sign());
      if (eff < threshold)
        continue;
      double weight = 1.0 / eff;
      nch += 1;
      nchCor += weight;
      if (track.sign() == 1) {
        fpos += 1;
        fposW += weight;
        fposW2 += weight * weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        fnegW += weight;
        fnegW2 += weight * weight;
      }
    } // tracks

    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    termpW = fposW * fposW - fposW2;
    termnW = fnegW * fnegW - fnegW2;
    posnegW = fposW * fnegW;

    histogramRegistry.fill(HIST("data/delta_eta_nch"), deltaEtaWidth, nch);
    histogramRegistry.fill(HIST("data/delta_eta_nchCor"), deltaEtaWidth, nchCor);
    histogramRegistry.fill(HIST("data/delta_eta_pos_cor"), deltaEtaWidth, fposW);
    histogramRegistry.fill(HIST("data/delta_eta_neg_cor"), deltaEtaWidth, fnegW);

    if (cfgEffNue) {
      histogramRegistry.fill(HIST("data/delta_eta_pos"), deltaEtaWidth, fposW);
      histogramRegistry.fill(HIST("data/delta_eta_neg"), deltaEtaWidth, fnegW);
      histogramRegistry.fill(HIST("data/delta_eta_termp"), deltaEtaWidth, termpW);
      histogramRegistry.fill(HIST("data/delta_eta_termn"), deltaEtaWidth, termnW);
      histogramRegistry.fill(HIST("data/delta_eta_pos_sq"), deltaEtaWidth, fposW * fposW);
      histogramRegistry.fill(HIST("data/delta_eta_neg_sq"), deltaEtaWidth, fnegW * fnegW);
      histogramRegistry.fill(HIST("data/delta_eta_posneg"), deltaEtaWidth, posnegW);
    } else {

      histogramRegistry.fill(HIST("data/delta_eta_pos"), deltaEtaWidth, fpos);
      histogramRegistry.fill(HIST("data/delta_eta_neg"), deltaEtaWidth, fneg);
      histogramRegistry.fill(HIST("data/delta_eta_termp"), deltaEtaWidth, termp);
      histogramRegistry.fill(HIST("data/delta_eta_termn"), deltaEtaWidth, termn);
      histogramRegistry.fill(HIST("data/delta_eta_pos_sq"), deltaEtaWidth, fpos * fpos);
      histogramRegistry.fill(HIST("data/delta_eta_neg_sq"), deltaEtaWidth, fneg * fneg);
      histogramRegistry.fill(HIST("data/delta_eta_posneg"), deltaEtaWidth, posneg);
    }

    const auto& mccolgen = coll.template mcCollision_as<aod::McCollisions>();

    if (std::abs(mccolgen.posZ()) >= vertexZcut)
      return;

    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);

    int posGen = 0, negGen = 0, posNegGen = 0, termNGen = 0, termPGen = 0, nchGen = 0;
    for (const auto& mcpart : mcpartgen) {
      if (!mcpart.isPhysicalPrimary())
        continue;

      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdgService->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }
      if (sign == 0)
        continue;
      if (std::abs(pid) != kElectron &&
          std::abs(pid) != kMuonMinus &&
          std::abs(pid) != kPiPlus &&
          std::abs(pid) != kKPlus &&
          std::abs(pid) != kProton)
        continue;

      if (std::fabs(mcpart.eta()) >= etaCut)
        continue;
      if ((mcpart.pt() <= ptMinCut) || (mcpart.pt() >= ptMaxCut))
        continue;

      double mcEta = mcpart.eta();
      if (mcEta < deta1 || mcEta > deta2)
        continue;

      histogramRegistry.fill(HIST("gen/delta_eta_eta"), mcpart.eta());

      nchGen += 1;
      if (sign == 1) {
        posGen += 1;
      }
      if (sign == -1) {
        negGen += 1;
      }
    }

    termPGen = posGen * (posGen - 1);
    termNGen = negGen * (negGen - 1);
    posNegGen = posGen * negGen;

    histogramRegistry.fill(HIST("gen/delta_eta_pos"), deltaEtaWidth, posGen);
    histogramRegistry.fill(HIST("gen/delta_eta_neg"), deltaEtaWidth, negGen);
    histogramRegistry.fill(HIST("gen/delta_eta_termp"), deltaEtaWidth, termPGen);
    histogramRegistry.fill(HIST("gen/delta_eta_termn"), deltaEtaWidth, termNGen);
    histogramRegistry.fill(HIST("gen/delta_eta_pos_sq"), deltaEtaWidth, posGen * posGen);
    histogramRegistry.fill(HIST("gen/delta_eta_neg_sq"), deltaEtaWidth, negGen * negGen);
    histogramRegistry.fill(HIST("gen/delta_eta_posneg"), deltaEtaWidth, posNegGen);
    histogramRegistry.fill(HIST("gen/delta_eta_nch"), deltaEtaWidth, nchGen);

    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    histogramRegistry.fill(HIST("subsample/delta_eta/gen/pos"), deltaEtaWidth, sampleIndex, posGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/neg"), deltaEtaWidth, sampleIndex, negGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/termp"), deltaEtaWidth, sampleIndex, termPGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/termn"), deltaEtaWidth, sampleIndex, termNGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/pos_sq"), deltaEtaWidth, sampleIndex, posGen * posGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/neg_sq"), deltaEtaWidth, sampleIndex, negGen * negGen);
    histogramRegistry.fill(HIST("subsample/delta_eta/gen/posneg"), deltaEtaWidth, sampleIndex, posNegGen);

  } // void

  SliceCache cache;
  Preslice<aod::McParticles> mcTrack = aod::mcparticle::mcCollisionId;

  // process function for Data Run3
  void processDataRun3(MyCollisionRun3 const& coll, MyTracks const& tracks)
  {
    calculationData<kRun3>(coll, tracks);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1);
      float etaMax = 0.1f * (ii + 1);

      calculationDeltaEta<kRun3>(coll, tracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun3, "Process for Run3 DATA", false);

  // process function for Data Run2
  void processDataRun2(MyCollisionRun2 const& coll, MyTracks const& tracks)
  {
    calculationData<kRun2>(coll, tracks);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1); // -0.1, -0.2, ..., -0.8
      float etaMax = 0.1f * (ii + 1);  // +0.1, +0.2, ..., +0.8

      calculationDeltaEta<kRun2>(coll, tracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun2, "Process for Run2 DATA", false);

  // process function for MC Run3

  void processMcRun3(MyMCCollisionRun3 const& coll, MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun3>(coll, inputTracks, mcCollisions, mcParticles);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1);
      float etaMax = 0.1f * (ii + 1);
      calculationMcDeltaEta<kRun3>(coll, inputTracks, mcCollisions, mcParticles, etaMin, etaMax);
    }
  }
  PROCESS_SWITCH(NetchargeFluctuations, processMcRun3, "Process reconstructed", true);

  // process function for MC Run2

  void processMcRun2(MyMCCollisionRun2 const& coll, MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun2>(coll, inputTracks, mcCollisions, mcParticles);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1);
      float etaMax = 0.1f * (ii + 1);
      calculationMcDeltaEta<kRun2>(coll, inputTracks, mcCollisions, mcParticles, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processMcRun2, "Process reconstructed", false);
};

// struct
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    {adaptAnalysisTask<NetchargeFluctuations>(cfgc)}};
}
