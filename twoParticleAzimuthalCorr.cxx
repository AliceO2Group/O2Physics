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

/// \file twoParticleAzimuthalCorr.cxx
/// \brief Two-particle azimuthal correlation for pp, O-O, p-O, Ne-Ne and Pb-Pb collisions
/// \author Mintu Haldar (mintu.haldar@cern.ch)
/// \since  March/12/2026

// O2 includes
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};
R__LOAD_LIBRARY(libO2PhysicsPWGCFCore)

// STEP 2
// Example task illustrating how to mix elements of different partitions and different events + process switches

struct twoParticleAzimuthalCorr {

  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgCutVtxZ, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 20.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.0, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, true, "cut time intervals with dead ITS staves")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, bool, false, "Verbose output")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, true, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, true, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgSoloPtTrack, bool, false, "Skip trigger tracks that are alone in their pT bin for same process")
  O2_DEFINE_CONFIGURABLE(cfgSingleSoloPtTrack, bool, false, "Skip associated tracks that are alone in their pT bin for same process, works only if cfgSoloPtTrack is enabled")
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultT0CCutLow = nullptr;
    TF1* fMultT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
    TF1* fMultMultV0ACutLow = nullptr;
    TF1* fMultMultV0ACutHigh = nullptr;
    TF1* fT0AV0AMean = nullptr;
    TF1* fT0AV0ASigma = nullptr;
  } cfgFuncParas;
  /// Function to aid in calculating delta-phi
  // /// \param phi1 first phi value
  // /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }
  SliceCache cache;

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = true;

  // Define the outputs
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Defining configurables

  Configurable<float> ConfMinNSigmaTPCCut{"ConfMinNSigmaTPCCut", 3., "N-sigma TPC cut"};
  Configurable<float> ConfChargeCut{"ConfChargeCut", 0., "N-sigma TPC cut"};

  Configurable<float> associatedMinPt{"associatedMinPt", 0.2f, "NSassociatedMinPt"};
  Configurable<float> associatedMaxPt{"associatedMaxPt", 6.0f, "associatedMaxPt"};
  Configurable<float> triggerMinPt{"triggerMinPt", 6.0f, "triggerMinPt"};
  // --- Mixing configuration
  Configurable<int> nMixEvents{"nMixEvents", 5, "Number of previous events to mix with"};
  Configurable<int> nMixZBins{"nMixZBins", 10, "Number of vertex-Z bins for mixing"};
  Configurable<float> minMixZ{"minMixZ", -10.0f, "Minimum vertex Z for mixing"};
  Configurable<float> maxMixZ{"maxMixZ", +10.0f, "Maximum vertex Z for mixing"};

  // lightweight track struct used for storing tracks into the mixing buffer
  struct SimpleTrack {
    float phi;
    float eta;
    float pt;
  };

  // mixBuffer[zBin] is a deque of previous events; each event is a vector<SimpleTrack>
  std::vector<std::deque<std::vector<SimpleTrack>>> mixBuffer;

  // Defining filters
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVtxZ);
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > 0.2f) && (aod::track::pt < 20.0f);

  // Applying filters
  using MyTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>;
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  using FilteredTracksWithMCLabels = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>>;

  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using FilteredMcParticles = soa::Filtered<aod::McParticles>;
  using McParticlesFull = aod::McParticles;

  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVtxZ;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;

  using SmallGroupMcCollisions = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;

  Partition<MyFilteredTracks> positive = aod::track::signed1Pt > ConfChargeCut;
  Partition<MyFilteredTracks> negative = aod::track::signed1Pt < ConfChargeCut;

  Partition<MyFilteredTracks> triggerTracks = aod::track::pt > triggerMinPt;
  Partition<MyFilteredTracks> associatedTracks = aod::track::pt > associatedMinPt&& aod::track::pt <= associatedMaxPt;

  Preslice<MyTracks> perCol = aod::track::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> collisionPerMCCollision = aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -0.5 * o2::constants::math::PI, +1.5 * o2::constants::math::PI}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 13, 16, 20}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {200, 0.0f, 20.0f}, "pt axis"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for efficiency histograms"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0A>;

  // define global variables
  TRandom3* gRandom = new TRandom3();
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    if (cfgCentTableUnavailable && !cfgSelCollByNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }

    const AxisSpec axisPhi{72, 0.0, 2 * o2::constants::math::PI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    histos.add("ptAssoHistogram", "ptAssoHistogram", kTH1D, {axisPtQA});
    histos.add("ptTrigHistogram", "ptTrigHistogram", kTH1D, {axisPtQA});

    histos.add("correlationFunction", "correlationFunction", kTH1D, {axisDeltaPhi});
    histos.add("correlationFunctionO2", "correlationFunctionO2", kTH1D, {axisDeltaPhi});

    histos.add("correlationFunction2d", "correlationFunction2d", kTH2F, {axisDeltaPhi, axisDeltaEta});

    // mixed-event histograms
    histos.add("correlationFunctionMixed", "correlationFunctionMixed", kTH1F, {axisDeltaPhi});
    histos.add("correlationFunction2dMixed", "correlationFunction2dMixed", kTH2F, {axisDeltaPhi, axisDeltaEta});
    histos.add("associatedPt", "associatedPt", kTH1D, {axisPtQA});

    // prepare mixing buffers (one deque per z-bin)
    mixBuffer.resize(nMixZBins);

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx", ";Z (cm)", kTH1F, {axisVertex});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEta", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {axisPtQA});
    histos.add("hNsigmaTPCP", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{35, 0.5, 4.}, {100, -5., 5.}});
    histos.add("hChargePos", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hChargeNeg", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hInvariantMass", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});
    histos.add("hInvariantMassMixed", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});
    histos.add("hInvariantMassMixedInterface", ";M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 1.0}});

    // Event Counter
    if (doprocessSame && cfgUseAdditionalEventCut) {
      histos.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{12, 0, 12}}});
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "occupancy");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(11, "MultCorrelation");
      histos.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(12, "cfgEvSelV0AT0ACut");
    }

    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    // Make histograms to check the distributions after cuts
    if (doprocessSame) {
      histos.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      histos.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
      histos.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      histos.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      histos.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
      histos.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      histos.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      histos.add("CentralityWeighted", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      histos.add("Centrality_used", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}}); // histogram to see how many events are in the same and mixed event
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("zVtx_used", "zVtx_used", {HistType::kTH1D, {axisVertex}});
      histos.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      histos.add("SE_correlation", "", {HistType::kTHnSparseF, {{axisDeltaPhi, axisDeltaEta, axisPtTrigger, axisPtAssoc, axisCentrality}}});
      histos.add("ME_correlation", "", {HistType::kTHnSparseF, {{axisDeltaPhi, axisDeltaEta, axisPtTrigger, axisPtAssoc, axisCentrality}}});
    }

    //  if (doprocessMCSame) {
    //   LOGF(fatal, "Full simulation and on-the-fly processing of same event not supported");
    // }
    // if (doprocessMCMixed) {
    //   LOGF(fatal, "Full simulation and on-the-fly processing of mixed event not supported");
    // }
    if (doprocessMCSame) {
      histos.add("MCTrue/MCeventcount", "MCeventcount", {HistType::kTH1F, {{5, 0, 5, "bin"}}}); // histogram to see how many events are in the same and mixed event
      histos.get<TH1>(HIST("MCTrue/MCeventcount"))->GetXaxis()->SetBinLabel(2, "same all");
      histos.get<TH1>(HIST("MCTrue/MCeventcount"))->GetXaxis()->SetBinLabel(3, "same reco");
      histos.get<TH1>(HIST("MCTrue/MCeventcount"))->GetXaxis()->SetBinLabel(4, "mixed all");
      histos.get<TH1>(HIST("MCTrue/MCeventcount"))->GetXaxis()->SetBinLabel(5, "mixed reco");
      histos.add("MCTrue/MCCentrality", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}});
      histos.add("MCTrue/MCNch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("MCTrue/MCzVtx", "MCzVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("MCTrue/MCPhi", "MCPhi", {HistType::kTH1D, {axisPhi}});
      histos.add("MCTrue/MCEta", "MCEta", {HistType::kTH1D, {axisEta}});
      histos.add("MCTrue/MCpT", "MCpT", {HistType::kTH1D, {axisPtTrigger}});
      histos.add("MCTrue/MCTrig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      histos.add("MCTrue/MCdeltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      histos.add("MCTrue/MCdeltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
      histos.add("MCdeltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      histos.add("MCdeltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
      histos.add("Trig_hist_MC", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      histos.add("SE_correlation_MC", "", {HistType::kTHnSparseF, {{axisDeltaPhi, axisDeltaEta, axisPtTrigger, axisPtAssoc}}});
      histos.add("ME_correlation_MC", "", {HistType::kTHnSparseF, {{axisDeltaPhi, axisDeltaEta, axisPtTrigger, axisPtAssoc}}});
    }
    if (doprocessMCEfficiency) {
      histos.add("MCEffeventcount", "bin", {HistType::kTH1F, {{5, 0, 5, "bin"}}});
      histos.get<TH1>(HIST("MCEffeventcount"))->GetXaxis()->SetBinLabel(1, "All");
      histos.get<TH1>(HIST("MCEffeventcount"))->GetXaxis()->SetBinLabel(2, "MC");
      histos.get<TH1>(HIST("MCEffeventcount"))->GetXaxis()->SetBinLabel(3, "Reco Primary");
      histos.get<TH1>(HIST("MCEffeventcount"))->GetXaxis()->SetBinLabel(4, "Reco All");
      histos.get<TH1>(HIST("MCEffeventcount"))->GetXaxis()->SetBinLabel(5, "Fake");
    }

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    LOGF(info, "End of init");
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgCutITSclu));
  }

  template <typename TTrack>
  bool genTrackSelected(TTrack track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (!track.producedByGenerator()) {
      return false;
    }
    if (std::abs(track.eta()) > cfgCutEta) {
      return false;
    }
    if (std::abs(track.pt()) < cfgCutPtMin || std::abs(track.pt()) > cfgCutPtMax) {
      return false;
    }
    return true;
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = reinterpret_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    if (cfgCentralityWeight.value.empty() == false) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgCentralityWeight.value.c_str(), (void*)mCentralityWeight);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  bool getCentralityWeight(float& weightCent, const float centrality)
  {
    float weight = 1.;
    if (mCentralityWeight)
      weight = mCentralityWeight->GetBinContent(mCentralityWeight->FindBin(centrality));
    else
      weight = 1.0;
    if (weight == 0)
      return false;
    weightCent = weight;
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision const& collision, TTracks const& tracks) // function to fill the yield and etaphi histograms.
  {
    // float weff1 = 1;
    float vtxz = collision.posZ();
    histos.fill(HIST("hZvtx"), vtxz);

    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      float weff1 = 1.f;
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), vtxz))
        continue;
      histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      histos.fill(HIST("Eta"), track1.eta());
      histos.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      histos.fill(HIST("pT"), track1.pt());
      histos.fill(HIST("pTCorrected"), track1.pt(), weff1);

      // LOGF(info, "///////////////////////////////////////////////////////////////////////////////////");
    }
  }

  template <typename TTrack, typename TTrackAssoc>
  float getdeltaPhi(TTrack const& track1, TTrackAssoc const& track2, float radius, int magField)
  {
    float charge1 = track1.sign();
    float charge2 = track2.sign();

    float phi1 = track1.phi();
    float phi2 = track2.phi();

    float pt1 = track1.pt();
    float pt2 = track2.pt();

    int fbSign = (magField > 0) ? 1 : -1;

    /// \param phi1 first phi value
    /// \param phi2 second phi value

    Double_t getdeltaPhi = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);
    if (getdeltaPhi < -TMath::Pi() / 2.) {
      getdeltaPhi += 2. * TMath::Pi();
    }
    if (getdeltaPhi > 3 * TMath::Pi() / 2.) {
      getdeltaPhi -= 2. * TMath::Pi();
    }
    return getdeltaPhi;
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    if (system == SameEvent) {
      if (!cfgCentTableUnavailable) {
        histos.fill(HIST("Centrality_used"), cent);
      }
    }

    // Fill centrality histogram
    histos.fill(HIST("Centrality"), cent);

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;
      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;

        if (!cfgUsePtOrder && track1.globalIndex() == track2.globalIndex())
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -0.5 * o2::constants::math::PI);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getdeltaPhi(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getdeltaPhi(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getdeltaPhi(track1, track2, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }
        // LOGF(info, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4");

        // fill the right sparse and histograms
        if (system == SameEvent) {

          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("SE_correlation"), deltaPhi, deltaEta, track1.pt(), track2.pt(), cent);

        } else if (system == MixedEvent) {

          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("ME_correlation"), deltaPhi, deltaEta, track1.pt(), track2.pt(), cent);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillMCCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {
      if (step >= CorrelationContainer::kCFStepTrackedOnlyPrim && !track1.isPhysicalPrimary())
        continue;
      if (!genTrackSelected(track1))
        continue;

      if (system == SameEvent && (doprocessMCSame))
        histos.fill(HIST("MCTrue/MCTrig_hist"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);

      for (auto const& track2 : tracks2) {

        if (step >= CorrelationContainer::kCFStepTrackedOnlyPrim && !track2.isPhysicalPrimary())
          continue;
        if (!genTrackSelected(track2))
          continue;

        if (!cfgUsePtOrder && track1.globalIndex() == track2.globalIndex())
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        // fill the right sparse and histograms
        if (system == SameEvent) {
          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("MCdeltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("SE_correlation_MC"), deltaPhi, deltaEta, track1.pt(), track2.pt());

          if (doprocessMCSame)
            histos.fill(HIST("MCTrue/MCdeltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {
          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("MCdeltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("ME_correlation_MC"), deltaPhi, deltaEta, track1.pt(), track2.pt());

          if (doprocessMCMixed)
            histos.fill(HIST("MCTrue/MCdeltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float centrality, const bool fillCounter)
  {
    histos.fill(HIST("hEventCountSpecific"), 0.5);
    if (cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && cfgEvSelkNoSameBunchPileup)
      histos.fill(HIST("hEventCountSpecific"), 1.5);
    if (cfgEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEvSelkNoITSROFrameBorder)
      histos.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEvSelkNoTimeFrameBorder)
      histos.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && cfgEvSelkIsGoodZvtxFT0vsPV)
      histos.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && cfgEvSelkNoCollInTimeRangeStandard)
      histos.fill(HIST("hEventCountSpecific"), 5.5);
    if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && cfgEvSelkIsGoodITSLayersAll)
      histos.fill(HIST("hEventCountSpecific"), 6.5);
    if (cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && cfgEvSelkNoCollInRofStandard)
      histos.fill(HIST("hEventCountSpecific"), 7.5);
    if (cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && cfgEvSelkNoHighMultCollInPrevRof)
      histos.fill(HIST("hEventCountSpecific"), 8.5);
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh))
      return 0;
    if (fillCounter && cfgEvSelOccupancy)
      histos.fill(HIST("hEventCountSpecific"), 9.5);

    auto multNTracksPV = collision.multNTracksPV();
    if (cfgEvSelMultCorrelation) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled && !cfgCentTableUnavailable) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(centrality))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(centrality))
          return 0;
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled && !cfgCentTableUnavailable) {
        if (multTrk < cfgFuncParas.fMultT0CCutLow->Eval(centrality))
          return 0;
        if (multTrk > cfgFuncParas.fMultT0CCutHigh->Eval(centrality))
          return 0;
      }
      if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
        if (multTrk < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
          return 0;
        if (multTrk > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
          return 0;
      }
      if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(multTrk))
          return 0;
        if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(multTrk))
          return 0;
      }
    }
    if (fillCounter && cfgEvSelMultCorrelation)
      histos.fill(HIST("hEventCountSpecific"), 10.5);

    // V0A T0A 5 sigma cut
    float sigma = 5.0;
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && cfgEvSelV0AT0ACut)
      histos.fill(HIST("hEventCountSpecific"), 11.5);

    return 1;
  }

  void processSame(MyFilteredCollisions::iterator const& collision, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    float weightCent = 1.0f;
    float weff1 = 1.0f;
    float vtxz = collision.posZ();
    if (!cfgCentTableUnavailable) {
      cent = getCentrality(collision);
    }
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true))
      return;
    loadCorrection(bc.timestamp());
    if (!cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      histos.fill(HIST("Centrality"), cent);
      histos.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    auto groupPositive = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    histos.fill(HIST("hZvtx"), collision.posZ());

    for (auto track : groupPositive) {
      histos.fill(HIST("hChargePos"), track.sign());
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      // histos.fill(HIST("hNsigmaTPCP"), track.p(), track.tpcNSigmaPi());
    }

    for (auto track : groupNegative) {
      histos.fill(HIST("hChargeNeg"), track.sign());
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      // histos.fill(HIST("hNsigmaTPCP"), track.p(), track.tpcNSigmaPi());
    }

    for (auto const& track : tracks) {
      histos.fill(HIST("pT"), track.pt());
    }

    // LOGF(info, "--------------------tracks size--------------------: %d", tracks.size());

    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), vtxz))
        continue;
      histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      histos.fill(HIST("Eta"), track1.eta());
      histos.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      histos.fill(HIST("pT"), track1.pt());
      histos.fill(HIST("pTCorrected"), track1.pt(), weff1);

      // LOGF(info, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
    }

    histos.fill(HIST("zVtx"), SameEvent); // because its same event i put it in the 1 bin

    auto assoTracksThisCollision = associatedTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto trigTracksThisCollision = triggerTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto& track : assoTracksThisCollision) {
      // LOGF(info, "================================================================================");
      histos.fill(HIST("ptAssoHistogram"), track.pt());
    }
    for (auto& track : trigTracksThisCollision) {
      histos.fill(HIST("ptTrigHistogram"), track.pt());
    }

    for (auto& trigger : trigTracksThisCollision) {
      for (auto& associated : assoTracksThisCollision) {
        histos.fill(HIST("correlationFunction"), ComputeDeltaPhi(trigger.phi(), associated.phi()));
      }
    }

    for (auto& [trigger, associated] : combinations(o2::soa::CombinationsFullIndexPolicy(trigTracksThisCollision, assoTracksThisCollision))) {
      histos.fill(HIST("correlationFunctionO2"), ComputeDeltaPhi(trigger.phi(), associated.phi()));
      histos.fill(HIST("correlationFunction2d"), ComputeDeltaPhi(trigger.phi(), associated.phi()), trigger.eta() - associated.eta());
    }

    fillYield(collision, trigTracksThisCollision);

    same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);

    // if (!cfgSoloPtTrack) {
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(trigTracksThisCollision, assoTracksThisCollision, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()), cent, weightCent);
    // }
    // else {
    //   fillCorrelationsExcludeSoloTracks<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), getMagneticField(bc.timestamp()), cent, weightCent);
    // }
  }
  PROCESS_SWITCH(twoParticleAzimuthalCorr, processSame, "Process same event", false);

  // the process for filling the mixed events
  void processMixed(MyFilteredCollisions const& collisions, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](MyFilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<MyFilteredCollisions, MyFilteredTracks, MyFilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      histos.fill(HIST("zVtx"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      if (!cfgCentTableUnavailable)
        getCentralityWeight(weightCent, cent1);

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()), cent1, eventWeight * weightCent);
    }
  }

  PROCESS_SWITCH(twoParticleAzimuthalCorr, processMixed, "Process mixed events", false);

  int getSpecies(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus: // pion
        return 0;
      case PDG_t::kKPlus: // Kaon
        return 1;
      case PDG_t::kProton: // proton
        return 2;
      default: // NOTE. The efficiency histogram is hardcoded to contain 4 species. Anything special will have the last slot.
        return 3;
    }
  }

  void processMCEfficiency(FilteredMcCollisions::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, McParticlesFull const& mcParticles, FilteredTracksWithMCLabels const& tracks)
  {
    histos.fill(HIST("MCEffeventcount"), 0.5);
    if (cfgSelCollByNch && (mcParticles.size() < cfgCutMultMin || mcParticles.size() >= cfgCutMultMax)) {
      return;
    }
    // Primaries
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("MCEffeventcount"), 1.5);
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), 0., mcCollision.posZ());
      }
    }
    for (const auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCol, collision.globalIndex());
      if (cfgVerbosity) {
        LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
        LOGF(info, "  which has %d tracks", groupedTracks.size());
      }

      for (const auto& track : groupedTracks) {
        if (track.has_mcParticle()) {

          auto mcId = track.mcParticleId();               //  get index
          auto mcParticle = mcParticles.iteratorAt(mcId); //  correct binding

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("MCEffeventcount"), 2.5);
            same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), 0., mcCollision.posZ());
          }
          histos.fill(HIST("MCEffeventcount"), 3.5);
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), 0., mcCollision.posZ());
        } else {
          // fake track
          histos.fill(HIST("MCEffeventcount"), 4.5);
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, track.eta(), track.pt(), 0, 0., mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(twoParticleAzimuthalCorr, processMCEfficiency, "MC: Extract efficiencies", true);

  // LOGF(info, "#######################################################################");

  void processMCSame(FilteredMcCollisions::iterator const& mcCollision, FilteredMcParticles const& mcParticles, SmallGroupMcCollisions const& collisions)
  {
    if (cfgVerbosity) {
      LOGF(info, "processMCSame. MC collision: %d, particles: %d, collisions: %d", mcCollision.globalIndex(), mcParticles.size(), collisions.size());
    }

    LOGF(info, "======================================================================");

    float cent = -1;
    if (!cfgCentTableUnavailable) {
      for (const auto& collision : collisions) {
        cent = getCentrality(collision);
      }
    }

    if (cfgSelCollByNch && (mcParticles.size() < cfgCutMultMin || mcParticles.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    histos.fill(HIST("MCTrue/MCeventcount"), SameEvent); // because its same event i put it in the 1 bin
    if (!cfgCentTableUnavailable)
      histos.fill(HIST("MCTrue/MCCentrality"), cent);
    histos.fill(HIST("MCTrue/MCNch"), mcParticles.size());
    histos.fill(HIST("MCTrue/MCzVtx"), mcCollision.posZ());
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("MCTrue/MCPhi"), mcParticle.phi());
        histos.fill(HIST("MCTrue/MCEta"), mcParticle.eta());
        histos.fill(HIST("MCTrue/MCpT"), mcParticle.pt());
      }
    }

    // auto assoTracksThisCollisionMC = associatedTracks->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    // auto trigTracksThisCollisionMC = triggerTracks->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    same->fillEvent(mcParticles.size(), CorrelationContainer::kCFStepAll);

    fillMCCorrelations<CorrelationContainer::kCFStepAll>(mcParticles, mcParticles, mcCollision.posZ(), SameEvent, 1.0f);
    // fillMCCorrelations<CorrelationContainer::kCFStepAll>(trigTracksThisCollisionMC, assoTracksThisCollisionMC, mcCollision.posZ(), SameEvent, 1.0f);

    if (collisions.size() == 0) {
      return;
    }

    histos.fill(HIST("MCTrue/MCeventcount"), 2.5);
    same->fillEvent(mcParticles.size(), CorrelationContainer::kCFStepTrackedOnlyPrim);
    fillMCCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(mcParticles, mcParticles, mcCollision.posZ(), SameEvent, 1.0f);
    // fillMCCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(trigTracksThisCollisionMC, assoTracksThisCollisionMC, mcCollision.posZ(), SameEvent, 1.0f);
  }
  PROCESS_SWITCH(twoParticleAzimuthalCorr, processMCSame, "Process MC same event", true);

  void processMCMixed(FilteredMcCollisions const& mcCollisions, FilteredMcParticles const& mcParticles, SmallGroupMcCollisions const& collisions)
  {
    auto getTracksSize = [&mcParticles, this](FilteredMcCollisions::iterator const& mcCollision) {
      auto associatedTracks = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, o2::aod::mccollision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<FilteredMcCollisions, FilteredMcParticles, FilteredMcParticles, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, mcCollisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, collision1.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "Found %d related collisions", groupedCollisions.size());
      }
      float cent = -1;
      if (!cfgCentTableUnavailable) {
        for (const auto& collision : groupedCollisions) {
          cent = getCentrality(collision);
        }
      }

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && groupedCollisions.size() != 0 && (cent < cfgCutCentMin || cent >= cfgCutCentMax))
        continue;

      histos.fill(HIST("MCTrue/MCeventcount"), MixedEvent); // fill the mixed event in the 3 bin
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      fillMCCorrelations<CorrelationContainer::kCFStepAll>(tracks1, tracks2, collision1.posZ(), MixedEvent, eventWeight);

      if (groupedCollisions.size() == 0) {
        continue;
      }

      histos.fill(HIST("MCTrue/MCeventcount"), 4.5);
      fillMCCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(tracks1, tracks2, collision1.posZ(), MixedEvent, eventWeight);
    }
  }
  PROCESS_SWITCH(twoParticleAzimuthalCorr, processMCMixed, "Process MC mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<twoParticleAzimuthalCorr>(cfgc)};
  return workflow;
}
