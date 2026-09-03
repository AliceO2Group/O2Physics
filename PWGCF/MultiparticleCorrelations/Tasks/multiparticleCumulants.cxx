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

/// \file multiparticleCumulants.cxx
/// \brief Task for producing multiparticle cumulants
/// \author Pei-Ying Kuan, TU München, pei-ying.kuan@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TCollection.h>
#include <TComplex.h>
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TParticlePDG.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <RtypesCore.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator;
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

using namespace std;

// *) Define enums:
enum EnRlMc {
  eRl = 0,
  eMc
};

enum EnRecSim {
  eRec = 0,
  eSim,
  eRecAndSim
};

enum EnProcess {
  eProcessRec = 0, // Run 3, only reconstructed
  eProcessRecSim,  // Run 3, both reconstructed and simulated
  eProcessSim,     // Run 3, only simulated
  eProcess_N
};

enum EnEventHistograms {
  eCent,
  eMult,
  eVertexX,
  eVertexY,
  eVertexZ,
  eNumContrib,
  eEventHistograms_N
};

static constexpr std::array<const char*, eEventHistograms_N> EventHistNames = {
  "Centrality",
  "Multiplicity",
  "VertexX",
  "VertexY",
  "VertexZ",
  "NumContrib"};

enum EnParticleHistograms {
  ePt,
  ePhi,
  eParticleHistograms_N
};

static constexpr std::array<const char*, eParticleHistograms_N> ParticleHistNames = {
  "Pt",
  "Phi"};

enum EnQAHistograms {
  eQACent,
  eQAMultNumContrib,
  eQAHistograms_N
};

enum EnCorrHistograms {
  eCorrCent,
  eCorrMult,
  eCorrNumContribMult,
  eCorrHistograms_N
};

static constexpr std::array<const char*, eCorrHistograms_N> CorrHistNames = {
  "Centrality",
  "Multiplicity",
  "NContributionMultiplicity"};

enum EnCorrBound {
  eCentABBound,
  eCentACBound,
  eCentBCBound,
  eMultABBound,
  eMultACBound,
  eMultBCBound,
  eNumContribMultBound,
  eCorrBound_N
};

enum EnCentEstm {
  eCentFT0C,
  eCentFT0M,
  eCentFV0A,
  eCentEstm_N
};

static constexpr std::array<const char*, eCentEstm_N> CentEstmNames = {
  "FT0C",
  "FT0M",
  "FV0A"};

enum EnMultEstm {
  eMultFT0C,
  eMultFT0M,
  eMultFV0A,
  eMultEstm_N
};

static constexpr std::array<const char*, eMultEstm_N> MultEstmNames = {
  "FT0C",
  "FT0M",
  "FV0A"};

enum EnEstmCorr {
  eFT0CFT0M,
  eFT0CFV0A,
  eFT0MFV0A,
  eEstmCorr_N
};

static constexpr std::array<const char*, eEstmCorr_N> EstmCorrName = {
  "FT0CFT0M",
  "FT0CFV0A",
  "FT0MFV0A"};

enum EnBeforeAfter {
  eBefore,
  eAfter,
  eBeforeAfter_N
};

static constexpr std::array<const char*, eBeforeAfter_N * 2> BeforeAfterNames = {
  "before",
  "after",
  "Before",
  "After"};

static constexpr int NumHarmonics = 3;                                                               // n = 2, 3, 4
static constexpr int NumTwoPCorrBins = NumHarmonics;                                                 // <v2^2>, <v3^2>, <v4^2>, in this case 3
static constexpr int NumFourPCorrBins = (NumHarmonics * (NumHarmonics + 1)) / 2;                     // <v2^2 v2^2>, <v2^2 v3^2>, <v2^2 v4^2>, <v3^2 v3^2>, <v3^2 v4^2>, <v4^2 v4^2>, with <vn^2 vn^2>, in this case 6
static constexpr int NumSixPCorrBins = (NumHarmonics * (NumHarmonics - 1) * (NumHarmonics - 2)) / 6; // <v2^2 v3^2 v4^2>, without <vn^2 vn^2 vm^2>, in this case 1
static constexpr std::array<std::array<int, 2>, NumFourPCorrBins> FourPHarmonicIndex{{{2, 2}, {2, 3}, {2, 4}, {3, 3}, {3, 4}, {4, 4}}};

// *) Main task:
struct MultiparticleCumulants { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name";
  OutputObj<TList> fBaseList{sBaseListName.Data(), OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  // *) Service:
  Service<ccdb::BasicCCDBManager> ccdb{}; // support for offline callibration data base
  Service<o2::framework::O2DatabasePDG> pdg{};

  // *) Define configurables:
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};
  Configurable<std::string> cfCentEstm{"cfCentEstm", "FT0M", "centrality estimator: FT0C, FT0M, FV0A"};
  Configurable<std::string> cfMultEstm{"cfMultEstm", "FV0A", "multiplicity estimator: FT0C, FT0M, FV0A"};

  Configurable<bool> cfQASwitch{"cfQASwitch", true, "quality assurance switch"};
  Configurable<bool> cfWeightSwitch{"cfWeightSwitch", true, "weight switch"};
  Configurable<bool> cfPrintSwitch{"cfPrintSwitch", true, "printing result switch"};

  // *) Event cut switches
  Configurable<bool> cfGeneralCutSwitch{"cfGeneralCutSwitch", true, "NoCollInTimeRangeStandard, NoCollInRofStandard, NoSameBunchPileup, IsVertexITSTPC, IsGoodITSLayersAll, IsGoodZvtxFT0vsPV, NoHighMultCollInPrevRof"};
  Configurable<bool> cfVertexZCutSwitch{"cfVertexZCutSwitch", true, "vertex z cut switch"};
  Configurable<bool> cfSel8CutSwitch{"cfSel8CutSwitch", true, "Sel8 cut switch"};
  Configurable<bool> cfCentCutSwitch{"cfCentCutSwitch", true, "centrality cut switch"};
  Configurable<bool> cfNumContribCutSwitch{"cfNumContribCutSwitch", true, "NContribution cut switch"};
  Configurable<bool> cfCentCorrCutSwitch{"cfCentCorrCutSwitch", true, "centrality correlation outlier cut switch"};
  Configurable<bool> cfMultCorrCutSwitch{"cfMultCorrCutSwitch", true, "multiplicity correlation outlier cut switch"};
  Configurable<bool> cfNumContribMultCorrCutSwitch{"cfNumContribMultCorrCutSwitch", true, "NContribution vs. multiplicity correlation outlier cut switch"};

  // *) Particle cut switches
  Configurable<bool> cfPtCutSwitch{"cfPtCutSwitch", true, "Pt cut switch"};
  Configurable<bool> cfEtaCutSwitch{"cfEtaCutSwitch", true, "Eta cut switch"};
  Configurable<bool> cfSignCutSwitch{"cfSignCutSwitch", true, "Charge cut switch"};
  Configurable<bool> cfTpcNClsFoundCutSwitch{"cfTpcNClsFoundCutSwitch", true, ""};
  Configurable<bool> cfDCAXYCutSwitch{"cfDCAXYCutSwitch", true, ""};
  Configurable<bool> cfDCAZCutSwitch{"cfDCAZCutSwitch", true, ""};
  Configurable<bool> cfEtaGapSwitch{"cfEtaGapSwitch", true, "Eta gap switch"};

  // *) Event cut
  Configurable<std::vector<float>> cfVertexZCut{"cfVertexZCut", {-10., 10.}, "vertex z position range: {min, max}[cm]"};
  Configurable<std::vector<float>> cfCentCut{"cfCentCut", {10., 20.}, "centrality range: {min, max}[%]"};
  Configurable<std::vector<float>> cfNumContribCut{"cfNumContribCut", {0, 3000.}, "NContribution range: {min, max}"};
  // Configurable<std::vector<float>> cfCentCorrCut{"cfCentCorrCut", {1., 10.}, "parameters of limits in centralities 2D histograms, (x-t)/m < y < mx+t: {m, t}"};
  Configurable<std::vector<float>> cfCentABCorrCut{"cfCentABCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0C vs. FT0M centrality 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfCentACCorrCut{"cfCentACCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0C vs. FT0M centrality 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfCentBCCorrCut{"cfCentBCCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0C vs. FT0M centrality 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfMultABCorrCut{"cfMultABCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0C vs. FT0M multiplicities 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfMultACCorrCut{"cfMultACCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0C vs. FV0A multiplicities 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfMultBCCorrCut{"cfMultBCCorrCut", {1.4, 300., 1.4, 300.}, "parameters of limits in FT0M vs. FV0A multiplicities 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};
  Configurable<std::vector<float>> cfNumContribMultCorrCut{"cfNumContribMultCorrCut", {0.035, 100., 0.024, -600.}, "parameters of limits in NContribution vs. multiplicity 2D histograms, cx+d < y < ax+b: {a,b,c,d}"};

  // *) Particle cut
  Configurable<std::vector<float>> cfPtCut{"cfPtCut", {0.2, 5.0}, "Pt range: {min, max}[GeV], with convention: min < Pt < max"};
  Configurable<std::vector<float>> cfEtaCut{"cfEtaCut", {-0.8, 0.8}, "Eta range: {min, max}, with convention: min < Eta < max"};
  Configurable<std::vector<int>> cfSignCut{"cfSignCut", {1, 0, 1}, "sign of charge, 1 to keep and 0 to discard, {negative, neutral, positive}"};
  Configurable<std::vector<float>> cfTpcNClsFoundCut{"cfTpcNClsFoundCut", {70., 160.}, "range of found TPC clusters for this track geometry: {min, max}"};
  Configurable<std::vector<float>> cfDCAXYCut{"cfDCAXYCut", {-3.2, 3.2}, "range of distance-of-closest-approach (DCA) of the extrapolated track to the primary position in XY-direction: {min, max}[cm]"};
  Configurable<std::vector<float>> cfDCAZCut{"cfDCAZCut", {-2.4, 2.4}, "range of distance-of-closest-approach (DCA) of the extrapolated track to the primary position in Z-direction: {min, max}[cm]"};
  Configurable<float> cfEtaGap{"cfEtaGap", 1., "|dEta| > gap"};

  // *) Others
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/scratch3/go52dab/O2tutorial/tutorial3-6/weights.root", "path to external ROOT file which holds all particle weights in O2 format"};

  // *) Bins
  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {1000, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfCentBins{"cfCentBins", {100, 0., 100.}, "nCenBins, cenMin, cenMax"};
  Configurable<std::vector<float>> cfFT0CMultBins{"cfFT0CMultBins", {100, 0., 60000.}, "nMultBins, MultMin, MultMax"};
  Configurable<std::vector<float>> cfFT0MMultBins{"cfFT0MMultBins", {100, 0., 250000.}, "nMultBins, MultMin, MultMax"};
  Configurable<std::vector<float>> cfFV0AMultBins{"cfFV0AMultBins", {100, 0., 200000.}, "nMultBins, MultMin, MultMax"};
  Configurable<std::vector<float>> cfVerXBins{"cfVerXBins", {100, -0.05, 0.05}, "nVerXBins, VerXMin, VerXMax"};
  Configurable<std::vector<float>> cfVerYBins{"cfVerYBins", {100, -0.05, 0.05}, "nVerYBins, VerYMin, VerYMax"};
  Configurable<std::vector<float>> cfVerZBins{"cfVerZBins", {100, -50., 50.}, "nVerZBins, VerZMin, VerZMax"};
  Configurable<std::vector<float>> cfNumContribBins{"cfNumContribBins", {100, 0., 5000.}, "nNumContribBins, NumContribMin, NumContribMax"};
  Configurable<std::vector<double>> cfTwoParticleCorrBins{"cfTwoParticleCorrBins", {100, 0., 2e-3}, "nTwoParticleCorrBins, TwoParticleCorrMin, TwoParticleCorrMax"};
  Configurable<std::vector<double>> cfFourParticleCorrBins{"cfFourParticleCorrBins", {100, 0., 1e-7}, "nFourParticleCorrBins, FourParticleCorrMin, FourParticleCorrMax"};
  Configurable<std::vector<double>> cfSixParticleCorrBins{"cfSixParticleCorrBins", {100, 0., 1e-11}, "nSixParticleCorrBins, SixParticleCorrMin, SixParticleCorrMax"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    std::array<bool, eProcess_N> fProcess = {{false}}; // Set what to process. See enum EnProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;

    std::string fCentEstm = "FT0M";
    std::string fMultEstm = "FV0A";

    bool fPrintSwitch = true;

    bool fGeneralCutSwitch = true;
    bool fVertexZCutSwitch = true;
    bool fSel8CutSwitch = true;
    bool fCentCutSwitch = true;
    bool fNumContribCutSwitch = true;
    bool fCentCorrCutSwitch = true;
    bool fMultCorrCutSwitch = true;
    bool fNumContribMultCorrCutSwitch = true;

    bool fPtCutSwitch = true;
    bool fEtaCutSwitch = true;
    bool fSignCutSwitch = true;
    bool fTpcNClsFoundCutSwitch = true;
    bool fDCAXYCutSwitch = true;
    bool fDCAZCutSwitch = true;
    bool fEtaGapSwitch = true;

    std::vector<float> fVertexZCut = {-10., 10.};
    std::vector<float> fCentCut = {10., 20.};
    std::vector<float> fNumContribCut = {0, 3000.};
    std::vector<std::vector<std::vector<float>>> fCentAllCorrCut = {{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}, {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}, {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};
    std::vector<std::vector<std::vector<float>>> fMultAllCorrCut = {{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}, {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}, {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};
    std::vector<float> fNumContribMultCorrCut = {0.024, -600., 0.035, 100.};
    std::vector<float> fPtCut = {0.2, 5.0};
    std::vector<float> fEtaCut = {-0.8, 0.8};
    std::vector<int> fSignCut = {1, 0, 1};
    std::vector<float> fTpcNClsFoundCut = {70., 160.};
    std::vector<float> fDCAXYCut = {-3.2, 3.2};
    std::vector<float> fDCAZCut = {-2.4, 2.4};
    float fEtaGap = 1.;

    std::vector<float> fPtBins = {1000, 0., 100.};
    std::vector<float> fPhiBins = {1000, 0., o2::constants::math::TwoPI};

    std::vector<float> fCentBins = {100, 0., 100.};
    std::vector<float> fMultBins = {100, 0, 5000};
    std::vector<float> fFT0CMultBins = {100, 0., 60000.};
    std::vector<float> fFT0MMultBins = {100, 0., 250000.};
    std::vector<float> fFV0AMultBins = {100, 0., 200000.};

    std::vector<float> fVerXBins = {100, -0.05, 0.05};
    std::vector<float> fVerYBins = {100, -0.05, 0.05};
    std::vector<float> fVerZBins = {100, -50., 50.};
    std::vector<float> fNumContribBins = {100, 0, 5000};
    std::vector<double> fTwoParticleCorrBins = {100, 0, 2e-3};
    std::vector<double> fFourParticleCorrBins = {100, 0, 1e-7};
    std::vector<double> fSixParticleCorrBins = {100, 0, 1e-11};

    std::string fFileWithWeights = "/scratch3/go52dab/O2tutorial/analysis_code/weights.root";

  } tc;

  struct ParticleHistograms {
    TList* fParticleHistogramsList = nullptr;
    std::array<std::array<std::array<TH1F*, 2>, eBeforeAfter_N>, eParticleHistograms_N> fParticleHistograms{}; // [type][before/after cut][rec/sim]
  } pc;

  struct EventHistograms {
    TList* fEventHistogramsList = nullptr;
    std::array<std::array<std::array<TH1F*, 2>, eBeforeAfter_N>, eEventHistograms_N> fEventHistograms{}; // [type][before/after cut][rec/sim]
  } ev;

  struct QAHistograms {
    bool fQASwitch = kTRUE;
    TList* fQAHistogramsList = nullptr;
    std::array<std::array<TH2F*, 2>, eQAHistograms_N> fQAHistograms{}; // [type][before/after cut]
  } qa;

  struct CorrHistograms {
    TList* fCorrHistogramsList = nullptr;
    std::array<std::array<std::array<std::array<TH2F*, eBeforeAfter_N>, eMultEstm_N>, eMultEstm_N>, eCorrHistograms_N - 1> fCorrHistograms{}; // [mult/cent][type][type][before/after cut]
    std::array<std::array<TF1*, 2>, eCorrBound_N> fCorrBounds{};                                                                              // [centAB/centAC/centBC/multAB/multAC/multBC/NumContrib vs. mult][upper/lower]
  } cr;

  struct WeightHistograms {
    bool fWeightSwitch = kTRUE;
    TList* fWeightHistogramsList = nullptr;
    // Fill phi/pt histograms with that run number:
    std::map<int, TH1F*> fPtRealByRunMap; // MC rec (too lazy to change name) data pt histograms, valid if processMonteCarlo
    std::map<int, TH1F*> fPtMCByRunMap;   // MC sim (too lazy to change name) data pt histograms, valid if processMonteCarlo
    std::map<int, TH1F*> fPhiByRunMap;    // Phi histograms, valid if processRealData
    // Make weight histograms locally. Upload weight histograms to CCDB:
    std::vector<TH1F*> fWeightHistograms;         // Get all weight histograms with that run number
    std::map<int, TH1F*> fPhiWeightHistogramsMap; // Get phi weight histograms
    std::map<int, TH1F*> fPtWeightHistogramsMap;  // Get pt weight histograms
    // Null weight histograms. Use them when no weight histograms found in the given run number:
    TH1F* fDummyPhiWeightHistogram = nullptr;
    TH1F* fDummyPtWeightHistogram = nullptr;
  } wt;

  struct MultiparticleCorrelationCalculation {
    int h1 = 0;
    int h2 = 0;
    int h3 = 0;
    int h4 = 0;
    int h5 = 0;
    int h6 = 0;
    int h7 = 0;
    int h8 = 0;
    // Book Q-vector components:
    static constexpr int MaxCorrelator = 6;                                                     // <<m>>, in this case 6
    static constexpr int MaxHarmonic = ((MaxCorrelator * (MaxCorrelator + 1)) / 2 - 1) * 2 + 1; // (2 + 3 + ... + n) * 2 + 1, in this case need (2 + 3 + 4) * 2 + 1 = 19 at least for SC(2,3,4)
    static constexpr int MaxPower = MaxCorrelator + 1;

    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorBefore;
    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorAfter;
    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorBeforeA; // Q-vector with eta gap
    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorBeforeB; // Q-vector with eta gap
    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorAfterA;  // Q-vector with eta gap
    std::array<std::array<TComplex, MaxPower>, MaxHarmonic> fQvectorAfterB;  // Q-vector with eta gap
  } mcc;

  struct MultiparticleCorrelationProfile {
    TList* fMultiparticleCorrelationProfilesList = nullptr;
    TList* fMultiparticleCorrelationBackupList = nullptr;
    std::map<int, TList*> fMultiparticleCorrelationByRunMap;
    std::array<TProfile*, eBeforeAfter_N> fTwoParticleCorrelationProfiles{}; // [cut]
    std::array<TProfile*, eBeforeAfter_N> fFourParticleCorrelationProfiles{};
    std::array<TProfile*, eBeforeAfter_N> fSixParticleCorrelationProfiles{};
    std::array<TProfile*, eBeforeAfter_N> fTwoParticleCorrelationGapProfiles{};
    std::array<std::array<std::array<TH1D*, NumTwoPCorrBins>, eBeforeAfter_N>, eBeforeAfter_N> fTwoParticleCorrelationHistograms{}; // [cut][event weight][n]
    std::array<std::array<std::array<TH1D*, NumFourPCorrBins>, eBeforeAfter_N>, eBeforeAfter_N> fFourParticleCorrelationHistograms{};
    std::array<std::array<std::array<TH1D*, NumSixPCorrBins>, eBeforeAfter_N>, eBeforeAfter_N> fSixParticleCorrelationHistograms{};
    std::array<TH1D*, NumTwoPCorrBins> fTwoParticleCorrelationGapHistograms{};
  } mc;

  struct EventByEventQuantities {
    int fRunNumber = 0;
    float fReferenceMultiplicity = 0.;
    std::array<float, eMultEstm_N> fMultiplicityMaxEbye{};
    float fCentrality = 0.;
    float fCentralitySim = 0.;
    float fImpactParameter = 0.;
    float fNumContrib = 0.;
    std::array<std::array<float, NumTwoPCorrBins>, eBeforeAfter_N> fTwoParticleCorrelationEbye{}; // <<2>> [before, after][<v2^2>, <v3^2>, <v4^2>]
    std::array<std::array<float, NumTwoPCorrBins>, eBeforeAfter_N> fTwoParticleCorrelationMinEbye{};
    std::array<std::array<float, NumTwoPCorrBins>, eBeforeAfter_N> fTwoParticleCorrelationMaxEbye{};
    std::array<std::array<float, NumFourPCorrBins>, eBeforeAfter_N> fFourParticleCorrelationEbye{}; // <<4>> [before, after][<v2^2 v2^2>, <v2^2 v3^2>, <v2^2 v4^2>, <v3^2 v3^2>, <v3^2 v4^2>, <v4^2 v4^2>]
    std::array<std::array<float, NumFourPCorrBins>, eBeforeAfter_N> fFourParticleCorrelationMinEbye{};
    std::array<std::array<float, NumFourPCorrBins>, eBeforeAfter_N> fFourParticleCorrelationMaxEbye{};
    std::array<std::array<float, NumSixPCorrBins>, eBeforeAfter_N> fSixParticleCorrelationEbye{}; // <<6>> [before, after][<v2^2 v3^2 v4^2>]
    std::array<std::array<float, NumSixPCorrBins>, eBeforeAfter_N> fSixParticleCorrelationMinEbye{};
    std::array<std::array<float, NumSixPCorrBins>, eBeforeAfter_N> fSixParticleCorrelationMaxEbye{};
  } ebye;

  template <EnRecSim rs, typename T1, typename T2, typename T3, typename T4> // rm should be rec and sim
  bool ctEventCuts(T1 const& collision, T2 const& rlCollisionCentAll, T3 const& rlCollisionMultAll, T4 const& rlCollisionNumContrib)
  {
    bool pass = true;

    bool bGeneralCut = true;
    bool bVertexZCut = true;
    bool bSel8Cut = true;
    bool bCentCut = true;
    bool bNumContribCut = true;
    bool bCentCorrCut = true;
    bool bMultCorrCut = true;
    bool bNumContribMultCorrCut = true;

    // *) For rec event and sim event
    bVertexZCut = collision.posZ() < tc.fVertexZCut[1] &&
                  collision.posZ() > tc.fVertexZCut[0];

    // *) For rec event only
    if constexpr (rs == eRec) {
      // *) General cuts
      bGeneralCut = collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) &&
                    collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard) &&
                    collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup) &&
                    collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC) &&
                    collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) &&
                    collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) &&
                    collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof);
      // *) Sel8Cut
      bSel8Cut = collision.sel8();
      // *) CentCut
      bCentCut = ebye.fCentrality < tc.fCentCut[1] &&
                 ebye.fCentrality > tc.fCentCut[0];
      // *) NumContribCut
      bNumContribCut = ebye.fNumContrib < tc.fNumContribCut[1] &&
                       ebye.fNumContrib > tc.fNumContribCut[0];
      // *) CentCorrCut
      float iCent = 0.;
      float jCent = 0.;
      int nCentComb = 0;
      for (int i = 0; i < eCentEstm_N; i++) {
        iCent = rlCollisionCentAll[i];
        for (int j = i + 1; j < eCentEstm_N; j++) {
          jCent = rlCollisionCentAll[j];
          float upper = cr.fCorrBounds[eCentABBound + nCentComb][0]->Eval(iCent);
          float lower = cr.fCorrBounds[eCentABBound + nCentComb][1]->Eval(iCent);
          bCentCorrCut &= jCent >= lower && jCent <= upper;
          nCentComb += 1;
        }
      }
      // *) MultCorrCut
      float iMult = 0.;
      float jMult = 0.;
      int nMultComb = 0;
      for (int i = 0; i < eMultEstm_N; i++) {
        iMult = rlCollisionMultAll[i];
        for (int j = i + 1; j < eMultEstm_N; j++) {
          jMult = rlCollisionMultAll[j];
          float upper = cr.fCorrBounds[eMultABBound + nMultComb][0]->Eval(iMult);
          float lower = cr.fCorrBounds[eMultABBound + nMultComb][1]->Eval(iMult);
          bMultCorrCut &= jMult >= lower && jMult <= upper;
          nMultComb += 1;
        }
      }
      // *) NumContribMultCorrCut
      float upper = cr.fCorrBounds[eNumContribMultBound][0]->Eval(ebye.fReferenceMultiplicity);
      float lower = cr.fCorrBounds[eNumContribMultBound][1]->Eval(ebye.fReferenceMultiplicity);
      bNumContribMultCorrCut = rlCollisionNumContrib >= lower && rlCollisionNumContrib <= upper;
    }

    // *) For sim event only
    if constexpr (rs == eSim) {
      // *) CentCut
      bCentCut = ebye.fCentralitySim < tc.fCentCut[1] &&
                 ebye.fCentralitySim > tc.fCentCut[0];
    }

    // *) Combine all switches
    if (tc.fGeneralCutSwitch) {
      pass &= bGeneralCut;
    }
    if (tc.fVertexZCutSwitch) {
      pass &= bVertexZCut;
    }
    if (tc.fSel8CutSwitch) {
      pass &= bSel8Cut;
    }
    if (tc.fCentCutSwitch) {
      pass &= bCentCut;
    }
    if (tc.fNumContribCutSwitch) {
      pass &= bNumContribCut;
    }
    if (tc.fCentCorrCutSwitch) {
      pass &= bCentCorrCut;
    }
    if (tc.fMultCorrCutSwitch) {
      pass &= bMultCorrCut;
    }
    if (tc.fNumContribMultCorrCutSwitch) {
      pass &= bNumContribMultCorrCut;
    }

    return pass;
  }

  template <EnRecSim rs, typename T1>
  bool ctParticleCuts(T1 const& track)
  {
    bool pass = true;

    bool bPtCut = true;
    bool bEtaCut = true;
    bool bSignCut = true;
    bool bTpcNClsFoundCut = true;
    bool bDCAXYCut = true;
    bool bDCAZCut = true;

    // *) For rec event and sim event
    bPtCut = track.pt() < tc.fPtCut[1] && track.pt() > tc.fPtCut[0];
    bEtaCut = track.eta() < tc.fEtaCut[1] && track.eta() > tc.fEtaCut[0];

    // *) For rec event only
    if constexpr (rs == eRec) {
      bSignCut = (track.sign() == -1 && tc.fSignCut[0]) ||
                 (track.sign() == 0 && tc.fSignCut[1]) ||
                 (track.sign() == 1 && tc.fSignCut[2]);
      bTpcNClsFoundCut = track.tpcNClsFound() < tc.fTpcNClsFoundCut[1] &&
                         track.tpcNClsFound() > tc.fTpcNClsFoundCut[0];
      bDCAXYCut = track.dcaXY() < tc.fDCAXYCut[1] &&
                  track.dcaXY() > tc.fDCAXYCut[0];
      bDCAZCut = track.dcaZ() < tc.fDCAZCut[1] &&
                 track.dcaZ() > tc.fDCAZCut[0];
    }

    // *) For sim event only
    if constexpr (rs == eSim) {
      TParticlePDG* particle = pdg->GetParticle(track.pdgCode());
      if (!particle) {
        // LOGF(warning, "PDG code %d not found", track.pdgCode());
        bSignCut = false;
      } else {
        // LOGF(info, "PDG code %d found", track.pdgCode());
        float charge = particle->Charge();
        bSignCut = (charge < 0 && tc.fSignCut[0]) ||
                   (charge == 0 && tc.fSignCut[1]) ||
                   (charge > 0 && tc.fSignCut[2]);
      }
    }

    if (tc.fPtCutSwitch) {
      pass &= bPtCut;
    }
    if (tc.fEtaCutSwitch) {
      pass &= bEtaCut;
    }
    if (tc.fSignCutSwitch) {
      pass &= bSignCut;
    }
    if (tc.fTpcNClsFoundCutSwitch) {
      pass &= bTpcNClsFoundCut;
    }
    if (tc.fDCAXYCutSwitch) {
      pass &= bDCAXYCut;
    }
    if (tc.fDCAZCutSwitch) {
      pass &= bDCAZCut;
    }

    return pass;
  }

  TComplex mccQ(int n, int p, EnBeforeAfter eba)
  {
    // Using the fact that Q{-n,p} = Q{n,p}^*.
    if (eba == eBefore) {
      // Q-vector before cut:
      if (n >= 0) {
        return mcc.fQvectorBefore[n][p];
      }
      return TComplex::Conjugate(mcc.fQvectorBefore[-n][p]);
    }
    if (n >= 0) {
      return mcc.fQvectorAfterB[n][p];
    }
    return TComplex::Conjugate(mcc.fQvectorAfterB[-n][p]);
  }

  TComplex mccTwo(int n1, int n2, EnBeforeAfter eba)
  {
    return mccQ(n1, 1, eba) * mccQ(n2, 1, eba) - mccQ(n1 + n2, 2, eba);
  }

  template <std::size_t N>
  TComplex mccRecursion(int n, std::array<int, N> harmonic, EnBeforeAfter eba, int mult = 1, int skip = 0)
  {
    // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk).

    int nm1 = n - 1;
    TComplex c(mccQ(harmonic[nm1], mult, eba));
    if (nm1 == 0) {
      return c;
    }
    c *= mccRecursion(nm1, harmonic, eba);
    if (nm1 == skip) {
      return c;
    }

    int multp1 = mult + 1;
    int nm2 = n - 2;
    int counter1 = 0;
    int hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    TComplex c2(mccRecursion(nm1, harmonic, eba, multp1, nm2));
    int counter2 = n - 3;
    while (counter2 >= skip) {
      harmonic[nm2] = harmonic[counter1];
      harmonic[counter1] = hhold;
      ++counter1;
      hhold = harmonic[counter1];
      harmonic[counter1] = harmonic[nm2];
      harmonic[nm2] = hhold + harmonic[nm1];
      c2 += mccRecursion(nm1, harmonic, eba, multp1, counter2);
      --counter2;
    }
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;

    if (mult == 1) {
      return c - c2;
    }
    return c - static_cast<double>(mult) * c2;
  }

  TObject* getObjectFromList(TList* list, const char* objectName)
  {
    // Get TObject pointer from TList, even if it's in some nested TList.
    // Foreseen to be used to fetch histograms or profiles from files directly.
    // Some ideas taken from TCollection::ls()
    // If you have added histograms directly to files (without TList's), then you can fetch them directly with file->Get("hist-name").

    // Usage: TH1D* hist = (TH1D*)
    // getObjectFromList("some-valid-TList-pointer","some-object-name");

    // Insanity checks:
    if (!list) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (!objectName) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (0 == list->GetEntries()) {
      return nullptr;
    }

    // The object is in the current base list:
    TObject* objectFinal = list->FindObject(objectName); // the final object I am after
    if (objectFinal) {
      return objectFinal;
    }

    // Otherwise, search for the object recursively in the nested lists:
    TObject* objectIter = nullptr; // iterator object in the loop below
    TIter next(list);
    while ((objectIter = next())) // double round braces are to silence the warnings
    {
      if (TString(objectIter->ClassName()).EqualTo("TList")) {
        objectFinal = getObjectFromList(dynamic_cast<TList*>(objectIter), objectName);
        if (objectFinal) {
          return objectFinal;
        }
      }
    } // while(objectIter = next())

    return nullptr;

  } // TObject* getObjectFromList(TList* list, char* objectName)

  std::vector<TH1F*> getHistogramsWithWeights(const char* filePath, const char* runNumber)
  {
    // a) Return value:
    std::vector<TH1F*> histograms;
    TList* baseList = nullptr;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = nullptr; // nested list with run-wise TList's holding run-specific weights

    // c) Determine from filePath if the file in on a local machine, or in home dir AliEn, or in CCDB:
    //    Algorithm:
    //    If filePath begins with "/alice/data/CCDB/" then it's in home dir AliEn.
    //    If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB.
    //    Therefore, files in AliEn and CCDB must be specified with abs path, for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    TString filePathStr(filePath);

    if (filePathStr(0, 7) == "/alice/") {
      bFileIsInAliEn = true;
    } else if (filePathStr(0, 20) == "/alice-ccdb.cern.ch/") {
      bFileIsInCCDB = true;
    }

    if (bFileIsInAliEn) {
      // File you want to access is in your home dir in AliEn:
      const TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", ""); // do not forget to add #include <TGrid.h> to the preamble of your analysis task
      if (!alien) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ"); // yes, ROOT can open a file transparently, even if it's sitting in AliEn, with this specific syntax
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      weightsFile->GetObject("ccdb_object", baseList);
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      // Finally, from the top-level TList, get the desired nested TList => the technical problem here is that it can be nested at any level,
      // for thare there is a helper utility function GetObjectFromList(...) , see its implementation further below
      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(warning, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current runnumber = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
          histograms = {nullptr};
          return histograms;
        }
      }
    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      // My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/ => adapt for your case
      ccdb->setURL("https://alice-ccdb.cern.ch"); // to be able to use "ccdb" this object in your analysis task, see 4b/ below
      baseList = dynamic_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(warning, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current runnumber = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
          histograms = {nullptr};
          return histograms;
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)
    } else {
      // this is the local case, please handle this one now:
      // Check if the external ROOT file exists at specified path:

      if (gSystem->AccessPathName(filePath, kFileExists)) {
        LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath, kFileExists)), filePath = %s \033[0m", filePath);
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      TFile* weightsFile = TFile::Open(filePath, "READ");
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      weightsFile->GetObject("ccdb_object", baseList);

      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(warning, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current runnumber = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
          histograms = {nullptr};
          return histograms;
        }
      }
    }

    TIter next(listWithRuns);
    TObject* object = nullptr;

    while (true) {

      object = next();
      if (!object) {
        break;
      }

      auto* hist = dynamic_cast<TH1F*>(object);
      if (!hist) {
        continue;
      }

      hist->SetDirectory(nullptr);
      auto* histClone = dynamic_cast<TH1F*>(hist->Clone());
      if (!histClone) {
        LOGF(warning, "Failed to clone histogram %s", hist->GetName());
        histograms = {nullptr};
        return histograms;
      }
      histClone->SetTitle(Form("%s:%s", filePath, histClone->GetName()));
      histograms.push_back(histClone);
    }

    return histograms;
  }

  //* ) Define all member functions to be called in the main process* functions:
  template <EnRecSim rs, typename T1, typename T2>
  void runLoop(T1 const& collision, T2 const& tracks)
  {

    // Dry run:
    if (tc.fDryRun) {
      return;
    }

    // Book Q-vector arrays:
    for (int h = 0; h < mcc.MaxHarmonic; h++) {
      for (int p = 0; p < mcc.MaxPower; p++) {
        mcc.fQvectorBefore[h][p] = TComplex(0., 0.);
        mcc.fQvectorAfter[h][p] = TComplex(0., 0.);
        mcc.fQvectorBeforeA[h][p] = TComplex(0., 0.);
        mcc.fQvectorBeforeB[h][p] = TComplex(0., 0.);
        mcc.fQvectorAfterA[h][p] = TComplex(0., 0.);
        mcc.fQvectorAfterB[h][p] = TComplex(0., 0.);
      }
    }

    // Get run number:
    ebye.fRunNumber = collision.bc().runNumber();
    std::string stringRunNumber = std::to_string(ebye.fRunNumber);

    // Book phi histogram with this run number:
    if (!wt.fPhiByRunMap.contains(ebye.fRunNumber)) {
      auto* hPhi = new TH1F(Form("hPhi_run%d", ebye.fRunNumber), Form("phi distribution for run %d", ebye.fRunNumber), static_cast<int>(tc.fPhiBins[0]), tc.fPhiBins[1], tc.fPhiBins[2]);
      hPhi->SetDirectory(nullptr);
      wt.fPhiByRunMap.try_emplace(ebye.fRunNumber, hPhi);
      wt.fWeightHistogramsList->Add(wt.fPhiByRunMap[ebye.fRunNumber]);
    }

    // Book pt MC rec histogram with this run number:
    if (!wt.fPtRealByRunMap.contains(ebye.fRunNumber)) {
      auto* hPtReal = new TH1F(Form("hPtReal_run%d", ebye.fRunNumber), Form("pt MC rec distribution for run %d", ebye.fRunNumber), static_cast<int>(tc.fPtBins[0]), tc.fPtBins[1], tc.fPtBins[2]);
      hPtReal->SetDirectory(nullptr);
      wt.fPtRealByRunMap.try_emplace(ebye.fRunNumber, hPtReal);
      wt.fWeightHistogramsList->Add(wt.fPtRealByRunMap[ebye.fRunNumber]);
    }

    // Book pt MC sim histogram with this run number:
    if (!wt.fPtMCByRunMap.contains(ebye.fRunNumber)) {
      auto* hPtMC = new TH1F(Form("hPtMC_run%d", ebye.fRunNumber), Form("pt MC sim distribution for run %d", ebye.fRunNumber), static_cast<int>(tc.fPtBins[0]), tc.fPtBins[1], tc.fPtBins[2]);
      hPtMC->SetDirectory(nullptr);
      wt.fPtMCByRunMap.try_emplace(ebye.fRunNumber, hPtMC);
      wt.fWeightHistogramsList->Add(wt.fPtMCByRunMap[ebye.fRunNumber]);
    }

    // Get phi and pt weight histogram with this run number:
    if (wt.fWeightSwitch && !wt.fPhiWeightHistogramsMap.contains(ebye.fRunNumber)) {

      TH1F* phiWeightHist = dynamic_cast<TH1F*>(wt.fDummyPhiWeightHistogram->Clone(Form("wPhi_run%d", ebye.fRunNumber)));
      TH1F* ptWeightHist = dynamic_cast<TH1F*>(wt.fDummyPtWeightHistogram->Clone(Form("wPt_run%d", ebye.fRunNumber)));

      wt.fWeightHistograms = getHistogramsWithWeights(tc.fFileWithWeights.c_str(), stringRunNumber.c_str());

      for (auto const& hist : wt.fWeightHistograms) {
        if (!hist) {
          LOGF(warning, "Fail to loop weight histograms");
          continue;
        }
        TString histName = hist->GetName();
        if (histName.BeginsWith("wPhi")) {
          delete phiWeightHist;
          phiWeightHist = dynamic_cast<TH1F*>(hist->Clone(Form("wPhi_run%d", ebye.fRunNumber)));
          break;
        }
      }

      for (auto const& hist : wt.fWeightHistograms) {
        if (!hist) {
          LOGF(warning, "Fail to loop weight histograms");
          continue;
        }
        TString histName = hist->GetName();
        if (histName.BeginsWith("wPt")) {
          delete ptWeightHist;
          ptWeightHist = dynamic_cast<TH1F*>(hist->Clone(Form("wPt_run%d", ebye.fRunNumber)));
          break;
        }
      }

      phiWeightHist->SetDirectory(nullptr);
      wt.fPhiWeightHistogramsMap.try_emplace(ebye.fRunNumber, phiWeightHist);
      wt.fWeightHistogramsList->Add(wt.fPhiWeightHistogramsMap[ebye.fRunNumber]);

      ptWeightHist->SetDirectory(nullptr);
      wt.fPtWeightHistogramsMap.try_emplace(ebye.fRunNumber, ptWeightHist);
      wt.fWeightHistogramsList->Add(wt.fPtWeightHistogramsMap[ebye.fRunNumber]);
    }

    // Multiparticle correlation list with this run number:
    if (!mc.fMultiparticleCorrelationByRunMap.contains(ebye.fRunNumber)) {
      mc.fMultiparticleCorrelationByRunMap.try_emplace(ebye.fRunNumber, new TList());
      mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->SetName(Form("mcc_run%d", ebye.fRunNumber));
      mc.fMultiparticleCorrelationProfilesList->Add(mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]);

      mc.fMultiparticleCorrelationBackupList = nullptr;
      mc.fMultiparticleCorrelationBackupList = new TList();
      mc.fMultiparticleCorrelationBackupList->SetName("backup");
      mc.fMultiparticleCorrelationBackupList->SetOwner(kTRUE);
      mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->Add(mc.fMultiparticleCorrelationBackupList);

      for (int i = 0; i < eBeforeAfter_N; i++) { // before/after cut
        mc.fTwoParticleCorrelationProfiles[i] = nullptr;
        mc.fFourParticleCorrelationProfiles[i] = nullptr;
        mc.fSixParticleCorrelationProfiles[i] = nullptr;
        mc.fTwoParticleCorrelationGapProfiles[i] = nullptr;

        for (int j = 0; j < eBeforeAfter_N; j++) { // before/after weight
          mc.fTwoParticleCorrelationHistograms[i][j].fill(nullptr);
          mc.fFourParticleCorrelationHistograms[i][j].fill(nullptr);
          mc.fSixParticleCorrelationHistograms[i][j].fill(nullptr);
        }
      }

      // Define 2p profiles and histograms:
      for (int i = 0; i < eBeforeAfter_N; i++) { // before/after cut
        mc.fTwoParticleCorrelationProfiles[i] = new TProfile(Form("prof2%sCut", BeforeAfterNames[i + 2]), Form("2-p correlation %s cut", BeforeAfterNames[i]), NumTwoPCorrBins, 0., 1.);
        mc.fTwoParticleCorrelationProfiles[i]->Sumw2();
        mc.fTwoParticleCorrelationGapProfiles[i] = new TProfile(Form("prof2%sCutGapped", BeforeAfterNames[i + 2]), Form("2-p correlation %s cut, gapped", BeforeAfterNames[i]), NumTwoPCorrBins, 0., 1.);
        mc.fTwoParticleCorrelationProfiles[i]->Sumw2();
        mc.fTwoParticleCorrelationGapProfiles[i]->Sumw2();

        for (int k = 0; k < NumTwoPCorrBins; k++) { // v2, v3, v4
          mc.fTwoParticleCorrelationProfiles[i]->GetXaxis()->SetBinLabel(k + 1, Form("<v%d^2>", k + 2));
          mc.fTwoParticleCorrelationGapProfiles[i]->GetXaxis()->SetBinLabel(k + 1, Form("<v%d^2>", k + 2));
          if (static_cast<bool>(i)) {
            mc.fTwoParticleCorrelationGapHistograms[k] = new TH1D(Form("hist2v%dAfterCutWithGap", k + 2), Form("2-p correlation v%d^2 after cut with gap", k + 2), static_cast<int>(tc.fTwoParticleCorrBins[0]), tc.fTwoParticleCorrBins[1], tc.fTwoParticleCorrBins[2]);
            mc.fTwoParticleCorrelationGapHistograms[k]->Sumw2();
            mc.fTwoParticleCorrelationGapHistograms[k]->SetOption("HIST");
          }
          for (int j = 0; j < eBeforeAfter_N; j++) { // before/after weight
            mc.fTwoParticleCorrelationHistograms[i][j][k] = new TH1D(Form("hist2v%d%sCut%sWeight", k + 2, BeforeAfterNames[i + 2], BeforeAfterNames[j + 2]), Form("2-p correlation v%d^2 %s cut %s weight", k + 2, BeforeAfterNames[i], BeforeAfterNames[j]), static_cast<int>(tc.fTwoParticleCorrBins[0]), tc.fTwoParticleCorrBins[1], tc.fTwoParticleCorrBins[2]);
            if (static_cast<bool>(j)) {
              mc.fTwoParticleCorrelationHistograms[i][j][k]->Sumw2();
              mc.fTwoParticleCorrelationHistograms[i][j][k]->SetOption("HIST");
            }
            mc.fMultiparticleCorrelationBackupList->Add(mc.fTwoParticleCorrelationHistograms[i][j][k]);
          }
          if (static_cast<bool>(i)) {
            mc.fMultiparticleCorrelationBackupList->Add(mc.fTwoParticleCorrelationGapHistograms[k]);
          }
        }
        mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->Add(mc.fTwoParticleCorrelationProfiles[i]);
        mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->Add(mc.fTwoParticleCorrelationGapProfiles[i]);
      }

      // Define 4p profiles and histograms:
      for (int i = 0; i < eBeforeAfter_N; i++) { // before/after cut
        mc.fFourParticleCorrelationProfiles[i] = new TProfile(Form("prof4%sCut", BeforeAfterNames[i + 2]), Form("4-p correlation %s cut", BeforeAfterNames[i]), NumFourPCorrBins, 0., 1.);
        mc.fFourParticleCorrelationProfiles[i]->Sumw2();

        for (int k = 0; k < NumFourPCorrBins; k++) { // v2v2, v2v3, v2v4, v3v3, v3v4, v4v4
          mc.fFourParticleCorrelationProfiles[i]->GetXaxis()->SetBinLabel(k + 1, Form("<v%d^2v%d^2>", FourPHarmonicIndex[k][0], FourPHarmonicIndex[k][1]));

          for (int j = 0; j < eBeforeAfter_N; j++) { // before/after weight
            mc.fFourParticleCorrelationHistograms[i][j][k] = new TH1D(Form("hist4v%dv%d%sCut%sWeight", FourPHarmonicIndex[k][0], FourPHarmonicIndex[k][1], BeforeAfterNames[i + 2], BeforeAfterNames[j + 2]), Form("4-p correlation v%d^2v%d^2 %s cut %s weight", FourPHarmonicIndex[k][0], FourPHarmonicIndex[k][1], BeforeAfterNames[i], BeforeAfterNames[j]), static_cast<int>(tc.fFourParticleCorrBins[0]), tc.fFourParticleCorrBins[1], tc.fFourParticleCorrBins[2]);
            if (static_cast<bool>(j)) {
              mc.fFourParticleCorrelationHistograms[i][j][k]->Sumw2();
              mc.fFourParticleCorrelationHistograms[i][j][k]->SetOption("HIST");
            }
            mc.fMultiparticleCorrelationBackupList->Add(mc.fFourParticleCorrelationHistograms[i][j][k]);
          }
        }
        mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->Add(mc.fFourParticleCorrelationProfiles[i]);
      }

      // Define 6p profiles and histograms:
      for (int i = 0; i < eBeforeAfter_N; i++) { // before/after cut
        mc.fSixParticleCorrelationProfiles[i] = new TProfile(Form("prof6%sCut", BeforeAfterNames[i + 2]), Form("6-p correlation %s cut", BeforeAfterNames[i]), NumSixPCorrBins, 0., 1.);
        mc.fSixParticleCorrelationProfiles[i]->Sumw2();

        for (int k = 0; k < NumSixPCorrBins; k++) { // v2v3v4
          mc.fSixParticleCorrelationProfiles[i]->GetXaxis()->SetBinLabel(k + 1, "<v2^2v3^2v4^2>");

          for (int j = 0; j < eBeforeAfter_N; j++) { // before/after weight
            mc.fSixParticleCorrelationHistograms[i][j][k] = new TH1D(Form("hist6v2v3v4%sCut%sWeight", BeforeAfterNames[i + 2], BeforeAfterNames[j + 2]), Form("6-p correlation v2^2v3^2v4^2 %s cut %s weight", BeforeAfterNames[i], BeforeAfterNames[j]), static_cast<int>(tc.fSixParticleCorrBins[0]), tc.fSixParticleCorrBins[1], tc.fSixParticleCorrBins[2]);
            if (static_cast<bool>(j)) {
              mc.fSixParticleCorrelationHistograms[i][j][k]->Sumw2();
              mc.fSixParticleCorrelationHistograms[i][j][k]->SetOption("HIST");
            }
            mc.fMultiparticleCorrelationBackupList->Add(mc.fSixParticleCorrelationHistograms[i][j][k]);
          }
        }
        mc.fMultiparticleCorrelationByRunMap[ebye.fRunNumber]->Add(mc.fSixParticleCorrelationProfiles[i]);
      }
    }

    // Real data centrality:
    std::array<float, eCentEstm_N> rlCollisionCentAll = {
      collision.centFT0C(),
      collision.centFT0M(),
      collision.centFV0A()};

    float rlCollisionCent = 0.;

    for (int i = 0; i < eCentEstm_N; i++) {
      if (tc.fPrintSwitch) {
        LOGF(info, "%s Centrality: %f", CentEstmNames[i],
             rlCollisionCentAll[i]);
      }
      if (tc.fCentEstm == CentEstmNames[i]) {
        rlCollisionCent = rlCollisionCentAll[i];
      }
    }

    // Real data multiplicity:
    std::array<float, eMultEstm_N> rlCollisionMultAll = {
      static_cast<float>(collision.multFT0C()),
      static_cast<float>(collision.multFT0M()),
      static_cast<float>(collision.multFV0A())};
    float rlCollisionMult = 0.;

    for (int i = 0; i < eMultEstm_N; i++) {
      if (tc.fPrintSwitch) {
        LOGF(info, "%s Multiplicity: %f", MultEstmNames[i], rlCollisionMultAll[i]);
      }
      if (tc.fMultEstm == MultEstmNames[i]) {
        rlCollisionMult = rlCollisionMultAll[i];
      }
      if (rlCollisionMultAll[i] >= ebye.fMultiplicityMaxEbye[i]) {
        ebye.fMultiplicityMaxEbye[i] = rlCollisionMultAll[i];
        if (tc.fPrintSwitch) {
          LOGF(info, "max %s multiplicity so far = %f", MultEstmNames[i], ebye.fMultiplicityMaxEbye[i]);
        }
      }
    }

    // Real data nContrib:
    float rlCollisionNumContrib = 0.;
    rlCollisionNumContrib = static_cast<float>(collision.numContrib());

    // Event-by-event quantity:
    ebye.fCentrality = rlCollisionCent;
    ebye.fReferenceMultiplicity = rlCollisionMult;
    ebye.fNumContrib = rlCollisionNumContrib;

    // Print...
    if (tc.fPrintSwitch) {

      LOGF(info, "Run number: %d", ebye.fRunNumber);

      LOGF(info, "Centrality: %f", rlCollisionCent);
      LOGF(info, "Multiplicity: %f", static_cast<float>(rlCollisionMult));

      LOGF(info, "Vertex X position: %f", collision.posX());
      LOGF(info, "Vertex Y position: %f", collision.posY());
      LOGF(info, "Vertex Z position: %f", collision.posZ());

      LOGF(info, "NContributors: %f", static_cast<float>(rlCollisionNumContrib));
    }

    // If Rec or RecAndSim:
    if constexpr (rs == eRec || rs == eRecAndSim) {

      // Fill real event histograms before cut:
      ev.fEventHistograms[eCent][eRec][eBefore]->Fill(rlCollisionCent);
      ev.fEventHistograms[eMult][eRec][eBefore]->Fill(rlCollisionMult);
      ev.fEventHistograms[eVertexX][eRec][eBefore]->Fill(collision.posX());
      ev.fEventHistograms[eVertexY][eRec][eBefore]->Fill(collision.posY());
      ev.fEventHistograms[eVertexZ][eRec][eBefore]->Fill(collision.posZ());
      ev.fEventHistograms[eNumContrib][eRec][eBefore]->Fill(rlCollisionNumContrib);

      // Fill centrality correlation histograms before cut:
      if (tc.fCentCorrCutSwitch) {
        for (int i = 0; i < eCentEstm_N; i++) {
          for (int j = i + 1; j < eCentEstm_N; j++) {
            auto* h = cr.fCorrHistograms[eCorrCent][i][j][eBefore];
            if (!h) {
              LOGF(fatal, "Missing histogram cr.fCorrHistograms[eCorrCent][%d][%d][eBefore]", i, j);
            }
            if (rlCollisionCentAll[i] >= 0. && rlCollisionCentAll[j] >= 0.) {
              h->Fill(rlCollisionCentAll[i], rlCollisionCentAll[j]);
            }
          }
        }
      }

      // Fill multiplicity correlation histograms before cut:
      if (tc.fMultCorrCutSwitch) {
        for (int i = 0; i < eMultEstm_N; i++) {
          for (int j = i + 1; j < eMultEstm_N; j++) {
            auto* h = cr.fCorrHistograms[eCorrMult][i][j][eBefore];
            if (!h) {
              LOGF(fatal, "Missing histogram cr.fCorrHistograms[eCorrMult][%d][%d][eBefore]", i, j);
            }
            if (rlCollisionMultAll[i] >= 0. && rlCollisionMultAll[j] >= 0.) {
              h->Fill(rlCollisionMultAll[i], rlCollisionMultAll[j]);
            }
          }
        }
      }

      // Fill NContrib vs. multiplicity QA histogram before cut:
      if (qa.fQASwitch) {
        qa.fQAHistograms[eQAMultNumContrib][eBefore]->Fill(rlCollisionMult, rlCollisionNumContrib);
      }

      // If RecAndSim:
      if constexpr (rs == eRecAndSim) {

        if (!collision.has_mcCollision()) {
          if (tc.fPrintSwitch) {
            LOGF(warning, "No MC collision for this collision, skip...");
          }
        } else {
          // Define MC collision:
          auto mccollision = collision.mcCollision();

          // Define MC centrality:
          float mcCollisionCent = 0.;
          double mcCollisionCentDouble = 0.;
          double b = mccollision.impactParameter() * 1e-15; // convert fm to m
          double xs = 7.71e-28;                             // convert barn to m^2
          mcCollisionCentDouble = o2::constants::math::PI * b * b / xs * 100;
          mcCollisionCent = static_cast<float>(mcCollisionCentDouble);

          // Event-by-event quantity:
          ebye.fCentralitySim = mcCollisionCent;
          ebye.fImpactParameter = b;

          if (tc.fPrintSwitch) {
            LOGF(info, "mc impact param (fm): %f", mccollision.impactParameter());
            LOGF(info, "mc centrality: %f", mcCollisionCent);
          }

          // Fill MC event histograms before cut:
          ev.fEventHistograms[eCent][eSim][eBefore]->Fill(mcCollisionCent);
          ev.fEventHistograms[eVertexX][eSim][eBefore]->Fill(mccollision.posX());
          ev.fEventHistograms[eVertexY][eSim][eBefore]->Fill(mccollision.posY());
          ev.fEventHistograms[eVertexZ][eSim][eBefore]->Fill(mccollision.posZ());

          // Fill MC event histograms after cut:
          if (ctEventCuts<eSim>(mccollision, rlCollisionCentAll, rlCollisionMultAll, rlCollisionNumContrib)) {
            ev.fEventHistograms[eCent][eSim][eAfter]->Fill(mcCollisionCent);
            ev.fEventHistograms[eVertexX][eSim][eAfter]->Fill(mccollision.posX());
            ev.fEventHistograms[eVertexY][eSim][eAfter]->Fill(mccollision.posY());
            ev.fEventHistograms[eVertexZ][eSim][eAfter]->Fill(mccollision.posZ());
          }

          // Fill centrality QA histogram before/after cut:
          if (qa.fQASwitch) {
            qa.fQAHistograms[eQACent][eBefore]->Fill(rlCollisionCent, mcCollisionCent);
            if (ctEventCuts<eRec>(collision, rlCollisionCentAll, rlCollisionMultAll, rlCollisionNumContrib) &&
                ctEventCuts<eSim>(mccollision, rlCollisionCentAll, rlCollisionMultAll, rlCollisionNumContrib)) {
              qa.fQAHistograms[eQACent][eAfter]->Fill(rlCollisionCent, mcCollisionCent);
            }
          }
        }
      }

      // Fill real event histograms after cut
      if (ctEventCuts<eRec>(collision, rlCollisionCentAll, rlCollisionMultAll, rlCollisionNumContrib)) {
        ev.fEventHistograms[eCent][eRec][eAfter]->Fill(rlCollisionCent);
        ev.fEventHistograms[eMult][eRec][eAfter]->Fill(rlCollisionMult);
        ev.fEventHistograms[eVertexX][eRec][eAfter]->Fill(collision.posX());
        ev.fEventHistograms[eVertexY][eRec][eAfter]->Fill(collision.posY());
        ev.fEventHistograms[eVertexZ][eRec][eAfter]->Fill(collision.posZ());
        ev.fEventHistograms[eNumContrib][eRec][eAfter]->Fill(rlCollisionNumContrib);

        // Fill centrality correlation histograms after cut:
        if (tc.fCentCorrCutSwitch) {
          for (int i = 0; i < eCentEstm_N; i++) {
            for (int j = i + 1; j < eCentEstm_N; j++) {
              auto* h = cr.fCorrHistograms[eCorrCent][i][j][eAfter];
              if (!h) {
                LOGF(fatal, "Missing histogram cr.fCorrHistograms[eCorrCent][%d][%d][eAfter]", i, j);
              }
              h->Fill(rlCollisionCentAll[i], rlCollisionCentAll[j]);
            }
          }
        }

        // Fill multiplicity correlation histograms after cut:
        if (tc.fMultCorrCutSwitch) {
          for (int i = 0; i < eMultEstm_N; i++) {
            for (int j = i + 1; j < eMultEstm_N; j++) {
              auto* h = cr.fCorrHistograms[eCorrMult][i][j][eAfter];
              if (!h) {
                LOGF(fatal, "Missing histogram cr.fCorrHistograms[eCorrMult][%d][%d][eAfter]", i, j);
              }
              h->Fill(rlCollisionMultAll[i], rlCollisionMultAll[j]);
            }
          }
        }

        // Fill NContrib vs. multiplicity QA histogram after cut:
        if (qa.fQASwitch) {
          qa.fQAHistograms[eQAMultNumContrib][eAfter]->Fill(rlCollisionMult, rlCollisionNumContrib);
        }

      } else {
        // Fail the event cut, skip this collision:
        return;
      }
    }

    int nTracksBefore = tracks.size();
    int nTracksAfter = 0;
    int nTracksBeforeA = 0;
    int nTracksBeforeB = 0;
    int nTracksAfterA = 0;
    int nTracksAfterB = 0;

    // Calculate Q-vectors for available angles and weights:
    double dEta = 0.;
    double dPhi = 0.;         // particle angle
    double dPt = 0.;
    double wPhi = 1.;         // particle weight
    double wPt = 1.;
    double wPhiToPowerP = 1.; // particle weight raised to power p. wPhi is actually wPt*wPhi but I'm too lazy to change the name.

    // Main loop over particles:
    for (auto const& track : tracks) {
      // LOGF(info, "Track azimuthal angle: %f", track.phi());
      // LOGF(info, "Transverse momentum: %f", track.pt());

      if constexpr (rs == eRec || rs == eRecAndSim) {

        // Fill phi/pt real histogram with this run number:
        wt.fPhiByRunMap.at(ebye.fRunNumber)->Fill(track.phi());
        wt.fPtRealByRunMap.at(ebye.fRunNumber)->Fill(track.pt());

        // Fill track histograms before cut:
        pc.fParticleHistograms[ePt][eRec][eBefore]->Fill(track.pt());
        pc.fParticleHistograms[ePhi][eRec][eBefore]->Fill(track.phi());

        // Calculating Q-vector before cut:
        dEta = track.eta();
        dPhi = track.phi();
        dPt = track.pt();
        if (wt.fWeightSwitch) {
          if (!wt.fPhiWeightHistogramsMap.contains(ebye.fRunNumber) || !wt.fPhiWeightHistogramsMap.at(ebye.fRunNumber)) {
            LOGF(fatal, "Missing Phi weight histogram for run %d", ebye.fRunNumber);
          }
          if (!wt.fPtWeightHistogramsMap.contains(ebye.fRunNumber) || !wt.fPtWeightHistogramsMap.at(ebye.fRunNumber)) {
            LOGF(fatal, "Missing Pt weight histogram for run %d", ebye.fRunNumber);
          }

          auto* histPhi = wt.fPhiWeightHistogramsMap.at(ebye.fRunNumber);
          auto* histPt = wt.fPtWeightHistogramsMap.at(ebye.fRunNumber);
          wPhi = histPhi->GetBinContent(histPhi->GetXaxis()->FindBin(dPhi));
          wPt = histPt->GetBinContent(histPt->GetXaxis()->FindBin(dPt));

          wPhi *= wPt;
        }

        for (int h = 0; h < mcc.MaxHarmonic; h++) {
          for (int p = 0; p < mcc.MaxPower; p++) {
            if (wt.fWeightSwitch) {
              wPhiToPowerP = std::pow(wPhi, p);
            }
            mcc.fQvectorBefore[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
            if (tc.fEtaGapSwitch) {
              if (dEta < -tc.fEtaGap / 2.) {
                mcc.fQvectorBeforeA[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
              } else if (dEta > tc.fEtaGap / 2.) {
                mcc.fQvectorBeforeB[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
              }
            }
          }
        }

        if (tc.fEtaGapSwitch) {
          if (dEta < -tc.fEtaGap / 2.) {
            nTracksBeforeA += 1;
          } else if (dEta > tc.fEtaGap / 2.) {
            nTracksBeforeB += 1;
          }
        }

        if (ctParticleCuts<eRec>(track)) {

          // Fill particle histograms after cut:
          pc.fParticleHistograms[ePt][eRec][eAfter]->Fill(track.pt());
          pc.fParticleHistograms[ePhi][eRec][eAfter]->Fill(track.phi());

          // Calculating Q-vector after cut:
          for (int h = 0; h < mcc.MaxHarmonic; h++) {
            for (int p = 0; p < mcc.MaxPower; p++) {
              if (wt.fWeightSwitch) {
                wPhiToPowerP = std::pow(wPhi, p);
              }
              mcc.fQvectorAfter[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
              if (tc.fEtaGapSwitch) {
                if (dEta < -tc.fEtaGap / 2.) {
                  mcc.fQvectorAfterA[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
                } else if (dEta > tc.fEtaGap / 2.) {
                  mcc.fQvectorAfterB[h][p] += TComplex(wPhiToPowerP * std::cos(h * dPhi), wPhiToPowerP * std::sin(h * dPhi));
                }
              }
            }
          }

          nTracksAfter += 1;
          if (tc.fEtaGapSwitch) {
            if (dEta < -tc.fEtaGap / 2.) {
              nTracksAfterA += 1;
            } else if (dEta > tc.fEtaGap / 2.) {
              nTracksAfterB += 1;
            }
          }
        }
        // ...

        // ... and corresponding MC truth simulated:
        // See
        // https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
        // See
        // https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
        if constexpr (rs == eRecAndSim) {
          if (!track.has_mcParticle()) {
            if (tc.fPrintSwitch) {
              LOGF(warning, "  No MC particle for this track, skip...");
            }
          } else {
            // Corresponding MC truth simulated particle
            auto mcparticle = track.mcParticle();

            //  Fill pt MC sim histogram with this run number:
            wt.fPtMCByRunMap.at(ebye.fRunNumber)->Fill(mcparticle.pt());

            // Fill MC particle histograms before cut:
            pc.fParticleHistograms[ePt][eSim][eBefore]->Fill(mcparticle.pt());
            pc.fParticleHistograms[ePhi][eSim][eBefore]->Fill(mcparticle.phi());

            if (ctParticleCuts<eSim>(mcparticle)) {
              // Fill MC particle histograms after cut:
              pc.fParticleHistograms[ePt][eSim][eAfter]->Fill(mcparticle.pt());
              pc.fParticleHistograms[ePhi][eSim][eAfter]->Fill(mcparticle.phi());
            }
          }
        } // end of if constexpr (rs == eRecAndSim) {
      } // if constexpr (rs == eRec || rs == eRecAndSim) {
    } // end of for (int64_t i = 0; i < tracks.size(); i++) {

    // 2-p correlations:
    for (int i = 0; i < NumTwoPCorrBins; i++) {
      int k = 2;
      mcc.h1 = -(i + 2);
      mcc.h2 = i + 2;
      std::array<int, 2> harmonicsTwoNum = {mcc.h1, mcc.h2};
      std::array<int, 2> harmonicsTwoDen = {0, 0};

      // Before cut:
      TComplex twoRecursionBefore = mccRecursion(2, harmonicsTwoNum, eBefore) / mccRecursion(2, harmonicsTwoDen, eBefore).Re();
      double wTwoRecursionBefore = mccRecursion(2, harmonicsTwoDen, eBefore).Re();
      ebye.fTwoParticleCorrelationEbye[eBefore][i] = twoRecursionBefore.Re();

      if (nTracksBefore > k - 1 && wTwoRecursionBefore > 0.) {
        // Fill vn^2 profile and histogram:
        mc.fTwoParticleCorrelationProfiles[eBefore]->Fill(mc.fTwoParticleCorrelationProfiles[eBefore]->GetXaxis()->GetBinCenter(i + 1), ebye.fTwoParticleCorrelationEbye[eBefore][i], wTwoRecursionBefore);
        mc.fTwoParticleCorrelationHistograms[eBefore][eBefore][i]->Fill(ebye.fTwoParticleCorrelationEbye[eBefore][i]);
        mc.fTwoParticleCorrelationHistograms[eBefore][eAfter][i]->Fill(ebye.fTwoParticleCorrelationEbye[eBefore][i], wTwoRecursionBefore);

        // Find max and min:
        if (ebye.fTwoParticleCorrelationEbye[eBefore][i] > ebye.fTwoParticleCorrelationMaxEbye[eBefore][i]) {
          ebye.fTwoParticleCorrelationMaxEbye[eBefore][i] = ebye.fTwoParticleCorrelationEbye[eBefore][i];
          LOGF(info, "max v%d^2 so far = %e", i + 2, ebye.fTwoParticleCorrelationEbye[eBefore][i]);
        }
        if (ebye.fTwoParticleCorrelationEbye[eBefore][i] < ebye.fTwoParticleCorrelationMinEbye[eBefore][i]) {
          ebye.fTwoParticleCorrelationMinEbye[eBefore][i] = ebye.fTwoParticleCorrelationEbye[eBefore][i];
          LOGF(info, "min v%d^2 so far = %e", i + 2, ebye.fTwoParticleCorrelationEbye[eBefore][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksBefore = %d, wTwoRecursionBefore = %e", rlCollisionCent, nTracksBefore, wTwoRecursionBefore);
      }

      // After cut:
      TComplex twoRecursionAfter = mccRecursion(2, harmonicsTwoNum, eAfter) / mccRecursion(2, harmonicsTwoDen, eAfter).Re();
      double wTwoRecursionAfter = mccRecursion(2, harmonicsTwoDen, eAfter).Re();
      ebye.fTwoParticleCorrelationEbye[eAfter][i] = twoRecursionAfter.Re();

      if (nTracksAfter > k - 1 && wTwoRecursionAfter > 0.) {
        // Fill vn^2 profile and histogram:
        mc.fTwoParticleCorrelationProfiles[eAfter]->Fill(mc.fTwoParticleCorrelationProfiles[eAfter]->GetXaxis()->GetBinCenter(i + 1), ebye.fTwoParticleCorrelationEbye[eAfter][i], wTwoRecursionAfter);
        mc.fTwoParticleCorrelationHistograms[eAfter][eBefore][i]->Fill(ebye.fTwoParticleCorrelationEbye[eAfter][i]);
        mc.fTwoParticleCorrelationHistograms[eAfter][eAfter][i]->Fill(ebye.fTwoParticleCorrelationEbye[eAfter][i], wTwoRecursionAfter);

        // Find max and min:
        if (ebye.fTwoParticleCorrelationEbye[eAfter][i] > ebye.fTwoParticleCorrelationMaxEbye[eAfter][i]) {
          ebye.fTwoParticleCorrelationMaxEbye[eAfter][i] = ebye.fTwoParticleCorrelationEbye[eAfter][i];
          LOGF(info, "max v%d^2 so far = %e", i + 2, ebye.fTwoParticleCorrelationEbye[eAfter][i]);
        }
        if (ebye.fTwoParticleCorrelationEbye[eAfter][i] < ebye.fTwoParticleCorrelationMinEbye[eAfter][i]) {
          ebye.fTwoParticleCorrelationMinEbye[eAfter][i] = ebye.fTwoParticleCorrelationEbye[eAfter][i];
          LOGF(info, "min v%d^2 so far = %e", i + 2, ebye.fTwoParticleCorrelationEbye[eAfter][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksAfter = %d, wTwoRecursionAfter = %e", rlCollisionCent, nTracksAfter, wTwoRecursionAfter);
      }

      // Before cut, with gap:
      TComplex qvba = mcc.fQvectorBeforeA[mcc.h2][0];
      TComplex qvbb = mcc.fQvectorBeforeB[mcc.h2][0];
      double gapb = 0.;
      if (nTracksBeforeA * nTracksBeforeB > 0) {
        gapb = (qvba * TComplex::Conjugate(qvbb)).Re() / (nTracksBeforeA * nTracksBeforeB);
        mc.fTwoParticleCorrelationGapProfiles[eBefore]->Fill(mc.fTwoParticleCorrelationGapProfiles[eBefore]->GetXaxis()->GetBinCenter(i + 1), gapb, nTracksBeforeA * nTracksBeforeB);
        // mc.fTwoParticleCorrelationGapHistograms[i]->Fill(gapb, nTracksBeforeA * nTracksBeforeB);
      } else {
        LOGF(warning, "cent=%f, nTracksBeforeA = %d, nTracksBeforeB = %d", rlCollisionCent, nTracksBeforeA, nTracksBeforeB);
      }

      // After cut, with gap:
      TComplex qvaa = mcc.fQvectorAfterA[mcc.h2][0];
      TComplex qvab = mcc.fQvectorAfterB[mcc.h2][0];
      double gapa = 0.;
      if (nTracksAfterA * nTracksAfterB > 0) {
        gapa = (qvaa * TComplex::Conjugate(qvab)).Re() / (nTracksAfterA * nTracksAfterB);
        mc.fTwoParticleCorrelationGapProfiles[eAfter]->Fill(mc.fTwoParticleCorrelationGapProfiles[eAfter]->GetXaxis()->GetBinCenter(i + 1), gapa, nTracksAfterA * nTracksAfterB);
        mc.fTwoParticleCorrelationGapHistograms[i]->Fill(gapa, nTracksAfterA * nTracksAfterB);
      } else {
        LOGF(warning, "cent=%f, nTracksAfterA = %d, nTracksAfterB = %d", rlCollisionCent, nTracksAfterA, nTracksAfterB);
      }
    }

    //  4-p correlations:
    for (int i = 0; i < NumFourPCorrBins; i++) {
      int k = 4;
      mcc.h1 = -FourPHarmonicIndex[i][0];
      mcc.h2 = -FourPHarmonicIndex[i][1];
      mcc.h3 = FourPHarmonicIndex[i][0];
      mcc.h4 = FourPHarmonicIndex[i][1];
      std::array<int, 4> harmonicsFourNum = {mcc.h1, mcc.h2, mcc.h3, mcc.h4};
      std::array<int, 4> harmonicsFourDen = {0, 0, 0, 0};

      // Before cut:
      TComplex fourRecursionBefore = mccRecursion(4, harmonicsFourNum, eBefore) / mccRecursion(4, harmonicsFourDen, eBefore).Re();
      double wFourRecursionBefore = mccRecursion(4, harmonicsFourDen, eBefore).Re();
      ebye.fFourParticleCorrelationEbye[eBefore][i] = fourRecursionBefore.Re();

      if (nTracksBefore > k - 1 && wFourRecursionBefore > 0.) {
        // Fill vn^2vm^2 profile and histogram:
        mc.fFourParticleCorrelationProfiles[eBefore]->Fill(mc.fFourParticleCorrelationProfiles[eBefore]->GetXaxis()->GetBinCenter(i + 1), ebye.fFourParticleCorrelationEbye[eBefore][i], wFourRecursionBefore);
        mc.fFourParticleCorrelationHistograms[eBefore][eBefore][i]->Fill(ebye.fFourParticleCorrelationEbye[eBefore][i]);
        mc.fFourParticleCorrelationHistograms[eBefore][eAfter][i]->Fill(ebye.fFourParticleCorrelationEbye[eBefore][i], wFourRecursionBefore);

        // Find max and min:
        if (ebye.fFourParticleCorrelationEbye[eBefore][i] > ebye.fFourParticleCorrelationMaxEbye[eBefore][i]) {
          ebye.fFourParticleCorrelationMaxEbye[eBefore][i] = ebye.fFourParticleCorrelationEbye[eBefore][i];
          LOGF(info, "max v%d^2v%d^2 so far = %e", mcc.h3, mcc.h4, ebye.fFourParticleCorrelationEbye[eBefore][i]);
        }
        if (ebye.fFourParticleCorrelationEbye[eBefore][i] < ebye.fFourParticleCorrelationMinEbye[eBefore][i]) {
          ebye.fFourParticleCorrelationMinEbye[eBefore][i] = ebye.fFourParticleCorrelationEbye[eBefore][i];
          LOGF(info, "min v%d^2v%d^2 so far = %e", mcc.h3, mcc.h4, ebye.fFourParticleCorrelationEbye[eBefore][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksBefore = %d, wFourRecursionBefore = %e", rlCollisionCent, nTracksBefore, wFourRecursionBefore);
      }

      // After cut:
      TComplex fourRecursionAfter = mccRecursion(4, harmonicsFourNum, eAfter) / mccRecursion(4, harmonicsFourDen, eAfter).Re();
      double wFourRecursionAfter = mccRecursion(4, harmonicsFourDen, eAfter).Re();
      ebye.fFourParticleCorrelationEbye[eAfter][i] = fourRecursionAfter.Re();

      if (nTracksAfter > k - 1 && wFourRecursionAfter > 0.) {
        // Fill vn^2vm^2 profile and histogram:
        mc.fFourParticleCorrelationProfiles[eAfter]->Fill(mc.fFourParticleCorrelationProfiles[eAfter]->GetXaxis()->GetBinCenter(i + 1), ebye.fFourParticleCorrelationEbye[eAfter][i], wFourRecursionAfter);
        mc.fFourParticleCorrelationHistograms[eAfter][eBefore][i]->Fill(ebye.fFourParticleCorrelationEbye[eAfter][i]);
        mc.fFourParticleCorrelationHistograms[eAfter][eAfter][i]->Fill(ebye.fFourParticleCorrelationEbye[eAfter][i], wFourRecursionAfter);

        // Find max and min:
        if (ebye.fFourParticleCorrelationEbye[eAfter][i] > ebye.fFourParticleCorrelationMaxEbye[eAfter][i]) {
          ebye.fFourParticleCorrelationMaxEbye[eAfter][i] = ebye.fFourParticleCorrelationEbye[eAfter][i];
          LOGF(info, "max v%d^2v%d^2 so far = %e", mcc.h3, mcc.h4, ebye.fFourParticleCorrelationEbye[eAfter][i]);
        }
        if (ebye.fFourParticleCorrelationEbye[eAfter][i] < ebye.fFourParticleCorrelationMinEbye[eAfter][i]) {
          ebye.fFourParticleCorrelationMinEbye[eAfter][i] = ebye.fFourParticleCorrelationEbye[eAfter][i];
          LOGF(info, "min v%d^2v%d^2 so far = %e", mcc.h3, mcc.h4, ebye.fFourParticleCorrelationEbye[eAfter][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksAfter = %d, wFourRecursionAfter = %e", rlCollisionCent, nTracksAfter, wFourRecursionAfter);
      }
    }

    //  6-p correlations:
    for (int i = 0; i < NumSixPCorrBins; i++) {
      int k = 6;
      mcc.h1 = -2;
      mcc.h2 = -3;
      mcc.h3 = -4;
      mcc.h4 = 2;
      mcc.h5 = 3;
      mcc.h6 = 4;
      std::array<int, 6> harmonicsSixNum = {mcc.h1, mcc.h2, mcc.h3, mcc.h4, mcc.h5, mcc.h6};
      std::array<int, 6> harmonicsSixDen = {0, 0, 0, 0, 0, 0};

      // Before cut:
      TComplex sixRecursionBefore = mccRecursion(6, harmonicsSixNum, eBefore) / mccRecursion(6, harmonicsSixDen, eBefore).Re();
      double wSixRecursionBefore = mccRecursion(6, harmonicsSixDen, eBefore).Re();
      ebye.fSixParticleCorrelationEbye[eBefore][i] = sixRecursionBefore.Re();

      if (nTracksBefore > k - 1 && wSixRecursionBefore > 0.) {
        // Fill v2^2v3^2v4^2 profile and histogram:
        mc.fSixParticleCorrelationProfiles[eBefore]->Fill(mc.fSixParticleCorrelationProfiles[eBefore]->GetXaxis()->GetBinCenter(i + 1), ebye.fSixParticleCorrelationEbye[eBefore][i], wSixRecursionBefore);
        mc.fSixParticleCorrelationHistograms[eBefore][eBefore][i]->Fill(ebye.fSixParticleCorrelationEbye[eBefore][i]);
        mc.fSixParticleCorrelationHistograms[eBefore][eAfter][i]->Fill(ebye.fSixParticleCorrelationEbye[eBefore][i], wSixRecursionBefore);

        // Find max and min:
        if (ebye.fSixParticleCorrelationEbye[eBefore][i] > ebye.fSixParticleCorrelationMaxEbye[eBefore][i]) {
          ebye.fSixParticleCorrelationMaxEbye[eBefore][i] = ebye.fSixParticleCorrelationEbye[eBefore][i];
          LOGF(info, "max v%d^2v%d^2v%d^2 so far = %e", mcc.h4, mcc.h5, mcc.h6, ebye.fSixParticleCorrelationEbye[eBefore][i]);
        }
        if (ebye.fSixParticleCorrelationEbye[eBefore][i] < ebye.fSixParticleCorrelationMinEbye[eBefore][i]) {
          ebye.fSixParticleCorrelationMinEbye[eBefore][i] = ebye.fSixParticleCorrelationEbye[eBefore][i];
          LOGF(info, "min v%d^2v%d^2v%d^2 so far = %e", mcc.h4, mcc.h5, mcc.h6, ebye.fSixParticleCorrelationEbye[eBefore][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksBefore = %d, wSixRecursionBefore = %e", rlCollisionCent, nTracksBefore, wSixRecursionBefore);
      }

      // After cut:
      TComplex sixRecursionAfter = mccRecursion(6, harmonicsSixNum, eAfter) / mccRecursion(6, harmonicsSixDen, eAfter).Re();
      double wSixRecursionAfter = mccRecursion(6, harmonicsSixDen, eAfter).Re();
      ebye.fSixParticleCorrelationEbye[eAfter][i] = sixRecursionAfter.Re();

      if (nTracksAfter > k - 1 && wSixRecursionAfter > 0.) {
        // Fill v2^2v3^2v4^2 profile and histogram:
        mc.fSixParticleCorrelationProfiles[eAfter]->Fill(mc.fSixParticleCorrelationProfiles[eAfter]->GetXaxis()->GetBinCenter(i + 1), ebye.fSixParticleCorrelationEbye[eAfter][i], wSixRecursionAfter);
        mc.fSixParticleCorrelationHistograms[eAfter][eBefore][i]->Fill(ebye.fSixParticleCorrelationEbye[eAfter][i]);
        mc.fSixParticleCorrelationHistograms[eAfter][eAfter][i]->Fill(ebye.fSixParticleCorrelationEbye[eAfter][i], wSixRecursionAfter);

        // Find max and min:
        if (ebye.fSixParticleCorrelationEbye[eAfter][i] > ebye.fSixParticleCorrelationMaxEbye[eAfter][i]) {
          ebye.fSixParticleCorrelationMaxEbye[eAfter][i] = ebye.fSixParticleCorrelationEbye[eAfter][i];
          LOGF(info, "max v%d^2v%d^2v%d^2 so far = %e", mcc.h4, mcc.h5, mcc.h6, ebye.fSixParticleCorrelationEbye[eAfter][i]);
        }
        if (ebye.fSixParticleCorrelationEbye[eAfter][i] < ebye.fSixParticleCorrelationMinEbye[eAfter][i]) {
          ebye.fSixParticleCorrelationMinEbye[eAfter][i] = ebye.fSixParticleCorrelationEbye[eAfter][i];
          LOGF(info, "min v%d^2v%d^2v%d^2 so far = %e", mcc.h4, mcc.h5, mcc.h6, ebye.fSixParticleCorrelationEbye[eAfter][i]);
        }
      } else {
        LOGF(warning, "cent=%f, nTracksAfter = %d, wSixRecursionAfter = %e", rlCollisionCent, nTracksAfter, wSixRecursionAfter);
      }
    }
  }

  template <EnParticleHistograms histType, typename T1>
  void bookParticleHistograms(T1 const& lPcBins)
  {
    const auto& lPtBins = lPcBins[histType]; // define local array and initialize it from an array set in the configurables
    int nBinsPt = static_cast<int>(lPtBins[0]);
    float minPt = lPtBins[1];
    float maxPt = lPtBins[2];

    for (int ba = 0; ba < eBeforeAfter_N; ba++) {

      std::string nameRec = Form("fHist%s[eRec][%s cut]", ParticleHistNames[histType], BeforeAfterNames[ba]);
      std::string nameSim = Form("fHist%s[eSim][%s cut]", ParticleHistNames[histType], BeforeAfterNames[ba]);
      std::string nameRecfull = Form("%s distribution for reconstructed particles", ParticleHistNames[histType]);
      std::string nameSimfull = Form("%s distribution for simulated particles", ParticleHistNames[histType]);

      if (doprocessRec || doprocessRecSim) {
        pc.fParticleHistograms[histType][eRec][ba] = new TH1F(nameRec.c_str(), nameRecfull.c_str(), nBinsPt, minPt, maxPt);
        pc.fParticleHistograms[histType][eRec][ba]->GetXaxis()->SetTitle(ParticleHistNames[histType]);
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eRec][ba]);
      }

      if (doprocessSim || doprocessRecSim) {
        pc.fParticleHistograms[histType][eSim][ba] = new TH1F(nameSim.c_str(), nameSimfull.c_str(), nBinsPt, minPt, maxPt);
        pc.fParticleHistograms[histType][eSim][ba]->GetXaxis()->SetTitle(ParticleHistNames[histType]);
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eSim][ba]);
      }
    }
  }

  template <EnEventHistograms histType, typename T1>
  void bookEventHistograms(T1 const& lEvBins)
  {
    const auto& lCentBins = lEvBins[histType]; // define local array and initialize it from an array set in the configurables
    int nBinsCent = static_cast<int>(lCentBins[0]);
    float minCent = lCentBins[1];
    float maxCent = lCentBins[2];

    for (int ba = 0; ba < eBeforeAfter_N; ba++) {

      std::string nameRec = Form("fHist%s[eRec][%s cut]", EventHistNames[histType], BeforeAfterNames[ba]);
      std::string nameSim = Form("fHist%s[eSim][%s cut]", EventHistNames[histType], BeforeAfterNames[ba]);
      std::string nameRecfull;
      std::string nameSimfull;

      if constexpr (histType == eCent) {
        nameRecfull = Form("%s %s distribution for reconstructed events", tc.fCentEstm.c_str(), EventHistNames[histType]);
        nameSimfull = Form("%s %s distribution for simulated events", tc.fCentEstm.c_str(), EventHistNames[histType]);
      } else if constexpr (histType == eMult) {
        nameRecfull = Form("%s %s distribution for reconstructed events", tc.fMultEstm.c_str(), EventHistNames[histType]);
        nameSimfull = Form("%s %s distribution for simulated events", tc.fMultEstm.c_str(), EventHistNames[histType]);
      } else {
        nameRecfull = Form("%s distribution for reconstructed events", EventHistNames[histType]);
        nameSimfull = Form("%s distribution for simulated events", EventHistNames[histType]);
      }

      if (doprocessRec || doprocessRecSim) {
        ev.fEventHistograms[histType][eRec][ba] = new TH1F(nameRec.c_str(), nameRecfull.c_str(), nBinsCent, minCent, maxCent);
        ev.fEventHistograms[histType][eRec][ba]->GetXaxis()->SetTitle(EventHistNames[histType]);
        ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eRec][ba]);
      }

      if (doprocessSim || doprocessRecSim) {
        if constexpr (histType != eNumContrib && histType != eMult) {
          ev.fEventHistograms[histType][eSim][ba] = new TH1F(nameSim.c_str(), nameSimfull.c_str(), nBinsCent, minCent, maxCent);
          ev.fEventHistograms[histType][eSim][ba]->GetXaxis()->SetTitle(EventHistNames[histType]);
          ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eSim][ba]);
        } // No nContrib and multiplicity for processSim
      }
    }
  }

  template <EnEventHistograms histType, typename T1>
  void bookQAHistograms(T1 const& lQABins)
  {
    int nBinsCentX = 0;
    float minCentX = 0.;
    float maxCentX = 0.;
    int nBinsCentY = 0;
    float minCentY = 0.;
    float maxCentY = 0.;
    int nBinsCent = 0;
    float minCent = 0.;
    float maxCent = 0.;

    if (histType == eMult) {
      const auto& lCentBinsX = lQABins[1]; // MultBins
      nBinsCentX = static_cast<int>(lCentBinsX[0]);
      minCentX = lCentBinsX[1];
      maxCentX = lCentBinsX[2];
      const auto& lCentBinsY = lQABins[2]; // nContribBins
      nBinsCentY = static_cast<int>(lCentBinsY[0]);
      minCentY = lCentBinsY[1];
      maxCentY = lCentBinsY[2];
    } else {
      const auto& lCentBins = lQABins[histType];
      nBinsCent = static_cast<int>(lCentBins[0]);
      minCent = lCentBins[1];
      maxCent = lCentBins[2];
    }

    for (int ba = 0; ba < eBeforeAfter_N; ba++) {
      std::string name = Form("fHist%s[%s cut]", EventHistNames[histType], BeforeAfterNames[ba]);
      std::string namefull;

      if constexpr (histType == eCent) {
        namefull = Form("Quality assurance of %s %s", tc.fCentEstm.c_str(), EventHistNames[histType]);
      } else if constexpr (histType == eMult) {
        namefull = Form("Quality assurance of %s %s vs. NContributors", tc.fMultEstm.c_str(), EventHistNames[histType]);
      } else {
        namefull = Form("Quality assurance of %s", EventHistNames[histType]);
      }

      if constexpr (histType == eMult) {
        qa.fQAHistograms[histType][ba] = new TH2F(name.c_str(), namefull.c_str(), nBinsCentX, minCentX, maxCentX, nBinsCentY, minCentY, maxCentY);
        qa.fQAHistograms[histType][ba]->GetYaxis()->SetTitle("NContributors");
        qa.fQAHistograms[histType][ba]->GetXaxis()->SetTitle("Reference multiplicity");
      } else {
        qa.fQAHistograms[histType][ba] = new TH2F(name.c_str(), namefull.c_str(), nBinsCent, minCent, maxCent, nBinsCent, minCent, maxCent);
        qa.fQAHistograms[histType][ba]->GetYaxis()->SetTitle(Form("Simulated %s", EventHistNames[histType]));
        qa.fQAHistograms[histType][ba]->GetXaxis()->SetTitle(Form("Reconstructed %s", EventHistNames[histType]));
      }
      qa.fQAHistogramsList->Add(qa.fQAHistograms[histType][ba]);
    }
  }

  template <EnCorrHistograms histType, typename T1>
  void bookCorrHistograms(T1 const& lCrBins)
  {
    // (mult)(cent)(NContribMult) offset
    constexpr int BoundOffset = histType == eCorrCent ? eCentABBound : histType == eCorrMult ? eMultABBound
                                                                                             : eNumContribMultBound;

    // (mult)(cent)(NContribMult) book functions:
    if constexpr (histType == eCorrCent || histType == eCorrMult) {
      for (int i = 0; i < eEstmCorr_N; i++) {
        int indexBound = BoundOffset + i;
        cr.fCorrBounds[indexBound][0] = new TF1(Form("fCorrUpperBound%s%s", EstmCorrName[i], CorrHistNames[histType]), "[0] * x + [1]");
        cr.fCorrBounds[indexBound][1] = new TF1(Form("fCorrLowerBound%s%s", EstmCorrName[i], CorrHistNames[histType]), "[0] * x + [1]");
        cr.fCorrBounds[indexBound][0]->SetTitle(Form("%s %s upper bound", EstmCorrName[i], CorrHistNames[histType]));
        cr.fCorrBounds[indexBound][1]->SetTitle(Form("%s %s lower bound", EstmCorrName[i], CorrHistNames[histType]));
      }
    } else if constexpr (histType == eCorrNumContribMult) {
      int indexBound = BoundOffset;
      cr.fCorrBounds[indexBound][0] = new TF1(Form("fCorrUpperBound%s", CorrHistNames[histType]), "[0] * x + [1]");
      cr.fCorrBounds[indexBound][1] = new TF1(Form("fCorrLowerBound%s", CorrHistNames[histType]), "[0] * x + [1]");
      cr.fCorrBounds[indexBound][0]->SetTitle(Form("%s upper bound", CorrHistNames[histType]));
      cr.fCorrBounds[indexBound][1]->SetTitle(Form("%s lower bound", CorrHistNames[histType]));
    }

    int nBinsXCent = 0;
    float minXCent = 0.;
    float maxXCent = 0.;
    int nBinsYCent = 0;
    float minYCent = 0.;
    float maxYCent = 0.;
    int nBinsXMult = 0;
    float minXMult = 0.;
    float maxXMult = 0.;
    int nBinsYMult = 0;
    float minYMult = 0.;
    float maxYMult = 0.;
    // int nBinsNumContrib = 0;
    // float minNumContrib = 0.;
    float maxNumContrib = 0.;

    // (NContribMult) set functions and add to list:
    if constexpr (histType == eCorrNumContribMult) {

      int indexBins = static_cast<int>(histType) + static_cast<int>(eMultEstm_N);

      // (NContribMult) get bins:
      const auto& lMultXBins = lCrBins[indexBins];
      nBinsXMult = static_cast<int>(lMultXBins[0]);
      minXMult = lMultXBins[1];
      maxXMult = lMultXBins[2];
      const auto& lNumContribBins = lCrBins[indexBins - 1];
      // nBinsNumContrib = static_cast<int>(lNumContribBins[0]);
      // minNumContrib = lNumContribBins[1];
      maxNumContrib = lNumContribBins[2];

      // (NContribMult) set functions and add to list:
      int indexBound = BoundOffset;
      cr.fCorrBounds[indexBound][0]->SetRange(minXMult, maxXMult);
      cr.fCorrBounds[indexBound][1]->SetRange(minXMult, maxXMult);
      cr.fCorrBounds[indexBound][0]->SetMinimum(0);
      cr.fCorrBounds[indexBound][1]->SetMinimum(0);
      cr.fCorrBounds[indexBound][0]->SetMaximum(maxNumContrib);
      cr.fCorrBounds[indexBound][1]->SetMaximum(maxNumContrib);
      cr.fCorrBounds[indexBound][0]->SetParameters(tc.fNumContribMultCorrCut[0], tc.fNumContribMultCorrCut[1]);
      cr.fCorrBounds[indexBound][1]->SetParameters(tc.fNumContribMultCorrCut[2], tc.fNumContribMultCorrCut[3]);
      cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][0]);
      cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][1]);
    }

    // (mult)(cent) Book histograms, set functions and add to list:
    if constexpr (histType != eCorrNumContribMult) {
      for (int ba = 0; ba < eBeforeAfter_N; ba++) {
        std::string name;
        std::string namefull = Form("%s correlation %s cut", CorrHistNames[histType], BeforeAfterNames[ba]);

        int nEstm = 0;
        if constexpr (histType == eCorrCent) {
          nEstm = eCentEstm_N;
        } else if constexpr (histType == eCorrMult) {
          nEstm = eMultEstm_N;
        }

        int nComb = 0; // count combination AB, AC, BC

        for (int i = 0; i < nEstm; i++) {

          std::string titleX;
          if constexpr (histType == eCorrCent) {
            titleX = Form("%s %s", CentEstmNames[i], CorrHistNames[histType]);

          } else if constexpr (histType == eCorrMult) {
            titleX = Form("%s %s", MultEstmNames[i], CorrHistNames[histType]);
          }

          for (int j = i + 1; j < nEstm; j++) {

            std::string titleY;
            if constexpr (histType == eCorrCent) {
              name = Form("fHist%s[%s vs. %s][%s cut]", CorrHistNames[histType], CentEstmNames[i], CentEstmNames[j], BeforeAfterNames[ba]);
              titleY = Form("%s %s", CentEstmNames[j], CorrHistNames[histType]);

              // (cent) get bins:
              const auto& lCentXBins = lCrBins[0];
              nBinsXCent = static_cast<int>(lCentXBins[0]);
              minXCent = lCentXBins[1];
              maxXCent = lCentXBins[2];
              const auto& lCentYBins = lCrBins[0];
              nBinsYCent = static_cast<int>(lCentYBins[0]);
              minYCent = lCentYBins[1];
              maxYCent = lCentYBins[2];

              // (cent) set functions and add to list:
              if (static_cast<bool>(ba)) {
                int indexBound = BoundOffset + nComb;
                cr.fCorrBounds[indexBound][0]->SetRange(minXCent, maxXCent);
                cr.fCorrBounds[indexBound][1]->SetRange(minXCent, maxXCent);
                cr.fCorrBounds[indexBound][0]->SetParameters(tc.fCentAllCorrCut[i][j][0], tc.fCentAllCorrCut[i][j][1]);
                cr.fCorrBounds[indexBound][1]->SetParameters(tc.fCentAllCorrCut[i][j][2], tc.fCentAllCorrCut[i][j][3]);
                cr.fCorrBounds[indexBound][0]->SetMinimum(0);
                cr.fCorrBounds[indexBound][1]->SetMinimum(0);
                cr.fCorrBounds[indexBound][0]->SetMaximum(maxYCent);
                cr.fCorrBounds[indexBound][1]->SetMaximum(maxYCent);
                cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][0]);
                cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][1]);
              }

              // (cent) book corr 2D histogram:
              cr.fCorrHistograms[histType][i][j][ba] = new TH2F(name.c_str(), namefull.c_str(), nBinsXCent, minXCent, maxXCent, nBinsYCent, minYCent, maxYCent);

            } else if constexpr (histType == eCorrMult) {
              name = Form("fHist%s[%s vs. %s][%s cut]", CorrHistNames[histType], MultEstmNames[i], MultEstmNames[j], BeforeAfterNames[ba]);
              titleY = Form("%s %s", MultEstmNames[j], CorrHistNames[histType]);

              // (mult) get bins:
              const auto& lMultXBins = lCrBins[i + 1];
              nBinsXMult = static_cast<int>(lMultXBins[0]);
              minXMult = lMultXBins[1];
              maxXMult = lMultXBins[2];
              const auto& lMultYBins = lCrBins[j + 1];
              nBinsYMult = static_cast<int>(lMultYBins[0]);
              minYMult = lMultYBins[1];
              maxYMult = lMultYBins[2];

              // (mult) set functions and add to list:
              if (static_cast<bool>(ba)) {
                int indexBound = BoundOffset + nComb;
                cr.fCorrBounds[indexBound][0]->SetRange(minXMult, maxXMult);
                cr.fCorrBounds[indexBound][1]->SetRange(minXMult, maxXMult);
                cr.fCorrBounds[indexBound][0]->SetParameters(tc.fMultAllCorrCut[i][j][0], tc.fMultAllCorrCut[i][j][1]);
                cr.fCorrBounds[indexBound][1]->SetParameters(tc.fMultAllCorrCut[i][j][2], tc.fMultAllCorrCut[i][j][3]);
                cr.fCorrBounds[indexBound][0]->SetMinimum(0);
                cr.fCorrBounds[indexBound][1]->SetMinimum(0);
                cr.fCorrBounds[indexBound][0]->SetMaximum(maxYMult);
                cr.fCorrBounds[indexBound][1]->SetMaximum(maxYMult);
                cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][0]);
                cr.fCorrHistogramsList->Add(cr.fCorrBounds[indexBound][1]);
              }

              // (mult) book corr 2D histogram:
              cr.fCorrHistograms[histType][i][j][ba] = new TH2F(name.c_str(), namefull.c_str(), nBinsXMult, minXMult, maxXMult, nBinsYMult, minYMult, maxYMult);
            }

            // (cent)(mult) set corr histogram
            cr.fCorrHistograms[histType][i][j][ba]->GetYaxis()->SetTitle(titleY.c_str());
            cr.fCorrHistograms[histType][i][j][ba]->GetXaxis()->SetTitle(titleX.c_str());
            cr.fCorrHistogramsList->Add(cr.fCorrHistograms[histType][i][j][ba]);

            nComb += 1;
          }
        }
      }
    }
  }

  // *) Initialize and book all objects:
  void init(InitContext&)
  {

    // ... code to book and initialize all analysis objects ...

    // *) Set automatically what to process, from an implicit variable
    // "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[eProcessRec] = doprocessRec;
    tc.fProcess[eProcessRecSim] = doprocessRecSim;
    tc.fProcess[eProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;
    tc.fCentEstm = cfCentEstm;
    tc.fMultEstm = cfMultEstm;

    tc.fGeneralCutSwitch = cfGeneralCutSwitch;
    tc.fVertexZCutSwitch = cfVertexZCutSwitch;
    tc.fSel8CutSwitch = cfSel8CutSwitch;
    tc.fCentCutSwitch = cfCentCutSwitch;
    tc.fNumContribCutSwitch = cfNumContribCutSwitch;
    tc.fCentCorrCutSwitch = cfCentCorrCutSwitch;
    tc.fMultCorrCutSwitch = cfMultCorrCutSwitch;
    tc.fNumContribMultCorrCutSwitch = cfNumContribMultCorrCutSwitch;

    tc.fPtCutSwitch = cfPtCutSwitch;
    tc.fEtaCutSwitch = cfEtaCutSwitch;
    tc.fSignCutSwitch = cfSignCutSwitch;
    tc.fTpcNClsFoundCutSwitch = cfTpcNClsFoundCutSwitch;
    tc.fDCAXYCutSwitch = cfDCAXYCutSwitch;
    tc.fDCAZCutSwitch = cfDCAZCutSwitch;
    tc.fEtaGapSwitch = cfEtaGapSwitch;

    tc.fVertexZCut = cfVertexZCut;
    tc.fCentCut = cfCentCut;
    tc.fNumContribCut = cfNumContribCut;
    tc.fCentAllCorrCut[0][1] = cfCentABCorrCut;
    tc.fCentAllCorrCut[0][2] = cfCentACCorrCut;
    tc.fCentAllCorrCut[1][2] = cfCentBCCorrCut;
    tc.fMultAllCorrCut[0][1] = cfMultABCorrCut;
    tc.fMultAllCorrCut[0][2] = cfMultACCorrCut;
    tc.fMultAllCorrCut[1][2] = cfMultBCCorrCut;
    tc.fNumContribMultCorrCut = cfNumContribMultCorrCut;

    tc.fPtCut = cfPtCut;
    tc.fEtaCut = cfEtaCut;
    tc.fSignCut = cfSignCut;
    tc.fTpcNClsFoundCut = cfTpcNClsFoundCut;
    tc.fDCAXYCut = cfDCAXYCut;
    tc.fDCAZCut = cfDCAZCut;
    tc.fEtaGap = cfEtaGap;

    tc.fPtBins = cfPtBins;
    tc.fPhiBins = cfPhiBins;

    tc.fCentBins = cfCentBins;
    tc.fFT0CMultBins = cfFT0CMultBins;
    tc.fFT0MMultBins = cfFT0MMultBins;
    tc.fFV0AMultBins = cfFV0AMultBins;

    tc.fVerXBins = cfVerXBins;
    tc.fVerYBins = cfVerYBins;
    tc.fVerZBins = cfVerZBins;
    tc.fNumContribBins = cfNumContribBins;
    tc.fTwoParticleCorrBins = cfTwoParticleCorrBins;
    tc.fFourParticleCorrBins = cfFourParticleCorrBins;
    tc.fSixParticleCorrBins = cfSixParticleCorrBins;

    tc.fPrintSwitch = cfPrintSwitch;
    tc.fFileWithWeights = cfFileWithWeights;

    qa.fQASwitch = cfQASwitch;
    wt.fWeightSwitch = cfWeightSwitch;

    if (tc.fMultEstm == "FT0C") {
      tc.fMultBins = tc.fFT0CMultBins;
    } else if (tc.fMultEstm == "FT0M") {
      tc.fMultBins = tc.fFT0MMultBins;
    } else if (tc.fMultEstm == "FV0A") {
      tc.fMultBins = tc.fFV0AMultBins;
    } else {
      LOGF(fatal, "Invalid multiplicity estimator: %s", tc.fMultEstm.c_str());
    }

    // *) Book base list:
    auto* temp = new TList();
    temp->SetOwner(kTRUE);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    pc.fParticleHistogramsList = new TList();
    pc.fParticleHistogramsList->SetName("ParticleHistograms");
    pc.fParticleHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(pc.fParticleHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    ev.fEventHistogramsList = new TList();
    ev.fEventHistogramsList->SetName("EventHistograms");
    ev.fEventHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(ev.fEventHistogramsList);

    qa.fQAHistogramsList = new TList();
    qa.fQAHistogramsList->SetName("QualityAssuranceHistograms");
    qa.fQAHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(qa.fQAHistogramsList);

    wt.fWeightHistogramsList = new TList();
    wt.fWeightHistogramsList->SetName("WeightHistograms");
    wt.fWeightHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(wt.fWeightHistogramsList);

    cr.fCorrHistogramsList = new TList();
    cr.fCorrHistogramsList->SetName("CorrelationHistograms");
    cr.fCorrHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(cr.fCorrHistogramsList);

    mc.fMultiparticleCorrelationProfilesList = new TList();
    mc.fMultiparticleCorrelationProfilesList->SetName("MultiparticleCorrelationProfiles");
    mc.fMultiparticleCorrelationProfilesList->SetOwner(kTRUE);
    fBaseList->Add(mc.fMultiparticleCorrelationProfilesList);

    std::vector<std::vector<float>> lPcBins = {tc.fPtBins, tc.fPhiBins};
    std::vector<std::vector<float>> lEvBins = {tc.fCentBins, tc.fMultBins, tc.fVerXBins, tc.fVerYBins, tc.fVerZBins, tc.fNumContribBins};
    std::vector<std::vector<float>> lQABins = {tc.fCentBins, tc.fMultBins, tc.fNumContribBins};
    std::vector<std::vector<float>> lCrBins = {tc.fCentBins, tc.fFT0CMultBins, tc.fFT0MMultBins, tc.fFV0AMultBins, tc.fNumContribBins, tc.fMultBins};

    bookParticleHistograms<ePt>(lPcBins);
    bookParticleHistograms<ePhi>(lPcBins);
    bookEventHistograms<eCent>(lEvBins);
    bookEventHistograms<eMult>(lEvBins);
    bookEventHistograms<eVertexX>(lEvBins);
    bookEventHistograms<eVertexY>(lEvBins);
    bookEventHistograms<eVertexZ>(lEvBins);
    bookEventHistograms<eNumContrib>(lEvBins);
    bookQAHistograms<eCent>(lQABins);
    bookQAHistograms<eMult>(lQABins);
    bookCorrHistograms<eCorrCent>(lCrBins);
    bookCorrHistograms<eCorrMult>(lCrBins);
    bookCorrHistograms<eCorrNumContribMult>(lCrBins);

    wt.fDummyPhiWeightHistogram = new TH1F("fDummyPhiWeightHistogram", "Dummy phi weight histogram", static_cast<int>(tc.fPhiBins[0]), tc.fPhiBins[1], tc.fPhiBins[2]);
    for (int i = 1; i <= wt.fDummyPhiWeightHistogram->GetNbinsX(); i++) {
      wt.fDummyPhiWeightHistogram->SetBinContent(i, 1.);
    }
    wt.fDummyPtWeightHistogram = new TH1F("fDummyPtWeightHistogram", "Dummy pt weight histogram", static_cast<int>(tc.fPtBins[0]), tc.fPtBins[1], tc.fPtBins[2]);
    for (int i = 1; i <= wt.fDummyPtWeightHistogram->GetNbinsX(); i++) {
      wt.fDummyPtWeightHistogram->SetBinContent(i, 1.);
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...
    // *) Steer all analysis steps:
    runLoop<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCumulants, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    runLoop<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCumulants, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // runLoop<eSim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCumulants, processSim, "process only simulated data", false);

}; // struct MultiparticleCumulants {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiparticleCumulants>(cfgc)};
} // WorkflowSpec...
