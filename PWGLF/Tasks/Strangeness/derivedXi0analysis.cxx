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
//
/// \file derivedXi0Analysis.cxx
/// \brief Xi0 (--> Lambda pi0 --> Lambda gamma gamma) analysis task using strangeness derived data (produced with the Xi0 builder)
///
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI
//
// Xi0 analysis task
// =================
//
// This code loops over the Xi0Cores table produced by the sigma0builder
// (PWGLF/TableProducer/Strangeness/sigma0builder.cxx, fillXi0Tables = true)
// and produces some standard analysis output. It is meant to be run over
// strangeness derived data.
//
// The Xi0 candidate is built out of three V0s: two photon conversions
// (forming the pi0) and one Lambda. The daughter V0s are reached by
// de-referencing the V0Cores table through the Xi0Indices table.
//
// Three process functions are provided:
//   - processRealData   : real data
//   - processMonteCarlo : reconstructed information in MC
//   - processGenerated  : pure generated information in MC (from CascMCCores)
//
// N.B.: when running over reconstructed MC information, the rapidity and the
//       pT used both in the selections and in the histograms are the generated
//       ones, so that the numerator and the denominator of the
//       acceptance x efficiency share exactly the same definition.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//

#include "PWGLF/DataModel/LFSigmaTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// simple checkers, but ensure 64 bit integers
#define BITSET(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define BITCHECK(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

// Xi0 candidates: cores + collision reference + indices to the daughter V0s.
// N.B.: aod::Xi0Indices and aod::Xi0MCIndices must never be joined together,
//       since they share the very same column names (photon1Index, ...).
using Xi0Candidates = soa::Join<aod::Xi0Cores, aod::Xi0CollRefs, aod::Xi0Indices>;
using Xi0McCandidates = soa::Join<aod::Xi0Cores, aod::Xi0CollRefs, aod::Xi0Indices, aod::Xi0MCCores>;

using V0Candidates = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras>;
using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

enum CentEstimator {
  kCentFT0C = 0,
  kCentFT0M,
  kCentFT0CVariant1,
  kCentMFT,
  kCentNGlobal,
  kCentFV0A
};

struct DerivedXi0Analysis {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Xi0 proper decay length, in cm (PDG: tau = 2.90e-10 s)
  static constexpr float CtauXi0 = 8.71f;
  // number of QA histogram sets: before and after the candidate selections
  static constexpr int NSelectionStages = 2;

  //__________________________________________________
  // Selection bits: one per selection criterion, so that the full selection is
  // a single mask operation and the bookkeeping fits in a single histogram
  enum Selectionbits : int { // o2-linter: disable=name/enum (bit names follow derivedlambdakzeroanalysis)
    // photon selections (both photons of the pi0 have to pass)
    selPhotonV0Type = 0,
    selPhotonMass,
    selPhotonRapidity,
    selPhotonDauEta,
    selPhotonDCADauToPV,
    selPhotonDCADau,
    selPhotonRadius,
    selPhotonCosPA,
    selPhotonArmenteros,
    selPhotonPt,
    selPhotonTPCCrossedRows,
    selPhotonTPCPID,
    // Lambda selections, common to Lambda and anti-Lambda
    selLambdaV0Type,
    selLambdaRapidity,
    selLambdaDauEta,
    selLambdaDCAPosToPV,
    selLambdaDCANegToPV,
    selLambdaDCADau,
    selLambdaRadius,
    selLambdaCosPA,
    selLambdaLifetime,
    selLambdaTPCCrossedRows,
    selLambdaITSClusters,
    // Lambda selections, specific to one of the two mass hypotheses
    selLambdaMass,
    selAntiLambdaMass,
    selLambdaArmenteros,     // Armenteros alpha > 0 --> the positive daughter is the baryon
    selAntiLambdaArmenteros, // Armenteros alpha < 0 --> the negative daughter is the baryon
    selTPCPIDPositiveProton,
    selTPCPIDNegativePion,
    selTPCPIDNegativeProton,
    selTPCPIDPositivePion,
    // pi0 and Xi0 (cascade) selections
    selPi0Mass,
    selPi0Radius,
    selPi0CosPA,
    selPi0DCADau,
    selXi0Radius,
    selXi0CosPA,
    selXi0DCADau,
    selXi0DCAxyToPV,
    selXi0DCAzToPV,
    selXi0Lifetime,
    selXi0Rapidity,
    selXi0Pt,
    // MC tagging
    selConsiderXi0,
    selConsiderAntiXi0,
    selPhysPrimXi0,
    selPhysPrimAntiXi0,
  };

  uint64_t maskTopological = 0;
  uint64_t maskXi0Specific = 0;
  uint64_t maskAntiXi0Specific = 0;
  uint64_t maskSelectionXi0 = 0;
  uint64_t maskSelectionAntiXi0 = 0;

  //__________________________________________________
  // Event level
  Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0M, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal, 5:CentFV0A)"};
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};
    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
    Configurable<bool> useEvtSelInDenomEff{"useEvtSelInDenomEff", false, "Consider event selections in the recoed <-> gen collision association for the denominator (or numerator) of the acc. x eff. (or signal loss)?"};
    Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Apply Z-vtx cut on the PV of the generated collision?"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};
  } eventSelections;

  //__________________________________________________
  // Photon (from the pi0) selections, applied on the de-referenced V0Cores
  struct : ConfigurableGroup {
    std::string prefix = "photonSelections"; // JSON group name
    Configurable<int> v0TypeSelection{"v0TypeSelection", 7, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> maxMass{"maxMass", 0.1, "Max photon mass (GeV/c^2)"};
    Configurable<float> minRapidity{"minRapidity", -0.8, "Min photon rapidity"};
    Configurable<float> maxRapidity{"maxRapidity", 0.8, "Max photon rapidity"};
    Configurable<float> maxDauEta{"maxDauEta", 0.8, "Max |eta| of the daughter tracks"};
    Configurable<float> minDCADauToPV{"minDCADauToPV", 0.0, "Min DCA of the daughter tracks to PV (cm)"};
    Configurable<float> maxDCADau{"maxDCADau", 3.5, "Max DCA between the V0 daughters (cm)"};
    Configurable<float> minRadius{"minRadius", 0.0, "Min photon conversion radius (cm)"};
    Configurable<float> maxRadius{"maxRadius", 240.0, "Max photon conversion radius (cm)"};
    Configurable<float> minCosPA{"minCosPA", 0.8, "Min photon cosine of pointing angle"};
    Configurable<float> maxQt{"maxQt", 0.05, "Max qt (Armenteros-Podolanski) (GeV/c)"};
    Configurable<float> maxAlpha{"maxAlpha", 0.95, "Max |alpha| (Armenteros-Podolanski)"};
    Configurable<float> minPt{"minPt", 0.0, "Min photon pT (GeV/c)"};
    Configurable<float> maxPt{"maxPt", 50.0, "Max photon pT (GeV/c)"};
    Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 30, "Min number of TPC crossed rows of the daughter tracks"};
    Configurable<float> minTPCNSigmaEl{"minTPCNSigmaEl", -7, "Min TPC NSigma (electron) of the daughter tracks"};
    Configurable<float> maxTPCNSigmaEl{"maxTPCNSigmaEl", +7, "Max TPC NSigma (electron) of the daughter tracks"};
  } photonSelections;

  //__________________________________________________
  // Lambda selections, applied on the de-referenced V0Cores
  struct : ConfigurableGroup {
    std::string prefix = "lambdaSelections"; // JSON group name
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> massWindow{"massWindow", 0.015, "Lambda mass window around the PDG value (GeV/c^2)"};
    Configurable<float> minRapidity{"minRapidity", -0.8, "Min Lambda rapidity"};
    Configurable<float> maxRapidity{"maxRapidity", 0.8, "Max Lambda rapidity"};
    Configurable<float> maxDauEta{"maxDauEta", 0.8, "Max |eta| of the daughter tracks"};
    Configurable<float> minDCAPosToPV{"minDCAPosToPV", 0.05, "Min DCA of the positive daughter to PV (cm)"};
    Configurable<float> minDCANegToPV{"minDCANegToPV", 0.05, "Min DCA of the negative daughter to PV (cm)"};
    Configurable<float> maxDCADau{"maxDCADau", 1.0, "Max DCA between the V0 daughters (cm)"};
    Configurable<float> minRadius{"minRadius", 0.5, "Min Lambda decay radius (cm)"};
    Configurable<float> maxRadius{"maxRadius", 200.0, "Max Lambda decay radius (cm)"};
    Configurable<float> minCosPA{"minCosPA", 0.95, "Min Lambda cosine of pointing angle (w.r.t. the PV)"};
    Configurable<float> maxLifetime{"maxLifetime", 30.0, "Max Lambda proper lifetime, m*L/p (cm)"};
    Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Min number of TPC crossed rows of the daughter tracks"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "Min number of ITS clusters of the daughter tracks (leave negative if no selection desired)"};
    Configurable<float> maxTPCNSigmaPr{"maxTPCNSigmaPr", 5, "Max |TPC NSigma| (proton) of the baryon daughter"};
    Configurable<float> maxTPCNSigmaPi{"maxTPCNSigmaPi", 5, "Max |TPC NSigma| (pion) of the meson daughter"};
    // Lambda / anti-Lambda are told apart with the Armenteros-Podolanski alpha:
    // the baryon takes most of the momentum, so alpha > 0 for Lambda (positive
    // daughter = proton) and alpha < 0 for anti-Lambda (negative daughter = anti-proton)
    Configurable<float> minAbsAlpha{"minAbsAlpha", 0.25, "Min |alpha| (Armenteros-Podolanski) to tag the Lambda charge"};
    Configurable<float> maxAbsAlpha{"maxAbsAlpha", 1.0, "Max |alpha| (Armenteros-Podolanski)"};
  } lambdaSelections;

  //__________________________________________________
  // Xi0 selections, applied on the cascade itself
  struct : ConfigurableGroup {
    std::string prefix = "xi0Selections"; // JSON group name
    Configurable<float> minRapidity{"minRapidity", -0.5, "Min Xi0 rapidity"};
    Configurable<float> maxRapidity{"maxRapidity", 0.5, "Max Xi0 rapidity"};
    Configurable<float> minPt{"minPt", 0.0, "Min Xi0 pT (GeV/c)"};
    Configurable<float> maxPt{"maxPt", 50.0, "Max Xi0 pT (GeV/c)"};
    Configurable<float> pi0MassWindow{"pi0MassWindow", 0.035, "pi0 mass window around the PDG value (GeV/c^2)"};
    Configurable<float> minPi0Radius{"minPi0Radius", -1, "Min pi0 decay radius (cm), leave negative if no selection desired"};
    Configurable<float> maxPi0Radius{"maxPi0Radius", -1, "Max pi0 decay radius (cm), leave negative if no selection desired"};
    Configurable<float> minPi0CosPA{"minPi0CosPA", -2, "Min pi0 cosine of pointing angle, leave below -1 if no selection desired"};
    Configurable<float> maxDCAPi0Daughters{"maxDCAPi0Daughters", -1, "Max DCA between the two photons (cm), leave negative if no selection desired"};
    Configurable<float> minCascRadius{"minCascRadius", 0.5, "Min Xi0 decay radius (cm)"};
    Configurable<float> maxCascRadius{"maxCascRadius", 200.0, "Max Xi0 decay radius (cm)"};
    Configurable<float> minCascCosPA{"minCascCosPA", 0.98, "Min Xi0 cosine of pointing angle"};
    Configurable<float> maxDCACascDaughters{"maxDCACascDaughters", 1.0, "Max DCA between the Xi0 daughters (cm)"};
    Configurable<float> maxDCAxyCascToPV{"maxDCAxyCascToPV", 1.0, "Max |DCAxy| of the Xi0 to the PV (cm), leave negative if no selection desired"};
    Configurable<float> maxDCAzCascToPV{"maxDCAzCascToPV", -1, "Max |DCAz| of the Xi0 to the PV (cm), leave negative if no selection desired"};
    Configurable<float> maxLifetime{"maxLifetime", 3.0, "Max Xi0 proper lifetime, in units of c*tau"};
  } xi0Selections;

  //__________________________________________________
  // MC-specific selections
  struct : ConfigurableGroup {
    std::string prefix = "mcSelections"; // JSON group name
    Configurable<bool> doMCAssociation{"doMCAssociation", false, "Keep only candidates whose full decay chain is matched to a true Xi0/anti-Xi0"};
    Configurable<bool> requirePhysicalPrimary{"requirePhysicalPrimary", true, "Keep only physical primary Xi0, at reconstructed and generated level"};
  } mcSelections;

  //__________________________________________________
  // QA switches
  Configurable<bool> doEventQA{"doEventQA", true, "do event QA histograms"};
  Configurable<bool> doCandidateQA{"doCandidateQA", true, "do candidate-level QA histograms"};
  Configurable<bool> doMCQA{"doMCQA", true, "do MC-specific QA histograms (pT resolution, background composition)"};

  //__________________________________________________
  // Axes
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Centrality (%)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 15.0f, 20.0f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisXi0Mass{"axisXi0Mass", {200, 1.2f, 1.5f}, "#it{M}_{#Lambda#pi^{0}} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPi0Mass{"axisPi0Mass", {200, 0.0f, 0.4f}, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.07f, 1.17f}, "#it{M}_{p#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "Decay radius (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.9f, 1.0f}, "cos(#theta_{PA})"};
  ConfigurableAxis axisDCADau{"axisDCADau", {200, 0.0f, 5.0f}, "DCA between daughters (cm)"};
  ConfigurableAxis axisDCAToPV{"axisDCAToPV", {200, -5.0f, 5.0f}, "DCA to PV (cm)"};
  ConfigurableAxis axisLifetime{"axisLifetime", {200, 0.0f, 50.0f}, "Proper lifetime (cm)"};
  ConfigurableAxis axisRapidity{"axisRapidity", {200, -1.0f, 1.0f}, "#it{y}"};
  ConfigurableAxis axisArmAlpha{"axisArmAlpha", {200, -1.0f, 1.0f}, "#alpha (Armenteros-Podolanski)"};
  ConfigurableAxis axisArmQt{"axisArmQt", {200, 0.0f, 0.5f}, "#it{q}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "#it{N}_{ch}"};
  ConfigurableAxis axisPtResolution{"axisPtResolution", {200, -1.0f, 1.0f}, "(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20.0f, 20.0f}, "PV #it{z} (cm)"};
  ConfigurableAxis axisPDGCode{"axisPDGCode", {10001, -5000.5f, 5000.5f}, "PDG code"};
  ConfigurableAxis axisNCandidates{"axisNCandidates", {10, -0.5f, 9.5f}, "Number of candidates"};

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    //_______________________________________________
    // Assemble the selection masks once and for all
    maskTopological = 0;
    BITSET(maskTopological, selPhotonV0Type);
    BITSET(maskTopological, selPhotonMass);
    BITSET(maskTopological, selPhotonRapidity);
    BITSET(maskTopological, selPhotonDauEta);
    BITSET(maskTopological, selPhotonDCADauToPV);
    BITSET(maskTopological, selPhotonDCADau);
    BITSET(maskTopological, selPhotonRadius);
    BITSET(maskTopological, selPhotonCosPA);
    BITSET(maskTopological, selPhotonArmenteros);
    BITSET(maskTopological, selPhotonPt);
    BITSET(maskTopological, selPhotonTPCCrossedRows);
    BITSET(maskTopological, selPhotonTPCPID);
    BITSET(maskTopological, selLambdaV0Type);
    BITSET(maskTopological, selLambdaRapidity);
    BITSET(maskTopological, selLambdaDauEta);
    BITSET(maskTopological, selLambdaDCAPosToPV);
    BITSET(maskTopological, selLambdaDCANegToPV);
    BITSET(maskTopological, selLambdaDCADau);
    BITSET(maskTopological, selLambdaRadius);
    BITSET(maskTopological, selLambdaCosPA);
    BITSET(maskTopological, selLambdaLifetime);
    BITSET(maskTopological, selLambdaTPCCrossedRows);
    BITSET(maskTopological, selLambdaITSClusters);
    BITSET(maskTopological, selPi0Mass);
    BITSET(maskTopological, selPi0Radius);
    BITSET(maskTopological, selPi0CosPA);
    BITSET(maskTopological, selPi0DCADau);
    BITSET(maskTopological, selXi0Radius);
    BITSET(maskTopological, selXi0CosPA);
    BITSET(maskTopological, selXi0DCADau);
    BITSET(maskTopological, selXi0DCAxyToPV);
    BITSET(maskTopological, selXi0DCAzToPV);
    BITSET(maskTopological, selXi0Lifetime);
    BITSET(maskTopological, selXi0Rapidity);
    BITSET(maskTopological, selXi0Pt);

    maskXi0Specific = 0;
    BITSET(maskXi0Specific, selLambdaMass);
    BITSET(maskXi0Specific, selLambdaArmenteros);
    BITSET(maskXi0Specific, selTPCPIDPositiveProton);
    BITSET(maskXi0Specific, selTPCPIDNegativePion);
    BITSET(maskXi0Specific, selConsiderXi0);

    maskAntiXi0Specific = 0;
    BITSET(maskAntiXi0Specific, selAntiLambdaMass);
    BITSET(maskAntiXi0Specific, selAntiLambdaArmenteros);
    BITSET(maskAntiXi0Specific, selTPCPIDNegativeProton);
    BITSET(maskAntiXi0Specific, selTPCPIDPositivePion);
    BITSET(maskAntiXi0Specific, selConsiderAntiXi0);

    if (mcSelections.requirePhysicalPrimary) {
      BITSET(maskXi0Specific, selPhysPrimXi0);
      BITSET(maskAntiXi0Specific, selPhysPrimAntiXi0);
    }

    maskSelectionXi0 = maskTopological | maskXi0Specific;
    maskSelectionAntiXi0 = maskTopological | maskAntiXi0Specific;

    //_______________________________________________
    // Event selection bookkeeping
    auto hEventSelection = histos.add<TH1>("hEventSelection", "hEventSelection", kTH1D, {{23, -0.5f, +22.5f}});
    hEventSelection->GetXaxis()->SetBinLabel(1, "All collisions");
    hEventSelection->GetXaxis()->SetBinLabel(2, "sel8 cut");
    hEventSelection->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    hEventSelection->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    hEventSelection->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    hEventSelection->GetXaxis()->SetBinLabel(6, "posZ cut");
    hEventSelection->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
    hEventSelection->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    hEventSelection->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
    hEventSelection->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
    hEventSelection->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
    hEventSelection->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
    hEventSelection->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
    hEventSelection->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
    hEventSelection->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
    hEventSelection->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
    hEventSelection->GetXaxis()->SetBinLabel(17, "INEL>0");
    hEventSelection->GetXaxis()->SetBinLabel(18, "INEL>1");
    hEventSelection->GetXaxis()->SetBinLabel(19, "Below min occup.");
    hEventSelection->GetXaxis()->SetBinLabel(20, "Above max occup.");
    hEventSelection->GetXaxis()->SetBinLabel(21, "Below min IR");
    hEventSelection->GetXaxis()->SetBinLabel(22, "Above max IR");
    hEventSelection->GetXaxis()->SetBinLabel(23, "Selected collisions");

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});
    histos.add("hEventPVz", "hEventPVz", kTH1D, {axisVtxZ});
    if (doEventQA) {
      histos.add("EventQA/hCentralityVsNch", "hCentralityVsNch", kTH2D, {axisCentrality, axisNch});
      histos.add("EventQA/hCentralityVsPVz", "hCentralityVsPVz", kTH2D, {axisCentrality, axisVtxZ});
    }

    //_______________________________________________
    // Single candidate-selection bookkeeping histogram
    auto hSelections = histos.add<TH1>("hSelections", "hSelections", kTH1D, {{static_cast<int>(selPhysPrimAntiXi0) + 3, -0.5f, static_cast<float>(selPhysPrimAntiXi0) + 2.5f}});
    hSelections->GetXaxis()->SetBinLabel(1, "All");
    hSelections->GetXaxis()->SetBinLabel(selPhotonV0Type + 2, "#gamma V0 type");
    hSelections->GetXaxis()->SetBinLabel(selPhotonMass + 2, "#gamma mass");
    hSelections->GetXaxis()->SetBinLabel(selPhotonRapidity + 2, "#gamma rapidity");
    hSelections->GetXaxis()->SetBinLabel(selPhotonDauEta + 2, "#gamma dau. #eta");
    hSelections->GetXaxis()->SetBinLabel(selPhotonDCADauToPV + 2, "#gamma DCA dau. to PV");
    hSelections->GetXaxis()->SetBinLabel(selPhotonDCADau + 2, "#gamma DCA dau.");
    hSelections->GetXaxis()->SetBinLabel(selPhotonRadius + 2, "#gamma radius");
    hSelections->GetXaxis()->SetBinLabel(selPhotonCosPA + 2, "#gamma cosPA");
    hSelections->GetXaxis()->SetBinLabel(selPhotonArmenteros + 2, "#gamma Arm. pod.");
    hSelections->GetXaxis()->SetBinLabel(selPhotonPt + 2, "#gamma #it{p}_{T}");
    hSelections->GetXaxis()->SetBinLabel(selPhotonTPCCrossedRows + 2, "#gamma TPC rows");
    hSelections->GetXaxis()->SetBinLabel(selPhotonTPCPID + 2, "#gamma TPC PID");
    hSelections->GetXaxis()->SetBinLabel(selLambdaV0Type + 2, "#Lambda V0 type");
    hSelections->GetXaxis()->SetBinLabel(selLambdaRapidity + 2, "#Lambda rapidity");
    hSelections->GetXaxis()->SetBinLabel(selLambdaDauEta + 2, "#Lambda dau. #eta");
    hSelections->GetXaxis()->SetBinLabel(selLambdaDCAPosToPV + 2, "#Lambda DCA pos. to PV");
    hSelections->GetXaxis()->SetBinLabel(selLambdaDCANegToPV + 2, "#Lambda DCA neg. to PV");
    hSelections->GetXaxis()->SetBinLabel(selLambdaDCADau + 2, "#Lambda DCA dau.");
    hSelections->GetXaxis()->SetBinLabel(selLambdaRadius + 2, "#Lambda radius");
    hSelections->GetXaxis()->SetBinLabel(selLambdaCosPA + 2, "#Lambda cosPA");
    hSelections->GetXaxis()->SetBinLabel(selLambdaLifetime + 2, "#Lambda lifetime");
    hSelections->GetXaxis()->SetBinLabel(selLambdaTPCCrossedRows + 2, "#Lambda TPC rows");
    hSelections->GetXaxis()->SetBinLabel(selLambdaITSClusters + 2, "#Lambda ITS clusters");
    hSelections->GetXaxis()->SetBinLabel(selLambdaMass + 2, "#Lambda mass");
    hSelections->GetXaxis()->SetBinLabel(selAntiLambdaMass + 2, "#bar{#Lambda} mass");
    hSelections->GetXaxis()->SetBinLabel(selLambdaArmenteros + 2, "#Lambda Arm. pod.");
    hSelections->GetXaxis()->SetBinLabel(selAntiLambdaArmenteros + 2, "#bar{#Lambda} Arm. pod.");
    hSelections->GetXaxis()->SetBinLabel(selTPCPIDPositiveProton + 2, "TPC PID p");
    hSelections->GetXaxis()->SetBinLabel(selTPCPIDNegativePion + 2, "TPC PID #pi^{-}");
    hSelections->GetXaxis()->SetBinLabel(selTPCPIDNegativeProton + 2, "TPC PID #bar{p}");
    hSelections->GetXaxis()->SetBinLabel(selTPCPIDPositivePion + 2, "TPC PID #pi^{+}");
    hSelections->GetXaxis()->SetBinLabel(selPi0Mass + 2, "#pi^{0} mass");
    hSelections->GetXaxis()->SetBinLabel(selPi0Radius + 2, "#pi^{0} radius");
    hSelections->GetXaxis()->SetBinLabel(selPi0CosPA + 2, "#pi^{0} cosPA");
    hSelections->GetXaxis()->SetBinLabel(selPi0DCADau + 2, "#pi^{0} DCA dau.");
    hSelections->GetXaxis()->SetBinLabel(selXi0Radius + 2, "#Xi^{0} radius");
    hSelections->GetXaxis()->SetBinLabel(selXi0CosPA + 2, "#Xi^{0} cosPA");
    hSelections->GetXaxis()->SetBinLabel(selXi0DCADau + 2, "#Xi^{0} DCA dau.");
    hSelections->GetXaxis()->SetBinLabel(selXi0DCAxyToPV + 2, "#Xi^{0} DCA_{xy} to PV");
    hSelections->GetXaxis()->SetBinLabel(selXi0DCAzToPV + 2, "#Xi^{0} DCA_{z} to PV");
    hSelections->GetXaxis()->SetBinLabel(selXi0Lifetime + 2, "#Xi^{0} lifetime");
    hSelections->GetXaxis()->SetBinLabel(selXi0Rapidity + 2, "#Xi^{0} rapidity");
    hSelections->GetXaxis()->SetBinLabel(selXi0Pt + 2, "#Xi^{0} #it{p}_{T}");
    hSelections->GetXaxis()->SetBinLabel(selConsiderXi0 + 2, "True #Xi^{0}");
    hSelections->GetXaxis()->SetBinLabel(selConsiderAntiXi0 + 2, "True #bar{#Xi^{0}}");
    hSelections->GetXaxis()->SetBinLabel(selPhysPrimXi0 + 2, "Phys. prim. #Xi^{0}");
    hSelections->GetXaxis()->SetBinLabel(selPhysPrimAntiXi0 + 2, "Phys. prim. #bar{#Xi^{0}}");
    hSelections->GetXaxis()->SetBinLabel(selPhysPrimAntiXi0 + 3, "Cand. selected");

    //_______________________________________________
    // Main analysis output
    histos.add("h3dMassXi0", "h3dMassXi0", kTH3D, {axisCentrality, axisPt, axisXi0Mass});
    histos.add("h3dMassAntiXi0", "h3dMassAntiXi0", kTH3D, {axisCentrality, axisPt, axisXi0Mass});
    histos.add("h2dNbrOfXi0VsCentrality", "h2dNbrOfXi0VsCentrality", kTH2D, {axisCentrality, axisNCandidates});
    histos.add("h2dNbrOfAntiXi0VsCentrality", "h2dNbrOfAntiXi0VsCentrality", kTH2D, {axisCentrality, axisNCandidates});

    //_______________________________________________
    // Candidate QA, before and after the selections
    if (doCandidateQA) {
      for (int mode = 0; mode < NSelectionStages; mode++) {
        const std::string dir = (mode == 0) ? "QA/BeforeSel/" : "QA/AfterSel/";
        histos.add(dir + "h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisXi0Mass});
        histos.add(dir + "hMassPi0", "hMassPi0", kTH1D, {axisPi0Mass});
        histos.add(dir + "hMassLambda", "hMassLambda", kTH1D, {axisLambdaMass});
        histos.add(dir + "hPt", "hPt", kTH1D, {axisPt});
        histos.add(dir + "hRapidity", "hRapidity", kTH1D, {axisRapidity});
        histos.add(dir + "hCascRadius", "hCascRadius", kTH1D, {axisRadius});
        histos.add(dir + "hCascCosPA", "hCascCosPA", kTH1D, {axisCosPA});
        histos.add(dir + "hDCACascDaughters", "hDCACascDaughters", kTH1D, {axisDCADau});
        histos.add(dir + "hDCAxyCascToPV", "hDCAxyCascToPV", kTH1D, {axisDCAToPV});
        histos.add(dir + "hDCAzCascToPV", "hDCAzCascToPV", kTH1D, {axisDCAToPV});
        histos.add(dir + "hLifetime", "hLifetime", kTH1D, {axisLifetime});
        histos.add(dir + "hLambdaRadius", "hLambdaRadius", kTH1D, {axisRadius});
        histos.add(dir + "hLambdaCosPA", "hLambdaCosPA", kTH1D, {axisCosPA});
        histos.add(dir + "hPi0Radius", "hPi0Radius", kTH1D, {axisRadius});
        histos.add(dir + "hPi0CosPA", "hPi0CosPA", kTH1D, {axisCosPA});
        histos.add(dir + "hDCAPi0Daughters", "hDCAPi0Daughters", kTH1D, {axisDCADau});
        histos.add(dir + "h2dArmenterosLambda", "h2dArmenterosLambda", kTH2D, {axisArmAlpha, axisArmQt});
      }
    }

    //_______________________________________________
    // MC-specific histograms
    if (doprocessMonteCarlo && doMCQA) {
      histos.add("MCQA/h2dPtResolution", "h2dPtResolution", kTH2D, {axisPt, axisPtResolution});
      histos.add("MCQA/h2dPtRecoVsPtGen", "h2dPtRecoVsPtGen", kTH2D, {axisPt, axisPt});
      histos.add("MCQA/h3dMassXi0VsPtReco", "h3dMassXi0VsPtReco", kTH3D, {axisCentrality, axisPt, axisXi0Mass});
      histos.add("MCQA/h3dMassAntiXi0VsPtReco", "h3dMassAntiXi0VsPtReco", kTH3D, {axisCentrality, axisPt, axisXi0Mass});
      histos.add("MCQA/hXi0PDGCode", "hXi0PDGCode", kTH1D, {axisPDGCode});
      histos.add("MCQA/hXi0PDGCodeMother", "hXi0PDGCodeMother", kTH1D, {axisPDGCode});
      histos.add("MCQA/h2dPhoton1VsPhoton2PDGCode", "h2dPhoton1VsPhoton2PDGCode", kTH2D, {axisPDGCode, axisPDGCode});
      histos.add("MCQA/hLambdaPDGCode", "hLambdaPDGCode", kTH1D, {axisPDGCode});
    }

    //_______________________________________________
    // Generated-level histograms
    if (doprocessGenerated) {
      histos.add("Gen/hGenEvents", "hGenEvents", kTH2D, {axisNch, {2, -0.5f, +1.5f}});
      histos.add("Gen/hGenEventCentrality", "hGenEventCentrality", kTH1D, {axisCentrality});
      histos.add("Gen/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2D, {axisCentrality, {50, -0.5f, +49.5f}});
      histos.add("Gen/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2D, {axisCentrality, {50, -0.5f, +49.5f}});
      histos.add("Gen/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2D, {axisCentrality, axisNch});

      histos.add("Gen/h2dGenXi0", "h2dGenXi0", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenAntiXi0", "h2dGenAntiXi0", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenXi0VsMultMC", "h2dGenXi0VsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenAntiXi0VsMultMC", "h2dGenAntiXi0VsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenXi0_RecoedEvt", "h2dGenXi0_RecoedEvt", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenAntiXi0_RecoedEvt", "h2dGenAntiXi0_RecoedEvt", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/hGenXi0Rapidity", "hGenXi0Rapidity", kTH1D, {axisRapidity});
    }
  }

  //_______________________________________________
  // Check that all the bits of a mask are set in a bitmap
  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  //_______________________________________________
  // Centrality getter (Run 3)
  template <typename TCollision>
  float getCentralityRun3(TCollision const& collision)
  {
    if (centralityEstimator == kCentFT0C)
      return collision.centFT0C();
    else if (centralityEstimator == kCentFT0M)
      return collision.centFT0M();
    else if (centralityEstimator == kCentFT0CVariant1)
      return collision.centFT0CVariant1();
    else if (centralityEstimator == kCentMFT)
      return collision.centMFT();
    else if (centralityEstimator == kCentNGlobal)
      return collision.centNGlobal();
    else if (centralityEstimator == kCentFV0A)
      return collision.centFV0A();

    return -1.f;
  }

  //_______________________________________________
  // Check whether the collision passes our collision selections
  template <typename TCollision>
  bool isEventAccepted(TCollision const& collision, bool fillHists)
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds */);
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 18 /* Below min occupancy */);
    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 19 /* Above max occupancy */);

    // Fetch interaction rate only if required (in order to limit ccdb calls)
    float interactionRate = (fGetIR) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;
    if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 20 /* Below min IR */);
    if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 21 /* Above max IR */);

    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 22 /* selected collisions */);
      float centrality = getCentralityRun3(collision);
      histos.fill(HIST("hEventCentrality"), centrality);
      histos.fill(HIST("hEventPVz"), collision.posZ());
      if (doEventQA) {
        histos.fill(HIST("EventQA/hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
        histos.fill(HIST("EventQA/hCentralityVsPVz"), centrality, collision.posZ());
      }
    }

    return true;
  }

  //_______________________________________________
  // Photon (pi0 daughter) selection bitmap, computed on the de-referenced V0Cores entry.
  // The bitmaps of the two photons are AND-ed, so that a bit is set only if
  // both photons of the pi0 satisfy the corresponding criterion.
  template <typename TV0Object>
  uint64_t computePhotonBitmap(TV0Object const& gamma)
  {
    uint64_t bitMap = 0;

    if (photonSelections.v0TypeSelection < 0 || gamma.v0Type() == photonSelections.v0TypeSelection)
      BITSET(bitMap, selPhotonV0Type);

    if (gamma.mGamma() > 0 && gamma.mGamma() < photonSelections.maxMass)
      BITSET(bitMap, selPhotonMass);

    const float photonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);
    if (photonY > photonSelections.minRapidity && photonY < photonSelections.maxRapidity)
      BITSET(bitMap, selPhotonRapidity);

    if (std::abs(gamma.positiveeta()) < photonSelections.maxDauEta && std::abs(gamma.negativeeta()) < photonSelections.maxDauEta)
      BITSET(bitMap, selPhotonDauEta);

    if (std::abs(gamma.dcapostopv()) > photonSelections.minDCADauToPV && std::abs(gamma.dcanegtopv()) > photonSelections.minDCADauToPV)
      BITSET(bitMap, selPhotonDCADauToPV);

    if (std::abs(gamma.dcaV0daughters()) < photonSelections.maxDCADau)
      BITSET(bitMap, selPhotonDCADau);

    if (gamma.v0radius() > photonSelections.minRadius && gamma.v0radius() < photonSelections.maxRadius)
      BITSET(bitMap, selPhotonRadius);

    if (gamma.v0cosPA() > photonSelections.minCosPA)
      BITSET(bitMap, selPhotonCosPA);

    if (gamma.qtarm() < photonSelections.maxQt && std::abs(gamma.alpha()) < photonSelections.maxAlpha)
      BITSET(bitMap, selPhotonArmenteros);

    if (gamma.pt() > photonSelections.minPt && gamma.pt() < photonSelections.maxPt)
      BITSET(bitMap, selPhotonPt);

    auto posTrack = gamma.template posTrackExtra_as<DauTracks>();
    auto negTrack = gamma.template negTrackExtra_as<DauTracks>();

    if (posTrack.tpcCrossedRows() >= photonSelections.minTPCCrossedRows && negTrack.tpcCrossedRows() >= photonSelections.minTPCCrossedRows)
      BITSET(bitMap, selPhotonTPCCrossedRows);

    if (posTrack.tpcNSigmaEl() > photonSelections.minTPCNSigmaEl && posTrack.tpcNSigmaEl() < photonSelections.maxTPCNSigmaEl &&
        negTrack.tpcNSigmaEl() > photonSelections.minTPCNSigmaEl && negTrack.tpcNSigmaEl() < photonSelections.maxTPCNSigmaEl)
      BITSET(bitMap, selPhotonTPCPID);

    return bitMap;
  }

  //_______________________________________________
  // Lambda selection bitmap, computed on the de-referenced V0Cores entry
  template <typename TV0Object, typename TCollision>
  uint64_t computeLambdaBitmap(TV0Object const& lambda, TCollision const& collision)
  {
    uint64_t bitMap = 0;

    if (lambdaSelections.v0TypeSelection < 0 || lambda.v0Type() == lambdaSelections.v0TypeSelection)
      BITSET(bitMap, selLambdaV0Type);

    if (lambda.yLambda() > lambdaSelections.minRapidity && lambda.yLambda() < lambdaSelections.maxRapidity)
      BITSET(bitMap, selLambdaRapidity);

    if (std::abs(lambda.positiveeta()) < lambdaSelections.maxDauEta && std::abs(lambda.negativeeta()) < lambdaSelections.maxDauEta)
      BITSET(bitMap, selLambdaDauEta);

    if (std::abs(lambda.dcapostopv()) > lambdaSelections.minDCAPosToPV)
      BITSET(bitMap, selLambdaDCAPosToPV);

    if (std::abs(lambda.dcanegtopv()) > lambdaSelections.minDCANegToPV)
      BITSET(bitMap, selLambdaDCANegToPV);

    if (std::abs(lambda.dcaV0daughters()) < lambdaSelections.maxDCADau)
      BITSET(bitMap, selLambdaDCADau);

    if (lambda.v0radius() > lambdaSelections.minRadius && lambda.v0radius() < lambdaSelections.maxRadius)
      BITSET(bitMap, selLambdaRadius);

    if (lambda.v0cosPA() > lambdaSelections.minCosPA)
      BITSET(bitMap, selLambdaCosPA);

    const float lambdaLifetime = lambda.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (lambdaSelections.maxLifetime < 0 || lambdaLifetime < lambdaSelections.maxLifetime)
      BITSET(bitMap, selLambdaLifetime);

    // Lambda / anti-Lambda separation via the Armenteros-Podolanski alpha:
    // in the Lambda decay the baryon carries most of the longitudinal momentum
    if (lambda.alpha() > lambdaSelections.minAbsAlpha && lambda.alpha() < lambdaSelections.maxAbsAlpha)
      BITSET(bitMap, selLambdaArmenteros);
    if (lambda.alpha() < -lambdaSelections.minAbsAlpha && lambda.alpha() > -lambdaSelections.maxAbsAlpha)
      BITSET(bitMap, selAntiLambdaArmenteros);

    if (std::abs(lambda.mLambda() - o2::constants::physics::MassLambda0) < lambdaSelections.massWindow)
      BITSET(bitMap, selLambdaMass);
    if (std::abs(lambda.mAntiLambda() - o2::constants::physics::MassLambda0) < lambdaSelections.massWindow)
      BITSET(bitMap, selAntiLambdaMass);

    auto posTrack = lambda.template posTrackExtra_as<DauTracks>();
    auto negTrack = lambda.template negTrackExtra_as<DauTracks>();

    if (posTrack.tpcCrossedRows() >= lambdaSelections.minTPCCrossedRows && negTrack.tpcCrossedRows() >= lambdaSelections.minTPCCrossedRows)
      BITSET(bitMap, selLambdaTPCCrossedRows);

    if (lambdaSelections.minITSclusters < 0 ||
        (posTrack.itsNCls() >= lambdaSelections.minITSclusters && negTrack.itsNCls() >= lambdaSelections.minITSclusters))
      BITSET(bitMap, selLambdaITSClusters);

    // TPC PID, for both mass hypotheses
    if (std::abs(posTrack.tpcNSigmaPr()) < lambdaSelections.maxTPCNSigmaPr)
      BITSET(bitMap, selTPCPIDPositiveProton);
    if (std::abs(negTrack.tpcNSigmaPi()) < lambdaSelections.maxTPCNSigmaPi)
      BITSET(bitMap, selTPCPIDNegativePion);
    if (std::abs(negTrack.tpcNSigmaPr()) < lambdaSelections.maxTPCNSigmaPr)
      BITSET(bitMap, selTPCPIDNegativeProton);
    if (std::abs(posTrack.tpcNSigmaPi()) < lambdaSelections.maxTPCNSigmaPi)
      BITSET(bitMap, selTPCPIDPositivePion);

    return bitMap;
  }

  //_______________________________________________
  // pi0 and Xi0 (cascade) selection bitmap.
  // rapidity and pt are the generated ones when running over MC.
  template <typename TXi0Object, typename TCollision>
  uint64_t computeCascadeBitmap(TXi0Object const& xi0, TCollision const& collision, float rapidity, float pt)
  {
    uint64_t bitMap = 0;

    if (std::abs(xi0.pi0Mass() - o2::constants::physics::MassPi0) < xi0Selections.pi0MassWindow)
      BITSET(bitMap, selPi0Mass);

    if ((xi0Selections.minPi0Radius < 0 || xi0.radiusPi0() > xi0Selections.minPi0Radius) &&
        (xi0Selections.maxPi0Radius < 0 || xi0.radiusPi0() < xi0Selections.maxPi0Radius))
      BITSET(bitMap, selPi0Radius);

    if (xi0Selections.minPi0CosPA < -1 || xi0.pi0CosPA(collision.posX(), collision.posY(), collision.posZ()) > xi0Selections.minPi0CosPA)
      BITSET(bitMap, selPi0CosPA);

    if (xi0Selections.maxDCAPi0Daughters < 0 || xi0.dcadaughtersPi0() < xi0Selections.maxDCAPi0Daughters)
      BITSET(bitMap, selPi0DCADau);

    if (xi0.radius() > xi0Selections.minCascRadius && xi0.radius() < xi0Selections.maxCascRadius)
      BITSET(bitMap, selXi0Radius);

    if (xi0.cascCosPA(collision.posX(), collision.posY(), collision.posZ()) > xi0Selections.minCascCosPA)
      BITSET(bitMap, selXi0CosPA);

    if (xi0.dcadaughters() < xi0Selections.maxDCACascDaughters)
      BITSET(bitMap, selXi0DCADau);

    if (xi0Selections.maxDCAxyCascToPV < 0 || std::abs(xi0.dcaXYCascToPV()) < xi0Selections.maxDCAxyCascToPV)
      BITSET(bitMap, selXi0DCAxyToPV);

    if (xi0Selections.maxDCAzCascToPV < 0 || std::abs(xi0.dcaZCascToPV()) < xi0Selections.maxDCAzCascToPV)
      BITSET(bitMap, selXi0DCAzToPV);

    if (xi0Selections.maxLifetime < 0 || getProperLifetime(xi0, collision) < xi0Selections.maxLifetime * CtauXi0)
      BITSET(bitMap, selXi0Lifetime);

    if (rapidity > xi0Selections.minRapidity && rapidity < xi0Selections.maxRapidity)
      BITSET(bitMap, selXi0Rapidity);

    if (pt > xi0Selections.minPt && pt < xi0Selections.maxPt)
      BITSET(bitMap, selXi0Pt);

    return bitMap;
  }

  //_______________________________________________
  // MC association bitmap.
  // The builder sets the Xi0 PDG code only if the full decay chain is matched,
  // i.e. two true photons from the same true pi0 and a true Lambda sharing the
  // same true Xi0 mother, so a single PDG check is enough here.
  template <typename TXi0Object>
  uint64_t computeMCAssociation(TXi0Object const& xi0)
  {
    uint64_t bitMap = 0;

    if (xi0.pdgCode() == o2::constants::physics::Pdg::kXi0) {
      BITSET(bitMap, selConsiderXi0);
      if (xi0.isPhysicalPrimary())
        BITSET(bitMap, selPhysPrimXi0);
    }
    if (xi0.pdgCode() == -o2::constants::physics::Pdg::kXi0) {
      BITSET(bitMap, selConsiderAntiXi0);
      if (xi0.isPhysicalPrimary())
        BITSET(bitMap, selPhysPrimAntiXi0);
    }

    return bitMap;
  }

  //_______________________________________________
  // Xi0 proper lifetime, m*L/p, in cm
  template <typename TXi0Object, typename TCollision>
  float getProperLifetime(TXi0Object const& xi0, TCollision const& collision)
  {
    const float decayLength = RecoDecay::distance(std::array{collision.posX(), collision.posY(), collision.posZ()},
                                                  std::array{xi0.x(), xi0.y(), xi0.z()});
    const float totalMomentum = xi0.p();
    return (totalMomentum > 0) ? o2::constants::physics::MassXi0 * decayLength / totalMomentum : 1e6f;
  }

  //_______________________________________________
  // Fill the candidate QA histograms. mode = 0: before selections, 1: after selections
  template <int mode, typename TXi0Object, typename TV0Object, typename TCollision>
  void fillCandidateQA(TXi0Object const& xi0, TV0Object const& lambda, TCollision const& collision, float pt, float rapidity, float centrality, float massLambda)
  {
    static constexpr std::string_view MainDir[] = {"QA/BeforeSel", "QA/AfterSel"};

    histos.fill(HIST(MainDir[mode]) + HIST("/h3dMass"), centrality, pt, xi0.xi0Mass());
    histos.fill(HIST(MainDir[mode]) + HIST("/hMassPi0"), xi0.pi0Mass());
    histos.fill(HIST(MainDir[mode]) + HIST("/hMassLambda"), massLambda);
    histos.fill(HIST(MainDir[mode]) + HIST("/hPt"), pt);
    histos.fill(HIST(MainDir[mode]) + HIST("/hRapidity"), rapidity);
    histos.fill(HIST(MainDir[mode]) + HIST("/hCascRadius"), xi0.radius());
    histos.fill(HIST(MainDir[mode]) + HIST("/hCascCosPA"), xi0.cascCosPA(collision.posX(), collision.posY(), collision.posZ()));
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCACascDaughters"), xi0.dcadaughters());
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCAxyCascToPV"), xi0.dcaXYCascToPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCAzCascToPV"), xi0.dcaZCascToPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/hLifetime"), getProperLifetime(xi0, collision));
    histos.fill(HIST(MainDir[mode]) + HIST("/hLambdaRadius"), xi0.radiusLambda());
    histos.fill(HIST(MainDir[mode]) + HIST("/hLambdaCosPA"), xi0.lambdaCosPA(collision.posX(), collision.posY(), collision.posZ()));
    histos.fill(HIST(MainDir[mode]) + HIST("/hPi0Radius"), xi0.radiusPi0());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPi0CosPA"), xi0.pi0CosPA(collision.posX(), collision.posY(), collision.posZ()));
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCAPi0Daughters"), xi0.dcadaughtersPi0());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenterosLambda"), lambda.alpha(), lambda.qtarm());
  }

  //_______________________________________________
  // Main reconstructed-level analysis function, common to data and MC
  template <typename TCollisions, typename TXi0s, typename TV0s>
  void analyseRecoedXi0s(TCollisions const& collisions, TXi0s const& fullXi0s, TV0s const& fullV0s)
  {
    // Custom grouping: the Xi0 candidates are not sorted per collision
    std::vector<std::vector<int>> xi0grouped(collisions.size());
    for (const auto& xi0 : fullXi0s) {
      xi0grouped[xi0.straCollisionId()].push_back(xi0.globalIndex());
    }

    const int64_t nV0s = fullV0s.size();

    for (const auto& coll : collisions) {
      // Event selection
      if (!isEventAccepted(coll, true)) {
        continue;
      }
      const float centrality = getCentralityRun3(coll);

      int nXi0s = 0;
      int nAntiXi0s = 0;

      // Xi0 candidates loop
      for (size_t i = 0; i < xi0grouped[coll.globalIndex()].size(); i++) {
        auto xi0 = fullXi0s.rawIteratorAt(xi0grouped[coll.globalIndex()][i]);

        histos.fill(HIST("hSelections"), 0 /* all candidates */);

        //_______________________________________________
        // De-reference the three daughter V0s
        if (xi0.photon1Index() < 0 || xi0.photon1Index() >= nV0s ||
            xi0.photon2Index() < 0 || xi0.photon2Index() >= nV0s ||
            xi0.lambdaIndex() < 0 || xi0.lambdaIndex() >= nV0s) {
          continue;
        }
        auto photon1 = fullV0s.rawIteratorAt(xi0.photon1Index());
        auto photon2 = fullV0s.rawIteratorAt(xi0.photon2Index());
        auto lambda = fullV0s.rawIteratorAt(xi0.lambdaIndex());

        //_______________________________________________
        // Kinematics entering both the selections and the histograms.
        // Over MC, the generated rapidity and pT are used instead of the
        // reconstructed ones, so that the acceptance x efficiency numerator
        // and denominator are defined identically.
        const float ptRecoed = xi0.pt();
        float pt = ptRecoed;
        float rapidity = xi0.rapidity();
        bool hasMCInfo = false;
        if constexpr (requires { xi0.pdgCode(); }) {
          // the MC information is available only if the three daughter V0s
          // could be associated to a V0MCCore entry by the builder
          hasMCInfo = (xi0.photon1PDGCode() != 0 && xi0.photon2PDGCode() != 0 && xi0.lambdaPDGCode() != 0);
          if (hasMCInfo) {
            pt = xi0.mcpt();
            rapidity = xi0.rapidityMC();
          }
        }

        //_______________________________________________
        // Assemble the selection bitmap
        uint64_t selMap = computePhotonBitmap(photon1) & computePhotonBitmap(photon2);
        selMap |= computeLambdaBitmap(lambda, coll);
        selMap |= computeCascadeBitmap(xi0, coll, rapidity, pt);

        if constexpr (requires { xi0.pdgCode(); }) {
          selMap |= computeMCAssociation(xi0);

          if (doMCQA) {
            histos.fill(HIST("MCQA/hXi0PDGCode"), xi0.pdgCode());
            histos.fill(HIST("MCQA/hXi0PDGCodeMother"), xi0.pdgCodeMother());
            histos.fill(HIST("MCQA/h2dPhoton1VsPhoton2PDGCode"), xi0.photon1PDGCode(), xi0.photon2PDGCode());
            histos.fill(HIST("MCQA/hLambdaPDGCode"), xi0.lambdaPDGCode());
          }

          // disregard the MC association if not explicitly asked for
          if (!mcSelections.doMCAssociation) {
            BITSET(selMap, selConsiderXi0);
            BITSET(selMap, selConsiderAntiXi0);
            BITSET(selMap, selPhysPrimXi0);
            BITSET(selMap, selPhysPrimAntiXi0);
          }
        } else {
          // real data: every candidate is considered for both hypotheses,
          // the Armenteros bits alone tell Xi0 and anti-Xi0 apart
          BITSET(selMap, selConsiderXi0);
          BITSET(selMap, selConsiderAntiXi0);
          BITSET(selMap, selPhysPrimXi0);
          BITSET(selMap, selPhysPrimAntiXi0);
        }

        // Selection bookkeeping: one bin per criterion
        for (int bit = 0; bit <= static_cast<int>(selPhysPrimAntiXi0); bit++) {
          if (BITCHECK(selMap, bit)) {
            histos.fill(HIST("hSelections"), bit + 1);
          }
        }

        const bool passXi0 = verifyMask(selMap, maskSelectionXi0);
        const bool passAntiXi0 = verifyMask(selMap, maskSelectionAntiXi0);
        // the sign of the Armenteros alpha already tells the two hypotheses apart
        const float massLambda = (lambda.alpha() < 0) ? lambda.mAntiLambda() : lambda.mLambda();

        if (doCandidateQA) {
          fillCandidateQA<0>(xi0, lambda, coll, pt, rapidity, centrality, massLambda);
        }

        if (!passXi0 && !passAntiXi0) {
          continue;
        }
        histos.fill(HIST("hSelections"), static_cast<int>(selPhysPrimAntiXi0) + 2 /* candidate selected */);

        if (doCandidateQA) {
          fillCandidateQA<1>(xi0, lambda, coll, pt, rapidity, centrality, massLambda);
        }

        if (passXi0) {
          nXi0s++;
          histos.fill(HIST("h3dMassXi0"), centrality, pt, xi0.xi0Mass());
        }
        if (passAntiXi0) {
          nAntiXi0s++;
          histos.fill(HIST("h3dMassAntiXi0"), centrality, pt, xi0.xi0Mass());
        }

        //_______________________________________________
        // MC-specific QA
        if constexpr (requires { xi0.pdgCode(); }) {
          if (doMCQA) {
            if (passXi0) {
              histos.fill(HIST("MCQA/h3dMassXi0VsPtReco"), centrality, ptRecoed, xi0.xi0Mass());
            }
            if (passAntiXi0) {
              histos.fill(HIST("MCQA/h3dMassAntiXi0VsPtReco"), centrality, ptRecoed, xi0.xi0Mass());
            }
            if (hasMCInfo && pt > 0) {
              histos.fill(HIST("MCQA/h2dPtResolution"), pt, (ptRecoed - pt) / pt);
              histos.fill(HIST("MCQA/h2dPtRecoVsPtGen"), pt, ptRecoed);
            }
          }
        }
      }

      histos.fill(HIST("h2dNbrOfXi0VsCentrality"), centrality, nXi0s);
      histos.fill(HIST("h2dNbrOfAntiXi0VsCentrality"), centrality, nAntiXi0s);
    }
  }

  //_______________________________________________
  // Generated-level processing
  // Return the list of indices of the reconstructed collision associated to each MC collision
  template <typename TMCCollisions, typename TCollisions>
  std::vector<int> getListOfRecoCollIndices(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size(), -1);
    std::vector<std::vector<int>> groupedCollisions(mcCollisions.size());

    for (const auto& coll : collisions) {
      if (coll.straMCCollisionId() < 0) {
        continue;
      }
      groupedCollisions[coll.straMCCollisionId()].push_back(coll.globalIndex());
    }

    for (auto const& mcCollision : mcCollisions) {
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (size_t i = 0; i < groupedCollisions[mcCollision.globalIndex()].size(); i++) {
        auto collision = collisions.rawIteratorAt(groupedCollisions[mcCollision.globalIndex()][i]);

        if (eventSelections.useEvtSelInDenomEff && !isEventAccepted(collision, false)) {
          continue;
        }
        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
        }
      }
      listBestCollisionIdx[mcCollision.globalIndex()] = bestCollisionIndex;
    }
    return listBestCollisionIdx;
  }

  //_______________________________________________
  // Event selection applied at generated level
  template <typename TMCCollision>
  bool isGeneratedEventAccepted(TMCCollision const& mcCollision)
  {
    if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
      return false;
    }
    if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
      return false;
    }
    return true;
  }

  //_______________________________________________
  // Generated-level processing: generated event properties (event loss/splitting)
  template <typename TMCCollisions, typename TCollisions>
  void fillGeneratedEventProperties(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<std::vector<int>> groupedCollisions(mcCollisions.size());
    for (const auto& coll : collisions) {
      if (coll.straMCCollisionId() < 0) {
        continue;
      }
      groupedCollisions[coll.straMCCollisionId()].push_back(coll.globalIndex());
    }

    for (auto const& mcCollision : mcCollisions) {
      if (!isGeneratedEventAccepted(mcCollision)) {
        continue;
      }

      histos.fill(HIST("Gen/hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events */);

      bool atLeastOne = false;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (size_t i = 0; i < groupedCollisions[mcCollision.globalIndex()].size(); i++) {
        auto collision = collisions.rawIteratorAt(groupedCollisions[mcCollision.globalIndex()][i]);

        if (!isEventAccepted(collision, false)) {
          continue;
        }
        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = getCentralityRun3(collision);
        }
        nCollisions++;
        atLeastOne = true;
      }

      histos.fill(HIST("Gen/hCentralityVsNcoll_beforeEvSel"), centrality, static_cast<int>(groupedCollisions[mcCollision.globalIndex()].size()));
      histos.fill(HIST("Gen/hCentralityVsNcoll_afterEvSel"), centrality, nCollisions);
      histos.fill(HIST("Gen/hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());

      if (atLeastOne) {
        histos.fill(HIST("Gen/hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* at least 1 rec. event */);
        histos.fill(HIST("Gen/hGenEventCentrality"), centrality);
      }
    }
  }

  //_______________________________________________
  // Main generated-level analysis function.
  // The generated Xi0 are read from the CascMCCores table.
  template <typename TMCCollisions, typename TCollisions, typename TCascMCs>
  void analyseGeneratedXi0s(TMCCollisions const& mcCollisions, TCollisions const& collisions, TCascMCs const& cascMCCores)
  {
    fillGeneratedEventProperties(mcCollisions, collisions);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions);

    for (auto const& cascMC : cascMCCores) {
      if (std::abs(cascMC.pdgCode()) != o2::constants::physics::Pdg::kXi0) {
        continue;
      }
      if (!cascMC.has_straMCCollision()) {
        continue;
      }
      if (mcSelections.requirePhysicalPrimary && !cascMC.isPhysicalPrimary()) {
        continue;
      }

      // N.B.: cascdata::RapidityMC only knows about the Xi- and Omega- masses,
      //       so the Xi0 rapidity is computed here explicitly
      const float ptMC = cascMC.ptMC();
      const float rapidityMC = RecoDecay::y(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}, o2::constants::physics::MassXi0);

      histos.fill(HIST("Gen/hGenXi0Rapidity"), rapidityMC);

      // same rapidity window as at reconstructed level
      if (rapidityMC < xi0Selections.minRapidity || rapidityMC > xi0Selections.maxRapidity) {
        continue;
      }

      auto mcCollision = cascMC.template straMCCollision_as<TMCCollisions>();
      if (!isGeneratedEventAccepted(mcCollision)) {
        continue;
      }

      float centrality = 100.5f;
      const int bestCollisionIdx = listBestCollisionIdx[mcCollision.globalIndex()];
      if (bestCollisionIdx > -1) {
        auto collision = collisions.rawIteratorAt(bestCollisionIdx);
        centrality = getCentralityRun3(collision);

        if (cascMC.pdgCode() > 0) {
          histos.fill(HIST("Gen/h2dGenXi0_RecoedEvt"), centrality, ptMC);
        } else {
          histos.fill(HIST("Gen/h2dGenAntiXi0_RecoedEvt"), centrality, ptMC);
        }
      }

      if (cascMC.pdgCode() > 0) {
        histos.fill(HIST("Gen/h2dGenXi0"), centrality, ptMC);
        histos.fill(HIST("Gen/h2dGenXi0VsMultMC"), mcCollision.multMCNParticlesEta05(), ptMC);
      } else {
        histos.fill(HIST("Gen/h2dGenAntiXi0"), centrality, ptMC);
        histos.fill(HIST("Gen/h2dGenAntiXi0VsMultMC"), mcCollision.multMCNParticlesEta05(), ptMC);
      }
    }
  }

  //_______________________________________________
  // Process functions
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps> const& collisions,
                       Xi0Candidates const& fullXi0s,
                       V0Candidates const& fullV0s,
                       DauTracks const&)
  {
    analyseRecoedXi0s(collisions, fullXi0s, fullV0s);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps, aod::StraCollLabels> const& collisions,
                         Xi0McCandidates const& fullXi0s,
                         V0Candidates const& fullV0s,
                         DauTracks const&)
  {
    analyseRecoedXi0s(collisions, fullXi0s, fullV0s);
  }

  void processGenerated(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions,
                        soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps, aod::StraCollLabels> const& collisions,
                        soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& cascMCCores)
  {
    analyseGeneratedXi0s(mcCollisions, collisions, cascMCCores);
  }

  PROCESS_SWITCH(DerivedXi0Analysis, processRealData, "process as if real data", true);
  PROCESS_SWITCH(DerivedXi0Analysis, processMonteCarlo, "process reconstructed information in MC", false);
  PROCESS_SWITCH(DerivedXi0Analysis, processGenerated, "process pure generated information in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DerivedXi0Analysis>(cfgc)};
}
