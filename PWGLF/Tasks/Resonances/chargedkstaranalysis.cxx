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
/// \brief this is a code for the CKS resonance
/// \author prottay das
/// \since 13/01/2024

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include "TF1.h"

#include <array>
#include <cmath>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

// For charged kstarpp analysis
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct chargedkstaranalysis {

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{
    "ccdb-no-later-than",
    std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch())
      .count(),
    "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080",
                                "url of the ccdb repository"};

  SliceCache cache;

  // For charged Kstarpp analysis use Resonance Initalizer and THnSparse
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Binning of the centrality axis"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis Etabins{"Etabins", {VARIABLE_WIDTH, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, "Eta Binning"};
  Configurable<int> cDCABinsQA{"cDCABinsQA", 150, "DCA binning"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};

  HistogramRegistry histos1{"histos1", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  HistogramRegistry rGenParticles{"genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rRecParticles{"recParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  /// PID Selections
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined

  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contributor
  // V0 selections
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.97, "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 1.0, "V0 daughter DCA Maximum"};
  // Competing V0 rejection
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.0043, "Mass window for competing Lambda0 rejection"};
  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};

  // Event mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};
  Configurable<int> cTpcNsigmaPionBinsQA{"cTpcNsigmaPionBinsQA", 140, "tpcNSigmaPi binning"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Confugrable for QA histograms
  Configurable<bool> QAbefore{"QAbefore", false, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", false, "QAafter"};
  Configurable<bool> QAv0{"QAv0", false, "QAv0"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f,
                                 "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f,
                                  "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.0f,
                                        "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.985f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 0.5f,
                                         "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 200.f,
                                         "Maximum transverse radius"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15,
                                      "Maximum V0 life time"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "DCA V0 to PV"};
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4,
                                     "n Sigma cut on KS0 mass"};
  Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};

  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f,
                                   "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f,
                                          "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{
    "ConfDaughDCAMin", 0.06f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 4,
                                       "PID selections for KS0 daughters"};

  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f,
                                  "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0,
                                   "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0,
                                   "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0,
                                        "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5,
                                     "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "cfgMultFT0"};
  Configurable<bool> cfgCentFT0C{"cfgCentFT0C", true, "cfgCentFT0C"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isMC{"isMC", true, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> TVXEvsel{"TVXEvsel", false, "Triggger selection"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    AxisSpec dcaxyAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec tpcNSigmaPiAxisQA = {cTpcNsigmaPionBinsQA, -7.0, 7.0, "N#sigma_{TPC}"};

    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec etaAxis = {Etabins, "#eta"};
    AxisSpec goodTrackCountAxis = {3, 0., 3., "Passed track = 1, Passed V0 = 2, Passed track and V0 = 3"};
    // register histograms
    histos1.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos1.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    // Multiplicity and accepted events QA
    histos1.add("QAbefore/collMult", "Collision multiplicity", HistType::kTH1F, {centAxis});
    // QA before
    histos1.add("QAbefore/pi_Eta", "Primary pion track eta", kTH1F, {etaAxis});
    histos1.add("QAbefore/k0s_Eta", "K0short track eta", kTH1F, {etaAxis});
    histos1.add("QAbefore/chargedkstarpmRapidity", "Reconstructed K*^{#pm} rapidity", kTH1F, {etaAxis});

    histos1.add("QAbefore/DCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxisQA});
    histos1.add("QAbefore/DCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxisQA});
    histos1.add("QAbefore/pT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos1.add("QAbefore/tpcNsigmaPionQA", "NsigmaTPC distribution of primary pion candidates", kTH2F, {ptAxisQA, tpcNSigmaPiAxisQA});

    // QA after
    histos1.add("QAAfter/DCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxisQA});
    histos1.add("QAAfter/DCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxisQA});
    histos1.add("QAAfter/pT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos1.add("QAAfter/tpcNsigmaPionQA", "NsigmaTPC distribution of primary pion candidates", kTH2F, {ptAxisQA, tpcNSigmaPiAxisQA});
    histos1.add("QAAfter/pi_Eta", "Primary pion track eta", kTH1F, {etaAxis});

    // Good tracks and V0 counts QA
    histos1.add("QAafter/hGoodTracksV0s", "Number of good track and V0 passed", kTH1F, {goodTrackCountAxis});
    histos1.add("chargedkstarinvmassUlikeSign", "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});
    histos1.add("chargedkstarinvmassMixedEvent", "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});

    // Mass vs Pt vs Multiplicity 3-dimensional histogram
    //    histos1.add("chargekstarMassPtMult", "Charged K*(892) mass vs pT vs V0 multiplicity distribution", kTH3F, {invMassAxis, ptAxis, centAxis});

    histos1.add("chargekstarMassPtMultPtUnlikeSign",
                "Invariant mass of CKS meson Unlike Sign", kTHnSparseF,
                {invMassAxis, ptAxis, centAxis}, true);
    histos1.add("chargekstarMassPtMultPtMixedEvent",
                "Invariant mass of CKS meson MixedEvent Sign", kTHnSparseF,
                {invMassAxis, ptAxis, centAxis}, true);

    // Axes
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f,
                                "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -10., 10., "vrtx_{Z} [cm]"};
    AxisSpec multAxis = {100, 0.0f, 100.0f, "Multiplicity"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec",
                        {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F,
                        {{200, 0.0f, 200.0f}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaPionTPC_before", "NsigmaPion TPC distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_before", "NsigmaPion TOF distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});

      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaxy_after", "Dcaxy distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hDcaz_after", "Dcaz distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
    }

    if (QAv0) {
      // K0s reconstruction
      histos.add(
        "hMassvsptvsmult", "hMassvsptvsmult",
        {HistType::kTHnSparseF, {{K0ShortMassAxis}, {ptAxis}, {multAxis}}},
        true);
      // K0s topological/PID cuts
      histos.add("hDCAV0Daughters", "hDCAV0Daughters",
                 {HistType::kTH1F, {{50, 0.0f, 5.0f}}});
      histos.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      histos.add("hV0CosPA", "hV0CosPA",
                 {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    }

    // histos.add("counter", "counter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});

    // CKStar histograms
    histos.add("h3CKSInvMassUnlikeSign",
               "Invariant mass of CKS meson Unlike Sign", kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);
    histos.add("h3CKSInvMassMixed", "Invariant mass of CKS meson Mixed",
               kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);

    if (isMC) {
      rGenParticles.add("hMC", "Gen MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      rGenParticles.add("hCentGen", "Gen MC Event centrality", kTH1F, {{100, 0.0f, 100.0f}});
      rRecParticles.add("hCentRec", "Rec MC Event centrality", kTH1F, {{100, 0.0f, 100.0f}});
      rRecParticles.add("hMCRec", "Rec MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      rGenParticles.add("hPtK0ShortGen", "hPtK0ShortGen", {HistType::kTH1F, {{ptAxis}}});
      // rGenParticles.add("hCKSGen", "hCKSGen", {HistType::kTH1F, {{ptAxis}}});
      // rRecParticles.add("hCKSRec", "hCKSRec", {HistType::kTH1F, {{ptAxis}}});

      rGenParticles.add("hCKSGen",
                        "Invariant mass of CKS meson Gen", kTHnSparseF,
                        {{200, 0.0, 20.0}, {100, 0.0f, 100.0f}}, true);
      rRecParticles.add("hCKSRec",
                        "Invariant mass of CKS meson Rec", kTHnSparseF,
                        {{200, 0.0, 20.0}, {100, 0.0f, 100.0f}}, true);
    }

    // Event selection cut additional - Alex
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
  }

  double massPi = TDatabasePDG::Instance()
                    ->GetParticle(kPiPlus)
                    ->Mass();
  double massK0s = TDatabasePDG::Instance()
                     ->GetParticle(kK0Short)
                     ->Mass();
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector CKSVector;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      // return 0;
    }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    // if (multTrk < fMultCutLow->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;

    return 1;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut &&
        (!candidate.isGlobalTrack() || !candidate.isPVContributor() ||
         candidate.itsNCls() < cfgITScluster)) {
      return false;
    }
    /*
    if (ismanualDCAcut &&
        (!candidate.isGlobalTrackWoDCA() || !candidate.isPVContributor() ||
         std::abs(candidate.dcaXY()) > cfgCutDCAxy ||
         std::abs(candidate.dcaZ()) > cfgCutDCAz ||
         candidate.itsNCls() < cfgITScluster || candidate.tpcNClsFound() < 70)) {
      return false;
    }
    */
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {

    if (candidate.hasTOF() &&
        (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() +
         candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) <
          (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }

    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    }

    /*
    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    }
    else if (candidate.hasTOF() &&
        std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
      return true;
      }
    */

    /*
    if (candidate.hasTOF() &&
        (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() +
         candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) <
          (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    */

    /*
    if (!isNoTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    */

    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate,
                   float multiplicity)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    if (TMath::Abs(candidate.yK0Short()) > 0.5) {
      return false;
    }

    const float qtarm = candidate.qtarm();
    const float alph = candidate.alpha();
    float arm = qtarm / alph;
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(),
                                             collision.posZ()) *
                    TDatabasePDG::Instance()
                      ->GetParticle(kK0Short)
                      ->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - cWidthKs0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + cWidthKs0 * cSigmaMassKs0;
    // float decayLength = candidate.distovertotmom(collision.posX(),
    // collision.posY(), collision.posZ()) *
    // RecoDecay::sqrtSumOfSquares(candidate.px(), candidate.py(),
    // candidate.pz());

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (fabs(CtauK0s) > cMaxV0LifeTime ||
        candidate.mK0Short() < lowmasscutks0 ||
        candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    if (arm < 0.2) {
      return false;
    }

    if (QAv0) {
      histos.fill(HIST("hLT"), CtauK0s);
      histos.fill(HIST("hMassvsptvsmult"), candidate.mK0Short(), candidate.pt(),
                  multiplicity);
      histos.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      histos.fill(HIST("hV0CosPA"), candidate.v0cosPA());
    }
    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge,
                            double nsigmaV0Daughter)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();

    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < 0.8)
      return false;

    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }
    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
      return false;
    }

    return true;
  }

  double massK0 = o2::constants::physics::MassK0Short;
  double massPicharged = o2::constants::physics::MassPionCharged;
  double massLambda0 = o2::constants::physics::MassLambda;
  double massAntiLambda0 = o2::constants::physics::MassLambda0Bar;
  // Fill histograms (main function)
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename V0sType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks, const V0sType& dV0s)
  {
    // auto multiplicity = collision.cent();
    auto multiplicity = collision.cent();
    histos1.fill(HIST("QAbefore/collMult"), multiplicity);
    TLorentzVector lDecayDaughter, lDecayV0, lResonance;

    for (auto track : dTracks) { // loop over all dTracks1
      // if (!trackCut(track1))
      // continue; // track selection and PID selection
      // trying to see the information without applying any cut yet //Let's I am trying to reconstruct the charged kstar it is V0s + pion
      histos1.fill(HIST("hEta"), track.eta());

      auto trackId = track.index();
      auto trackptPi = track.pt();
      auto tracketaPi = track.eta();

      histos1.fill(HIST("QAbefore/pi_Eta"), tracketaPi);

      if (!IsMix) {
        // DCA QA (before cuts)
        histos1.fill(HIST("QAbefore/DCAxy_pi"), track.dcaXY());
        histos1.fill(HIST("QAbefore/DCAz_pi"), track.dcaZ());
        // Pseudo-rapidity QA (before cuts)
        histos1.fill(HIST("QAbefore/pi_Eta"), tracketaPi);
        // pT QA (before cuts)
        histos1.fill(HIST("QAbefore/pT_pi"), trackptPi);
        // TPC PID (before cuts)
        histos1.fill(HIST("QAbefore/tpcNsigmaPionQA"), trackptPi, track.tpcNSigmaPi());
      }

      // apply the track cut
      if (!trackCutpp(track) || !selectionPIDpp(track))
        continue;

      histos1.fill(HIST("QAafter/hGoodTracksV0s"), 0.5);

      if (!IsMix) {
        // DCA QA (before cuts)
        histos1.fill(HIST("QAAfter/DCAxy_pi"), track.dcaXY());
        histos1.fill(HIST("QAAfter/DCAz_pi"), track.dcaZ());
        // Pseudo-rapidity QA (before cuts)
        histos1.fill(HIST("QAAfter/pi_Eta"), tracketaPi);
        // pT QA (before cuts)
        histos1.fill(HIST("QAAfter/pT_pi"), trackptPi);
        // TPC PID (before cuts)
        histos1.fill(HIST("QAAfter/tpcNsigmaPionQA"), trackptPi, track.tpcNSigmaPi());
      }

      for (auto& v0 : dV0s) {

        // Full index policy is needed to consider all possible combinations
        if (v0.indices()[0] == trackId || v0.indices()[1] == trackId)
          continue; // To avoid combining secondary and primary pions
                    //// Initialize variables
        // trk: Pion, v0: K0s
        // apply the track cut
        if (!V0Cut(v0))
          continue;
        histos1.fill(HIST("QAafter/hGoodTracksV0s"), 1.5);

        lDecayDaughter.SetXYZM(track.px(), track.py(), track.pz(), massPi);
        lDecayV0.SetXYZM(v0.px(), v0.py(), v0.pz(), massK0);
        lResonance = lDecayDaughter + lDecayV0;
        // Counting how many resonances passed
        histos1.fill(HIST("QAafter/hGoodTracksV0s"), 2.5);

        // Checking whether the mid-rapidity condition is met
        if (abs(lResonance.Rapidity()) > 0.5)
          continue;
        if constexpr (!IsMix) {
          histos1.fill(HIST("chargedkstarinvmassUlikeSign"), lResonance.M());
          // Reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos1.fill(HIST("chargekstarMassPtMultPtUnlikeSign"), lResonance.M(), lResonance.Pt(), multiplicity);

        } else {
          histos1.fill(HIST("chargedkstarinvmassMixedEvent"), lResonance.M());
          // Reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos1.fill(HIST("chargekstarMassPtMultPtMixedEvent"), lResonance.M(), lResonance.Pt(), multiplicity);
        }
      }
    }
  }

  template <typename T>
  bool selectionPIDpp(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      tpcPIDPassed = true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        tofPIDPassed = true;
      }
      if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <typename TrackType>
  bool trackCutpp(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.eta()) > ConfDaughEta)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }
  template <typename V0Type>
  bool V0Cut(const V0Type v0)
  {
    // V0 track cuts
    if (std::abs(v0.eta()) > ConfDaughEta)
      return false;
    if (v0.v0CosPA() < cV0MinCosPA)
      return false;
    if (v0.daughDCA() > cV0MaxDaughDCA)
      return false;

    // apply the competing V0 rejection cut (excluding Lambda0 candidates, massLambdaPDG = 1115.683 MeV/c2)

    if (std::abs(v0.mLambda() - massLambda0) < cV0MassWindow)
      return false;
    if (std::abs(v0.mAntiLambda() - massAntiLambda0) < cV0MassWindow)
      return false;

    return true;
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  Filter acceptanceFilter =
    (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) &&
                        (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  // using EventCandidates = soa::Filtered<
  // soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs,
  // aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                                  aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  using TrackCandidatesMC = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  using V0TrackCandidatesMC = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using V0TrackCandidate = aod::V0Datas;

  ConfigurableAxis axisVertex{
    "axisVertex",
    {20, -10, 10},
    "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{
    "axisMultiplicityClass",
    {2, 0, 100},
    "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{
    "axisMultiplicity",
    {2000, 0, 10000},
    "TPC multiplicity  for bin"};

  using BinningTypeTPCMultiplicity =
    ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // using BinningTypeVertexContributor =
  // ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  using BinningTypeCentralityM =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  BinningTypeVertexContributor binningOnPositions{
    {axisVertex, axisMultiplicityClass},
    true};

  Pair<EventCandidates, TrackCandidates, V0TrackCandidate,
       BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  /*
  SameKindPair<EventCandidates, TrackCandidates,
       BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  */

  void processSE(EventCandidates::iterator const& collision,
                 TrackCandidates const& tracks, aod::V0Datas const& V0s,
                 aod::BCs const&)

  /*
  void processSE(EventCandidates::iterator const& collision,
     TrackCandidates const& tracks, aod::BCs const&)
  */
  {

    if (!collision.sel8()) {
      return;
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> pions, kaons, kshorts;
    std::vector<ROOT::Math::PtEtaPhiMVector> pions2, kaons2;
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> PionSign = {};
    std::vector<int64_t> PioncollIndex = {};
    std::vector<int64_t> PionIndex2 = {};
    std::vector<int64_t> PionSign2 = {};
    std::vector<int64_t> PioncollIndex2 = {};

    std::vector<int64_t> KaonIndex = {};
    std::vector<int64_t> KaonSign = {};
    std::vector<int64_t> KaoncollIndex = {};
    std::vector<int64_t> KaonIndex2 = {};
    std::vector<int64_t> KaonSign2 = {};
    std::vector<int64_t> KaoncollIndex2 = {};

    std::vector<int64_t> V0collIndex = {};
    std::vector<int64_t> KshortPosDaughIndex = {};
    std::vector<int64_t> KshortNegDaughIndex = {};
    /*
    float multiplicity = 0.0f;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
      multiplicity = collision.centFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 0)
      multiplicity = collision.centFT0M();
    */
    float centrality = 0.0f;
    centrality = collision.centFT0C();

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }

    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), centrality);

    for (auto track1 : tracks) {
      /*
      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track1.tofNSigmaPi());
      }
      */
      if (!selectionPID(track1))
        continue; // for primary particle PID

      if (!selectionTrack(track1)) {
        continue;
      }

      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hDcaxy_after"), track1.dcaXY());
        histos.fill(HIST("hDcaz_after"), track1.dcaZ());
        // histos.fill(HIST("hNsigmaPionTPC_after"), track1.tpcNSigmaPi());
        // histos.fill(HIST("hNsigmaPionTOF_after"), track1.tofNSigmaPi());
      }

      ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(),
                                        massPi);
      pions.push_back(temp1);
      PionIndex.push_back(track1.globalIndex());
      // PionIndex.push_back(track1.index());
      PioncollIndex.push_back(track1.collisionId());
      PionSign.push_back(track1.sign());
    }
    /*
    for (auto track2 : tracks) {

      if (!selectionPID(track2))
        continue; // for primary particle PID

      if (!selectionTrack(track2)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(),
                                        massKa);
      kaons.push_back(temp2);
      //PionIndex.push_back(track1.globalIndex());
      KaonIndex.push_back(track2.index());
      KaoncollIndex.push_back(track2.collisionId());
      KaonSign.push_back(track2.sign());
    } // track loop ends
    */

    for (auto& v0 : V0s) {

      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();
      double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
      double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

      if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
        continue;
      }

      if (!SelectionV0(collision, v0, centrality)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector temp2(v0.pt(), v0.eta(), v0.phi(), massK0s);
      kshorts.push_back(temp2);
      V0collIndex.push_back(v0.collisionId());
      // KshortPosDaughIndex.push_back(postrack.index());
      // KshortNegDaughIndex.push_back(negtrack.index());
      KshortPosDaughIndex.push_back(postrack.globalIndex());
      KshortNegDaughIndex.push_back(negtrack.globalIndex());
    }

    if (pions.size() != 0 && kshorts.size() != 0) {
      // if (pions.size() != 0 && kaons.size() != 0) {
      //  if (pions.size() != 0 && pions2.size() != 0) {
      for (auto ipion = pions.begin(); ipion != pions.end(); ++ipion) {
        auto i1 = std::distance(pions.begin(), ipion);
        if (PionSign.at(i1) == 0)
          continue;
        for (auto ikshort = kshorts.begin(); ikshort != kshorts.end();
             ++ikshort) {
          // for (auto ikaon = kaons.begin(); ikaon != kaons.end();
          //  ++ikaon) {
          //  for (auto ikshort = pions2.begin(); ikshort != pions2.end();
          //     ++ikshort) {
          auto i3 = std::distance(kshorts.begin(), ikshort);
          // auto i3 = std::distance(kaons.begin(), ikaon);
          if (PionIndex.at(i1) == KshortPosDaughIndex.at(i3))
            continue;
          if (PionIndex.at(i1) == KshortNegDaughIndex.at(i3))
            continue;
          // if (KaonIndex.at(i3) <= PionIndex.at(i1))
          // continue;
          if (PioncollIndex.at(i1) != V0collIndex.at(i3))
            continue;

          // if (PionSign.at(i1) * KaonSign.at(i3) >= 0)
          // continue;

          CKSVector = pions.at(i1) + kshorts.at(i3);
          // CKSVector = pions.at(i1) + kaons.at(i3);
          if (TMath::Abs(CKSVector.Rapidity()) < 0.5) {
            histos.fill(HIST("h3CKSInvMassUnlikeSign"), centrality,
                        CKSVector.Pt(), CKSVector.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(chargedkstaranalysis, processSE, "Process Same event", false);

  void processME(EventCandidates const& /*collisions*/,
                 TrackCandidates const& /*tracks*/, V0TrackCandidate const& /*V0s*/)

  /*
    void processME(EventCandidates const& collisions,
    TrackCandidates const& tracks)*/
  {

    // histos.fill(HIST("counter"), 1.5);
    /*
      auto tracksTuple = std::make_tuple(tracks);
      //////// currently mixing the event with similar TPC multiplicity ////////
      BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
      SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, t
      racksTuple, &cache};
    */

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      // histos.fill(HIST("counter"), 2.5);

      auto centrality = c1.centFT0C();
      auto centrality2 = c2.centFT0C();

      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }

      if (TVXEvsel && (!c1.selection_bit(aod::evsel::kIsTriggerTVX) || !c2.selection_bit(aod::evsel::kIsTriggerTVX))) {
        continue;
      }

      if (additionalEvsel && !eventSelected(c1, centrality)) {
        continue;
      }
      if (additionalEvsel && !eventSelected(c2, centrality2)) {
        continue;
      }

      for (auto& [t1, t2] : o2::soa::combinations(
             o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        // histos.fill(HIST("counter"), 3.5);
        if (!selectionTrack(t1))
          continue;
        // histos.fill(HIST("counter"), 4.5);
        if (!selectionPID(t1))
          continue;
        // histos.fill(HIST("counter"), 5.5);
        if (t1.sign() == 0)
          continue;
        // histos.fill(HIST("counter"), 6.5);

        /*if (!selectionTrack(t2))
                continue;
          if (!selectionPID(t2))
                continue;
        */

        if (!SelectionV0(c2, t2, centrality2))
          continue;

        // histos.fill(HIST("counter"), 7.5);

        auto postrack = t2.template posTrack_as<TrackCandidates>();
        auto negtrack = t2.template negTrack_as<TrackCandidates>();
        double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
        double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
          continue;
        }
        // histos.fill(HIST("counter"), 8.5);
        if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
          continue;
        }
        // histos.fill(HIST("counter"), 9.5);

        // if (t1.sign() * t2.sign() >= 0)
        // continue;

        TLorentzVector pi;
        pi.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massPi);
        TLorentzVector KSh;
        KSh.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massK0s);

        TLorentzVector CKSmix = pi + KSh;

        if (TMath::Abs(CKSmix.Rapidity()) < 0.5) {
          histos.fill(HIST("h3CKSInvMassMixed"), centrality, CKSmix.Pt(),
                      CKSmix.M());
        }
      }
    }
  }

  PROCESS_SWITCH(chargedkstaranalysis, processME, "Process Mixed event", false);

  void processGenMC(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {

    if (std::abs(mcCollision.posZ()) < cutzvertex)
      rGenParticles.fill(HIST("hMC"), 0.5);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    auto cent = 0;
    for (const auto& collision : collisions) {
      // if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 1.5);
      /*if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
  }*/
      if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 2.5);

      cent = collision.centFT0C();
      rGenParticles.fill(HIST("hCentGen"), cent);
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();

    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }

    rGenParticles.fill(HIST("hMC"), 3.5);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 4.5);
      if (std::abs(mcParticle.pdgCode()) != 323) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 5.5);
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }

      rGenParticles.fill(HIST("hMC"), 6.5);
      auto daughts = false;
      auto daughtp = false;
      // int count = 0;
      for (auto kCurrentDaughter : kDaughters) {
        // LOG(info) << "Daughters PDG:\t" << count<<" "<<kCurrentDaughter.pdgCode();
        if (kCurrentDaughter.pdgCode() == 311) {
          auto kDaughter2 = kCurrentDaughter.daughters_as<aod::McParticles>();
          for (auto kCurrentDaughter2 : kDaughter2) {
            if (kCurrentDaughter2.pdgCode() == 310)
              daughts = true;
          }
        } else if (std::abs(kCurrentDaughter.pdgCode()) == 211) {
          if (kCurrentDaughter.isPhysicalPrimary() == 1)
            daughtp = true;
        }
        // count += 1;
      }
      rGenParticles.fill(HIST("hMC"), 7.5);
      if (daughtp && daughts) {
        rGenParticles.fill(HIST("hCKSGen"), mcParticle.pt(), cent);
      }
    }
  }

  void processRecMC(EventCandidatesMC::iterator const& collision,
                    TrackCandidatesMC const& tracks, V0TrackCandidatesMC const& V0s,
                    aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)

  {

    if (!collision.has_mcCollision()) {
      return;
    }
    // if (std::abs(collision.mcCollision().posZ()) > cutzvertex || !collision.sel8()) {
    if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      return;
    }

    rRecParticles.fill(HIST("hMCRec"), 0.5);

    // if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
      return;
    }
    if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }

    auto cent = 0;
    cent = collision.centFT0C();

    rRecParticles.fill(HIST("hMCRec"), 1.5);
    rRecParticles.fill(HIST("hCentRec"), cent);

    float centrality = 0.0f;
    auto oldindex = -999;

    for (auto track1 : tracks) {

      if (!selectionPID(track1))
        continue; // for primary particle PID

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (!selectionTrack(track1)) {
        continue;
      }

      auto mctrack1 = track1.mcParticle();

      if (!mctrack1.isPhysicalPrimary()) {
        continue;
      }

      for (auto& v0 : V0s) {

        if (!v0.has_mcParticle()) {
          continue;
        }

        auto postrack = v0.template posTrack_as<TrackCandidatesMC>();
        auto negtrack = v0.template negTrack_as<TrackCandidatesMC>();

        if (!postrack.has_mcParticle())
          continue;                     // Checking that the daughter tracks come from particles and are not fake
        if (!negtrack.has_mcParticle()) // Checking that the daughter tracks come from particles and are not fake
          continue;

        // auto posParticle = postrack.mcParticle();
        // auto negParticle = negtrack.mcParticle();

        double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
        double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
          continue;
        }

        if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
          continue;
        }

        if (!SelectionV0(collision, v0, centrality)) {
          continue;
        }

        auto mctrackv0 = v0.mcParticle();

        int track1PDG = std::abs(mctrack1.pdgCode());
        // int track2PDG = std::abs(mctrack2.pdgCode());
        int trackv0PDG = std::abs(mctrackv0.pdgCode());

        if (postrack.globalIndex() == track1.globalIndex())
          continue;
        if (negtrack.globalIndex() == track1.globalIndex())
          continue;

        rRecParticles.fill(HIST("hMCRec"), 2.5);

        if (track1PDG != 211) {
          continue;
        }
        // if (track2PDG != 321) {
        // continue;
        // }
        if (trackv0PDG != 310) {
          continue;
        }

        rRecParticles.fill(HIST("hMCRec"), 3.5);

        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          // for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrackv0.mothers_as<aod::McParticles>()) {

            rRecParticles.fill(HIST("hMCRec"), 4.5);
            // LOG(info) << "Initial Mothers PDG:\t" <<mothertrack1.pdgCode()<<" "<<mothertrack2.pdgCode();

            if (mothertrack2.pdgCode() != 311) // K0
              continue;

            for (auto& mothertrack3 : mothertrack2.mothers_as<aod::McParticles>()) {

              // LOG(info) << "final Mothers PDG:\t" <<mothertrack3.pdgCode();

              if (std::abs(mothertrack3.pdgCode()) != 323) {
                continue;
              }

              if (mothertrack3.pdgCode() != mothertrack1.pdgCode()) {
                continue;
              }

              // LOG(info) << "final Mothers PDG:\t" <<mothertrack3.pdgCode()<<" "<<mothertrack3.globalIndex()<<" "<<mothertrack1.globalIndex();

              rRecParticles.fill(HIST("hMCRec"), 5.5);

              if (mothertrack3.globalIndex() != mothertrack1.globalIndex()) {
                continue;
              }

              rRecParticles.fill(HIST("hMCRec"), 6.5);
              // rRecParticles.fill(HIST("hMCRec"), 7.5);

              // LOG(info) << "Mothers PDG:\t" <<mothertrack1.pdgCode()<<" "<<mothertrack2.pdgCode();

              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                continue;
              }

              oldindex = mothertrack1.globalIndex();

              if (std::abs(mothertrack1.y()) >= 0.5) {
                continue;
              }

              rRecParticles.fill(HIST("hMCRec"), 7.5);

              rRecParticles.fill(HIST("hCKSRec"), mothertrack1.pt(), cent);
            }
          }
        }
      }
    } // track loop ends
  }

  PROCESS_SWITCH(chargedkstaranalysis, processGenMC, "Process Gen event", false);
  PROCESS_SWITCH(chargedkstaranalysis, processRecMC, "Process Rec event", false);

  void processSEnew(aod::ResoCollision& collision, aod::ResoTracks const& resotracks, aod::ResoV0s const& resov0s)
  {
    // Fill the event counter
    histos1.fill(HIST("hVertexZ"), collision.posZ());
    fillHistograms<false, false>(collision, resotracks, resov0s); // Fill histograms, no MC, no mixing
  }
  PROCESS_SWITCH(chargedkstaranalysis, processSEnew, "Process Same event new", true);

  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMEnew(aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks, aod::ResoV0s const& resov0s)
  {
    auto tracksV0sTuple = std::make_tuple(resotracks, resov0s);
    auto V0sTuple = std::make_tuple(resov0s);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    Pair<aod::ResoCollisions, aod::ResoTracks, aod::ResoV0s, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksV0sTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, restrk1, c2, resov0s2] : pairs) {
      fillHistograms<false, true>(c1, restrk1, resov0s2);
    }
  }
  PROCESS_SWITCH(chargedkstaranalysis, processMEnew, "Process Mixed events new", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<chargedkstaranalysis>(cfgc)};
}
