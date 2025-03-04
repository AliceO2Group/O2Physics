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

/// \file lambda1520analysis.cxx
/// \brief This task reconstructs track-track decay lambda(1520) resonance candidate
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#include "TLorentzVector.h"
#include "TF1.h"
#include "TRandom3.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::constants::physics;

struct Lambda1520analysis {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<aod::ResoTracks> perResoCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  // Configurables
  // switches
  Configurable<bool> cEtaAssym{"cEtaAssym", false, "Turn on/off EtaAssym calculation"};
  Configurable<bool> isFilladditionalQA{"isFilladditionalQA", false, "Turn on/off additional QA plots"};
  Configurable<bool> cOldPIDcut{"cOldPIDcut", false, "Switch to turn on/off old PID cut to apply pt dependent cut"};
  Configurable<bool> fixedPIDcut{"fixedPIDcut", false, "Switch to turn on/off FIXED PID cut to apply pt dependent cut"};
  Configurable<bool> crejectPion{"crejectPion", false, "Switch to turn on/off pion contamination"};
  Configurable<bool> cDCAr7SigCut{"cDCAr7SigCut", false, "Track DCAr 7 Sigma cut to PV Maximum"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cTPCNClsFound{"cTPCNClsFound", false, "Switch to turn on/off TPCNClsFound cut"};
  Configurable<bool> additionalQAeventPlots{"additionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> additionalMEPlots{"additionalMEPlots", false, "Additional Mixed event plots"};

  // Pre-selection Track cuts
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
  Configurable<float> cMinTPCNClsFound{"cMinTPCNClsFound", 120, "minimum TPCNClsFound value for good track"};
  Configurable<int> cMinTPCncr{"cMinTPCncr", 70, "Minimum number of TPC X rows"};

  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1f, "Track DCAz cut to PV Maximum"};

  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};

  /// PID Selections
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"}; // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"}; // Reject missidentified particles when tpc bands merge
  Configurable<bool> cUseRejNsigma{"cUseRejNsigma", false, "Switch on/off track rejection method to improve purity"};
  Configurable<bool> tofAtHighPt{"tofAtHighPt", false, "Use TOF at high pT"};
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                    // By pass TOF PID selection
  Configurable<int> pidCutType{"pidCutType", 2, "pidCutType = 1 for square cut, 2 for circular cut"}; // By pass TOF PID selection

  // Kaon
  // Old PID use case
  Configurable<std::vector<double>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {3}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {3}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTPCTOFCombinedpTintv{"kaonTPCTOFCombinedpTintv", {999.}, "pT intervals for Kaon TPC-TOF PID cuts"};
  Configurable<std::vector<double>> kaonTPCTOFCombinedPIDcuts{"kaonTPCTOFCombinedPIDcuts", {3}, "nSigma list for Kaon TPC-TOF PID cuts"};
  Configurable<double> cMaxTPCnSigmaKaonVETO{"cMaxTPCnSigmaKaonVETO", 3.0, "TPC nSigma VETO cut for Kaon"}; // TPC

  // New PID use case
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};                // TPC
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"};   // Combined
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection

  // Proton
  // Old PID use case
  Configurable<std::vector<double>> protonTPCPIDpTintv{"protonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTPCPIDcuts{"protonTPCPIDcuts", {3}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDpTintv{"protonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDcuts{"protonTOFPIDcuts", {3}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> protonTPCTOFCombinedpTintv{"protonTPCTOFCombinedpTintv", {999.}, "pT intervals for Proton TPC-TOF PID cuts"};
  Configurable<std::vector<double>> protonTPCTOFCombinedPIDcuts{"protonTPCTOFCombinedPIDcuts", {3}, "nSigma list for Proton TPC-TOF PID cuts"};
  Configurable<double> cMaxTPCnSigmaProtonVETO{"cMaxTPCnSigmaProtonVETO", 3.0, "TPC nSigma VETO cut for Proton"}; // TPC

  // New PID use case
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"};     // Use only TOF track for Pion PID selection

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};

  // MC Event selection
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  // cuts on mother
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", false, "Enamble additional cuts on mother"};
  Configurable<double> cMaxPtMotherCut{"cMaxPtMotherCut", 10.0, "Maximum pt of mother cut"};
  Configurable<double> cMaxMinvMotherCut{"cMaxMinvMotherCut", 3.0, "Maximum Minv of mother cut"};
  Configurable<bool> cfgCutsOnDaughters{"cfgCutsOnDaughters", false, "Enamble additional cuts on daughters"};
  Configurable<int> cetaphiBins{"cetaphiBins", 400, "number of eta and phi bins"};
  Configurable<double> cMaxDeltaEtaCut{"cMaxDeltaEtaCut", 0.7, "Maximum deltaEta between daughters"};
  Configurable<double> cMaxDeltaPhiCut{"cMaxDeltaPhiCut", 1.5, "Maximum deltaPhi between daughters"};
  Configurable<bool> invmass1D{"invmass1D", false, "Invariant mass 1D"};
  Configurable<bool> cAdditionalMCPlots{"cAdditionalMCPlots", false, "Draw additional plots related to MC"};

  TRandom* rn = new TRandom();

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 4.6, 4.8, 4.9, 5.0, 5.5, 5.6, 6.0, 6.4, 6.5, 7.0, 7.2, 8.0, 9.0, 9.5, 9.6, 10.0, 11.0, 11.5, 12.0, 13.0, 14.0, 14.4, 15.0, 16.0, 18.0, 19.2, 20.}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, ""};
  ConfigurableAxis binsMass{"binsMass", {500, 1.3, 2.3}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {110, 0.0, 110.0}, "mult_{FT0M}"};
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {100, 60, 160}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5, 6.5}, ""};
  ConfigurableAxis binsnTPCSignal{"binsnTPCSignal", {1000, 0, 1000}, ""};
  ConfigurableAxis occupancybins{"occupancybins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning of the occupancy axis"};
  Configurable<bool> applyOccupancyCut{"applyOccupancyCut", false, "Apply occupancy cut"};
  Configurable<int> occupancyCut{"occupancyCut", 1000, "Mimimum Occupancy cut"};

  // Rotational background
  Configurable<bool> isCalcRotBkg{"isCalcRotBkg", true, "Calculate rotational background"};
  Configurable<int> rotationalcut{"rotationalcut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<int> cNofRotations{"cNofRotations", 3, "Number of random rotations in the rotational background"};

  void init(o2::framework::InitContext&)
  {
    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, "#eta"};
    AxisSpec axisMassLambda1520{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{V0M}"};
    AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    AxisSpec axisTPCXrow{binsTPCXrows, "#Xrows_{TPC}"};
    AxisSpec pidQAAxis = {binsnSigma, "#sigma"};
    AxisSpec axisTPCSignal = {binsnTPCSignal, ""};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec occupancyaxis = {occupancybins, "Occupancy [-40,100]"};

    if (additionalQAeventPlots) {
      // Test on Mixed event
      histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
      histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});

      // event histograms
      histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

      histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
    }

    if (doprocessData) {
      // Track QA before cuts
      //  --- Track
      histos.add("QA/QAbefore/Track/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPConly_Nsigma_ka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPConly_Nsigma_pr", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAbefore/Track/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAbefore/Track/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAbefore/Track/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAbefore/Track/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});

      if (isFilladditionalQA) {
        // TPC ncluster distirbutions
        histos.add("TPCncluster/TPCnclusterpr", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
        histos.add("TPCncluster/TPCnclusterka", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
        histos.add("TPCncluster/TPCnclusterPhipr", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
        histos.add("TPCncluster/TPCnclusterPhika", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

        // Multiplicity correlation calibrations
        histos.add("MultCalib/centglopr", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/centgloka", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/GloPVpr", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
        histos.add("MultCalib/GloPVka", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
      }

      // PID QA after cuts
      //  --- Kaon
      histos.add("QA/QAafter/Kaon/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Kaon/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Kaon/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Kaon/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Kaon/TPC_Signal_ka_all", "TPC Signal for Kaon;#it{p}_{T} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});

      //  --- Proton
      histos.add("QA/QAafter/Proton/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/dcaZ", "DCA_{Z} distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Proton/dcaXY", "DCA_{XY} momentum distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Proton/TPC_CR", "# TPC Xrows distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Proton/eta", "#eta distribution of Protons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Proton/TPC_Signal_pr_all", "TPC Signal for Proton;#it{p}_{T} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});

      //  Mass QA 1D for quick check
      if (invmass1D) {
        histos.add("Result/Data/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});   // K+ + Pr
        histos.add("Result/Data/lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr
      }
      // eta phi QA
      if (cfgCutsOnDaughters) {
        histos.add("QAbefore/deltaEta", "deltaEta of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
        histos.add("QAbefore/deltaPhi", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});

        histos.add("QAafter/deltaEta", "deltaEta of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
        histos.add("QAafter/deltaPhi", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});

        histos.add("QAafter/deltaEtaafter", "deltaEta of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
        histos.add("QAafter/deltaPhiafter", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
        histos.add("QAafter/EtaPrafter", "Eta of  proton candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
        histos.add("QAafter/PhiPrafter", "Phi of  proton candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});
        histos.add("QAafter/EtaKaafter", "Eta of kaon  candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
        histos.add("QAafter/PhiKaafter", "Phi of kaon  candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});
      }

      if (isCalcRotBkg) {
        histos.add("Result/Data/h3lambda1520InvMassRotation", "Invariant mass of #Lambda(1520) rotation", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520, occupancyaxis});
      }

      // 3d histogram
      histos.add("Result/Data/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
      histos.add("Result/Data/h3lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
    }
    if (doprocessME) {
      if (invmass1D) {
        histos.add("Result/Data/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      }
      histos.add("Result/Data/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});

      if (additionalMEPlots) {
        histos.add("Result/Data/lambda1520invmassME_DS", "Invariant mass of #Lambda(1520) mixed event DS", kTH1F, {axisMassLambda1520});
        histos.add("Result/Data/lambda1520invmassME_DSAnti", "Invariant mass of #Lambda(1520) mixed event DSAnti", kTH1F, {axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassME_DS", "Invariant mass of #Lambda(1520) mixed event DS", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassME_DSAnti", "Invariant mass of #Lambda(1520) mixed event DSAnti", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      }
    }

    if (cEtaAssym) {
      histos.add("Result/Data/hlambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", {HistType::kTH1F, {axisMassLambda1520}});

      histos.add("Result/Data/h3lambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      if (doprocessME) {
        histos.add("Result/Data/hlambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/hlambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", {HistType::kTH1F, {axisMassLambda1520}});

        histos.add("Result/Data/h3lambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      }
    }
    //}
    // MC QA
    if (doprocessMCTrue) {
      histos.add("Result/MC/Genlambda1520pt", "pT distribution of True MC #Lambda(1520)0", kTH3F, {mcLabelAxis, axisPt, axisMult});
      histos.add("Result/MC/Genantilambda1520pt", "pT distribution of True MC Anti-#Lambda(1520)0", kTH3F, {mcLabelAxis, axisPt, axisMult});
    }
    if (doprocessMC) {
      histos.add("QA/MC/trkDCAxy_pr", "DCAxy distribution of proton track candidates", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAz_pr", "DCAz distribution of proton track candidates", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/MC/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/MC/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});

      histos.add("Result/MC/h3lambda1520Recoinvmass", "Invariant mass of Reconstructed MC #Lambda(1520)0", kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/MC/h3antilambda1520Recoinvmass", "Invariant mass of Reconstructed MC Anti-#Lambda(1520)0", kTH3F, {axisMult, axisPt, axisMassLambda1520});
      if (cAdditionalMCPlots) {
        histos.add("Result/MC/lambda1520Reco", "pT distribution of Reconstructed MC #Lambda(1520)0", kTH2F, {axisPt, axisMult});
        histos.add("Result/MC/antilambda1520Reco", "pT distribution of Reconstructed MC Anti-#Lambda(1520)0", kTH2F, {axisPt, axisMult});
        histos.add("Result/MC/hlambda1520Recoinvmass", "Inv mass distribution of Reconstructed MC #Lambda(1520)", kTH1F, {axisMassLambda1520});
        histos.add("Result/MC/hantilambda1520Recoinvmass", "Inv mass distribution of Reconstructed MC Anti-#Lambda(1520)", kTH1F, {axisMassLambda1520});
      }
    }

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massKa = MassKaonCharged;
  double massPr = MassProton;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (cDCAr7SigCut) {
      if (std::abs(track.dcaXY()) > (0.004f + 0.0130f / (track.pt()))) // 7 - Sigma cut
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
        return false;
    }
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cTPCNClsFound && (track.tpcNClsFound() < cMinTPCNClsFound))
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (track.tpcNClsCrossedRows() < cMinTPCncr)
      return false;

    return true;
  }

  // PID selection new PID method
  template <typename T>
  bool selectionnewPIDProton(const T& candidate)
  {
    if (tofAtHighPt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaPr()) < cMaxTOFnSigmaProton)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaPr()) < cMaxTPCnSigmaKaon)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaPr()) < cMaxTPCnSigmaProton) {
        if (cUseRejNsigma) {
          if (std::abs(candidate.tpcNSigmaPi()) > cRejNsigmaTpc && std::abs(candidate.tpcNSigmaKa()) > cRejNsigmaTpc) {
            tpcPIDPassed = true;
          }
        } else if (!cUseRejNsigma) {
          tpcPIDPassed = true;
        }
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if ((nsigmaCutCombinedProton > 0) && ((candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (nsigmaCutCombinedProton * nsigmaCutCombinedProton))) {
          if (cUseRejNsigma) {
            if (std::abs(candidate.tofNSigmaPi()) > cRejNsigmaTof && std::abs(candidate.tofNSigmaKa()) > cRejNsigmaTof) {
              tofPIDPassed = true;
            }
          } else if (!cUseRejNsigma) {
            tofPIDPassed = true;
          }
        } else if ((nsigmaCutCombinedProton <= 0) && (std::abs(candidate.tofNSigmaPr()) < cMaxTOFnSigmaProton)) {
          if (cUseRejNsigma) {
            if (std::abs(candidate.tofNSigmaPi()) > cRejNsigmaTof && std::abs(candidate.tofNSigmaKa()) > cRejNsigmaTof) {
              tofPIDPassed = true;
            }
          } else if (!cUseRejNsigma) {
            tofPIDPassed = true;
          }
        }
      } else {
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionnewPIDKaon(const T& candidate)
  {
    if (tofAtHighPt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
        if (cUseRejNsigma) {
          if (std::abs(candidate.tpcNSigmaPi()) > cRejNsigmaTpc && std::abs(candidate.tpcNSigmaPr()) > cRejNsigmaTpc) {
            tpcPIDPassed = true;
          }
        } else if (!cUseRejNsigma) {
          tpcPIDPassed = true;
        }
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if ((nsigmaCutCombinedKaon > 0) && ((candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < (nsigmaCutCombinedKaon * nsigmaCutCombinedKaon))) {
          if (cUseRejNsigma) {
            if (std::abs(candidate.tofNSigmaPi()) > cRejNsigmaTof && std::abs(candidate.tofNSigmaPr()) > cRejNsigmaTof) {
              tofPIDPassed = true;
            }
          } else if (!cUseRejNsigma) {
            tofPIDPassed = true;
          }
        } else if ((nsigmaCutCombinedKaon <= 0) && (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon)) {
          if (cUseRejNsigma) {
            if (std::abs(candidate.tofNSigmaPi()) > cRejNsigmaTof && std::abs(candidate.tofNSigmaPr()) > cRejNsigmaTof) {
              tofPIDPassed = true;
            }
          } else if (!cUseRejNsigma) {
            tofPIDPassed = true;
          }
        }
      } else {
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }

  // PID selection old PID method
  template <typename T>
  bool selectionoldPIDProton(const T& candidate)
  {
    auto vProtonTPCPIDpTintv = static_cast<std::vector<double>>(protonTPCPIDpTintv);
    vProtonTPCPIDpTintv.insert(vProtonTPCPIDpTintv.begin(), cMinPtcut);
    auto vProtonTPCPIDcuts = static_cast<std::vector<double>>(protonTPCPIDcuts);
    auto vProtonTOFPIDpTintv = static_cast<std::vector<double>>(protonTOFPIDpTintv);
    auto vProtonTPCTOFCombinedpTintv = static_cast<std::vector<double>>(protonTPCTOFCombinedpTintv);
    auto vProtonTPCTOFCombinedPIDcuts = static_cast<std::vector<double>>(protonTPCTOFCombinedPIDcuts);
    auto vProtonTOFPIDcuts = static_cast<std::vector<double>>(protonTOFPIDcuts);
    auto lengthOfprotonTPCPIDpTintv = static_cast<int>(vProtonTPCPIDpTintv.size());
    auto lengthOfprotonTOFPIDpTintv = static_cast<int>(vProtonTOFPIDpTintv.size());
    auto lengthOfprotonTPCTOFCombinedPIDpTintv = static_cast<int>(vProtonTPCTOFCombinedpTintv.size());

    bool isTrk1Selected{true};

    // For Proton candidate:
    if (candidate.hasTOF()) {
      if (pidCutType == 1) {
        if (lengthOfprotonTOFPIDpTintv > 0) {
          if (candidate.pt() > vProtonTOFPIDpTintv[lengthOfprotonTOFPIDpTintv - 1]) {
            isTrk1Selected = false;
          } else {
            for (int i = 0; i < lengthOfprotonTOFPIDpTintv; i++) {
              if (candidate.pt() < vProtonTOFPIDpTintv[i]) {

                if (std::abs(candidate.tofNSigmaPr()) > vProtonTOFPIDcuts[i])
                  isTrk1Selected = false;
                if (std::abs(candidate.tpcNSigmaPr()) > cMaxTPCnSigmaProtonVETO)
                  isTrk1Selected = false;
              }
            }
          }
        }
      } else if (pidCutType == 2) {
        if (lengthOfprotonTPCTOFCombinedPIDpTintv > 0) {
          if (candidate.pt() > vProtonTPCTOFCombinedpTintv[lengthOfprotonTPCTOFCombinedPIDpTintv - 1]) {
            isTrk1Selected = false;
          } else {
            for (int i = 0; i < lengthOfprotonTPCTOFCombinedPIDpTintv; i++) {
              if (candidate.pt() < vProtonTPCTOFCombinedpTintv[i]) {
                if ((candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr()) > (vProtonTPCTOFCombinedPIDcuts[i] * vProtonTPCTOFCombinedPIDcuts[i]))
                  isTrk1Selected = false;
              }
            }
          }
        }
      }
    } else {
      if (lengthOfprotonTPCPIDpTintv > 0) {
        if (candidate.pt() > vProtonTPCPIDpTintv[lengthOfprotonTPCPIDpTintv - 1]) {
          isTrk1Selected = false;
        } else {
          for (int i = 0; i < lengthOfprotonTPCPIDpTintv; i++) {
            if (candidate.pt() > vProtonTPCPIDpTintv[i] && candidate.pt() < vProtonTPCPIDpTintv[i + 1]) {
              if (std::abs(candidate.tpcNSigmaPr()) > vProtonTPCPIDcuts[i])
                isTrk1Selected = false;
            }
          }
        }
      }
    }
    return isTrk1Selected;
  }

  template <typename T>
  bool selectionoldPIDKaon(const T& candidate)
  {
    auto vKaonTPCPIDpTintv = static_cast<std::vector<double>>(kaonTPCPIDpTintv);
    vKaonTPCPIDpTintv.insert(vKaonTPCPIDpTintv.begin(), cMinPtcut);
    auto vKaonTPCPIDcuts = static_cast<std::vector<double>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<double>>(kaonTOFPIDpTintv);
    auto vKaonTPCTOFCombinedpTintv = static_cast<std::vector<double>>(kaonTPCTOFCombinedpTintv);
    auto vKaonTPCTOFCombinedPIDcuts = static_cast<std::vector<double>>(kaonTPCTOFCombinedPIDcuts);
    auto vKaonTOFPIDcuts = static_cast<std::vector<double>>(kaonTOFPIDcuts);
    auto lengthOfkaonTPCPIDpTintv = static_cast<int>(vKaonTPCPIDpTintv.size());
    auto lengthOfkaonTOFPIDpTintv = static_cast<int>(vKaonTOFPIDpTintv.size());
    auto lengthOfkaonTPCTOFCombinedPIDpTintv = static_cast<int>(vKaonTPCTOFCombinedpTintv.size());

    bool isTrk2Selected{true};

    // For Kaon candidate:
    if (candidate.hasTOF()) {
      if (pidCutType == 1) {
        if (lengthOfkaonTOFPIDpTintv > 0) {
          if (candidate.pt() > vKaonTOFPIDpTintv[lengthOfkaonTOFPIDpTintv - 1]) {
            isTrk2Selected = false;
          } else {
            for (int i = 0; i < lengthOfkaonTOFPIDpTintv; i++) {
              if (candidate.pt() < vKaonTOFPIDpTintv[i]) {

                if (std::abs(candidate.tofNSigmaKa()) > vKaonTOFPIDcuts[i])
                  isTrk2Selected = false;
                if (std::abs(candidate.tpcNSigmaKa()) > cMaxTPCnSigmaKaonVETO)
                  isTrk2Selected = false;
              }
            }
          }
        }
      } else if (pidCutType == 2) {
        if (lengthOfkaonTPCTOFCombinedPIDpTintv > 0) {
          if (candidate.pt() > vKaonTPCTOFCombinedpTintv[lengthOfkaonTPCTOFCombinedPIDpTintv - 1]) {
            isTrk2Selected = false;
          } else {
            for (int i = 0; i < lengthOfkaonTPCTOFCombinedPIDpTintv; i++) {
              if (candidate.pt() < vKaonTPCTOFCombinedpTintv[i]) {
                if ((candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) > (vKaonTPCTOFCombinedPIDcuts[i] * vKaonTPCTOFCombinedPIDcuts[i]))
                  isTrk2Selected = false;
              }
            }
          }
        }
      }
    } else {
      if (lengthOfkaonTPCPIDpTintv > 0) {
        if (candidate.pt() > vKaonTPCPIDpTintv[lengthOfkaonTPCPIDpTintv - 1]) {
          isTrk2Selected = false;
        } else {
          for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
            if (candidate.pt() > vKaonTPCPIDpTintv[i] && candidate.pt() < vKaonTPCPIDpTintv[i + 1]) {
              if (std::abs(candidate.tpcNSigmaKa()) > vKaonTPCPIDcuts[i])
                isTrk2Selected = false;
            }
          }
        }
      }
    }
    return isTrk2Selected;
  }

  // newly added to test
  template <typename T>
  bool selectionPIDProtonFixed(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (candidate.pt() < 1.5 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0 && std::abs(candidate.tofNSigmaPr()) < 4.0) {
        return true;
      }
      if (candidate.pt() >= 1.5 && candidate.pt() < 2.0 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -3.0 && candidate.tofNSigmaPr() < 4.0) {
        return true;
      }
      if (candidate.pt() >= 2.0 && candidate.pt() < 2.5 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < 4.0) {
        return true;
      }
      if (candidate.pt() >= 2.5 && candidate.pt() < 3.0 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < 3.0) {
        return true;
      }
      if (candidate.pt() >= 3.0 && candidate.pt() < 4.0 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < 2.0) {
        return true;
      }
    }
    if (!candidate.hasTOF()) {
      if (candidate.pt() < 0.4 && std::abs(candidate.tpcNSigmaPr()) < 4.0) {
        return true;
      }
      if (candidate.pt() >= 0.4 && candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (candidate.pt() >= 0.5 && candidate.pt() < 0.7 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < 2.5) {
        return true;
      }
      if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.5 && candidate.tpcNSigmaPr() < 2.5) {
        return true;
      }
      if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < 2.5) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDKaonFixed(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (candidate.pt() < 0.8 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0 && std::abs(candidate.tofNSigmaKa()) < 4.0) {
        return true;
      }
      if (candidate.pt() >= 0.8 && candidate.pt() < 1.3 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0 && candidate.tofNSigmaKa() > -3.0 && candidate.tofNSigmaKa() < 4.0) {
        return true;
      }
      if (candidate.pt() >= 1.3 && candidate.pt() < 1.6 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0 && candidate.tofNSigmaKa() > -2.0 && candidate.tofNSigmaKa() < 3.0) {
        return true;
      }
      if (candidate.pt() >= 1.6 && candidate.pt() < 1.8 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0 && candidate.tofNSigmaKa() > -1.5 && candidate.tofNSigmaKa() < 2.5) {
        return true;
      }
      if (candidate.pt() >= 1.8 && candidate.pt() < 2.5 && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0 && candidate.tofNSigmaKa() > -1.0 && candidate.tofNSigmaKa() < 2.0) {
        return true;
      }
    }
    if (!candidate.hasTOF()) {
      if (candidate.pt() < 0.3 && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
        return true;
      }
      if (candidate.pt() >= 0.3 && candidate.pt() < 0.4 && candidate.tpcNSigmaKa() > -2.0 && candidate.tpcNSigmaKa() < 2.5) {
        return true;
      }
      if (candidate.pt() >= 0.4 && candidate.pt() < 0.5 && candidate.tpcNSigmaKa() > -1.0 && candidate.tpcNSigmaKa() < 2.5) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool rejectPion(const T& candidate)
  {
    if (candidate.pt() > 1.0 && candidate.pt() < 2.0 && !candidate.hasTOF() && candidate.tpcNSigmaPi() < 2) {
      return false;
    }
    return true;
  }

  template <bool IsData, bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();

    // LOG(info) << "Before pass, Collision index:" << collision.index() << "multiplicity: " << collision.cent() << std::endl;

    auto occupancyNo = collision.trackOccupancyInTimeRange();
    if (applyOccupancyCut && occupancyNo < occupancyCut) {
      return;
    }

    // Multiplicity correlation calibration plots
    if (isFilladditionalQA) {
      if constexpr (IsData) {
        histos.fill(HIST("MultCalib/centglopr"), multiplicity, dTracks1.size());
        histos.fill(HIST("MultCalib/centgloka"), multiplicity, dTracks2.size());
        histos.fill(HIST("MultCalib/GloPVpr"), dTracks1.size(), collision.multNTracksPV());
        histos.fill(HIST("MultCalib/GloPVka"), dTracks2.size(), collision.multNTracksPV());
      }
    }

    if (additionalQAeventPlots) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }
    // LOG(info) << "After pass, Collision index:" << collision.index() << "multiplicity: " << collision.cent() << std::endl;
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance, ldaughterRot, lresonanceRot;

    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (additionalQAeventPlots) {
        if constexpr (IsData) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else if (IsMix) {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
        }
      }

      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      //// Initialize variables
      // Trk1: Proton, Trk2: Kaon
      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();

      auto trk1ptPr = trk1.pt();
      auto trk1NSigmaPrTPC = trk1.tpcNSigmaPr();
      auto trk1NSigmaPrTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPr() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      auto deltaEta = std::abs(trk1.eta() - trk2.eta());
      auto deltaPhi = std::abs(trk1.phi() - trk2.phi());
      deltaPhi = (deltaPhi > o2::constants::math::PI) ? (o2::constants::math::TwoPI - deltaPhi) : deltaPhi;

      //// QA plots before the selection
      //  --- Track QA all
      if constexpr (IsData) {
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_pr"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTPC);
        }

        histos.fill(HIST("QA/QAbefore/Track/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAbefore/Track/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAbefore/Track/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAbefore/Track/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAbefore/Track/eta"), trk1.eta());
        if (cfgCutsOnDaughters) {
          histos.fill(HIST("QAbefore/deltaEta"), deltaEta);
          histos.fill(HIST("QAbefore/deltaPhi"), deltaPhi);
        }
      }

      //// Apply the pid selection
      if (cUseOnlyTOFTrackPr && !isTrk1hasTOF)
        continue;
      if (cUseOnlyTOFTrackKa && !isTrk2hasTOF)
        continue;
      if (crejectPion && rejectPion(trk2))
        continue;
      if (cOldPIDcut) {
        if (!selectionoldPIDProton(trk1) || !selectionoldPIDKaon(trk2))
          continue;
      } else if (!cOldPIDcut && !fixedPIDcut) {
        if (!selectionnewPIDProton(trk1) || !selectionnewPIDKaon(trk2))
          continue;
      } else if (fixedPIDcut) {
        if (!selectionPIDProtonFixed(trk1) || !selectionPIDKaonFixed(trk2))
          continue;
      }

      //// QA plots after the selection
      if constexpr (IsData) { //  --- PID QA Proton
        histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        histos.fill(HIST("QA/QAafter/Proton/TPC_Signal_pr_all"), trk1ptPr, trk1.tpcSignal());
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAafter/Proton/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAafter/Proton/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAafter/Proton/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAafter/Proton/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Proton/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAafter/Proton/eta"), trk1.eta());

        //  --- PID QA Kaon
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Signal_ka_all"), trk2ptKa, trk2.tpcSignal());
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly"), trk2ptKa, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/QAafter/Kaon/dcaZ"), trk2ptKa, trk2.dcaZ());
        histos.fill(HIST("QA/QAafter/Kaon/dcaXY"), trk2ptKa, trk2.dcaXY());
        histos.fill(HIST("QA/QAafter/Kaon/TPC_CR"), trk2ptKa, trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Kaon/pT"), trk2ptKa);
        histos.fill(HIST("QA/QAafter/Kaon/eta"), trk2.eta());
        if (cfgCutsOnDaughters) {
          histos.fill(HIST("QAafter/deltaEta"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhi"), deltaPhi);
        }
        if (isFilladditionalQA) {
          // TPCncluster distributions
          histos.fill(HIST("TPCncluster/TPCnclusterpr"), trk1.tpcNClsFound());
          histos.fill(HIST("TPCncluster/TPCnclusterka"), trk2.tpcNClsFound());
          histos.fill(HIST("TPCncluster/TPCnclusterPhipr"), trk1.tpcNClsFound(), trk1.phi());
          histos.fill(HIST("TPCncluster/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());
        }
      }

      // Apply kinematic cuts.
      if (cKinCuts) {
        TVector3 v1(trk1.px(), trk1.py(), trk1.pz());
        TVector3 v2(trk2.px(), trk2.py(), trk2.pz());
        float alpha = v1.Angle(v2);
        if (alpha > 1.4 && alpha < 2.4)
          continue;
      }

      //// Resonance reconstruction
      lDecayDaughter1.SetPtEtaPhiM(trk1.pt(), trk1.eta(), trk1.phi(), massPr);
      lDecayDaughter2.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (std::abs(lResonance.Rapidity()) > 0.5)
        continue;

      if (cfgCutsOnMother) {
        if (lResonance.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonance.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      if (cfgCutsOnDaughters) {
        if (deltaEta >= cMaxDeltaEtaCut)
          continue;
        if (deltaPhi >= cMaxDeltaPhiCut)
          continue;

        if constexpr (!IsMix) {
          histos.fill(HIST("QAafter/EtaPrafter"), trk1.eta());
          histos.fill(HIST("QAafter/PhiPrafter"), trk1.phi());
          histos.fill(HIST("QAafter/EtaKaafter"), trk2.eta());
          histos.fill(HIST("QAafter/PhiKaafter"), trk2.phi());
          histos.fill(HIST("QAafter/deltaEtaafter"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhiafter"), deltaPhi);
        }
      }

      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (IsData) {
          if (isCalcRotBkg) {
            for (int i = 0; i < cNofRotations; i++) {
              float theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalcut, o2::constants::math::PI + o2::constants::math::PI / rotationalcut);
              ldaughterRot.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi() + theta2, massKa); // for rotated background
              lresonanceRot = lDecayDaughter1 + ldaughterRot;
              histos.fill(HIST("Result/Data/h3lambda1520InvMassRotation"), multiplicity, lresonanceRot.Pt(), lresonanceRot.M(), occupancyNo);
            }
          }

          if (trk1.sign() < 0) {
            if (invmass1D) {
              histos.fill(HIST("Result/Data/lambda1520invmass"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (trk1.sign() > 0) {
            if (invmass1D) {
              histos.fill(HIST("Result/Data/antilambda1520invmass"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3antilambda1520invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          }
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignAside"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignCside"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        } else if (IsMix) {
          if (invmass1D) {
            histos.fill(HIST("Result/Data/lambda1520invmassME"), lResonance.M());
          }
          histos.fill(HIST("Result/Data/h3lambda1520invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedAside"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedCside"), multiplicity, lResonance.Pt(), lResonance.M());
          }
          if (additionalMEPlots) {
            if (trk1.sign() < 0) {
              histos.fill(HIST("Result/Data/lambda1520invmassME_DS"), lResonance.M());
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (trk1.sign() > 0) {
              histos.fill(HIST("Result/Data/lambda1520invmassME_DSAnti"), lResonance.M());
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          }
        }

        // MC
        if constexpr (IsMC) {
          // LOG(info) << "trk1 pdgcode: " << trk1.pdgCode() << "trk2 pdgcode: " << trk2.pdgCode() << std::endl;

          if (std::abs(trk1.pdgCode()) != 2212 || std::abs(trk2.pdgCode()) != 321)
            continue;
          if (trk1.motherId() != trk2.motherId()) // Same mother
            continue;
          if (std::abs(trk1.motherPDG()) != 102134)
            continue;

          // Track selection check.
          histos.fill(HIST("QA/MC/trkDCAxy_pr"), trk1ptPr, trk1.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAxy_ka"), trk2ptKa, trk2.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAz_pr"), trk1ptPr, trk1.dcaZ());
          histos.fill(HIST("QA/MC/trkDCAz_ka"), trk2ptKa, trk2.dcaZ());

          histos.fill(HIST("QA/MC/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          }
          histos.fill(HIST("QA/MC/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (trk1.motherPDG() > 0) {
            if (cAdditionalMCPlots) {
              histos.fill(HIST("Result/MC/lambda1520Reco"), lResonance.Pt(), multiplicity);
              histos.fill(HIST("Result/MC/hlambda1520Recoinvmass"), lResonance.M());
            }
            histos.fill(HIST("Result/MC/h3lambda1520Recoinvmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else {
            if (cAdditionalMCPlots) {
              histos.fill(HIST("Result/MC/antilambda1520Reco"), lResonance.Pt(), multiplicity);
              histos.fill(HIST("Result/MC/hantilambda1520Recoinvmass"), lResonance.M());
            }
            histos.fill(HIST("Result/MC/h3antilambda1520Recoinvmass"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      } else {
        if constexpr (IsData) {
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignAside"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignCside"), multiplicity, lResonance.Pt(), lResonance.M());
          }
          // Like sign pair ++
          if (trk1.sign() > 0) {
            if (invmass1D) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSPP"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSPP"), multiplicity, lResonance.Pt(), lResonance.M());
          } else { // Like sign pair --
            if (invmass1D) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSMM"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSMM"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processData(aod::ResoCollision const& collision,
                   aod::ResoTracks const& resotracks)
  {
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    fillHistograms<true, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(Lambda1520analysis, processData, "Process Event for data without partition", false);

  void processMC(ResoMCCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    if (!collision.isInAfterAllCuts() || (std::abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;
    fillHistograms<false, true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(Lambda1520analysis, processMC, "Process Event for MC Light without partition", false);

  void processMCTrue(ResoMCCols::iterator const& collision, aod::ResoMCParents const& resoParents)
  {
    auto multiplicity = collision.cent();
    // Not related to the real collisions
    for (const auto& part : resoParents) {    // loop over all MC particles
      if (std::abs(part.pdgCode()) != 102134) // Lambda1520(0)
        continue;
      if (std::abs(part.y()) > 0.5) // rapidity cut
        continue;
      bool pass1 = std::abs(part.daughterPDG1()) == 321 || std::abs(part.daughterPDG2()) == 321;   // At least one decay to Kaon
      bool pass2 = std::abs(part.daughterPDG1()) == 2212 || std::abs(part.daughterPDG2()) == 2212; // At least one decay to Proton

      if (!pass1 || !pass2) // If we have both decay products
        continue;

      if (collision.isVtxIn10()) // INEL10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 0, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 0, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isInSel8()) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 1, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 1, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isTriggerTVX()) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 2, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 2, part.pt(), multiplicity);
      }
      if (collision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 3, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 3, part.pt(), multiplicity);
      }
    }
  }
  PROCESS_SWITCH(Lambda1520analysis, processMCTrue, "Process Event for MC only", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(o2::aod::ResoCollisions const& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      fillHistograms<false, false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(Lambda1520analysis, processME, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lambda1520analysis>(cfgc)};
}
