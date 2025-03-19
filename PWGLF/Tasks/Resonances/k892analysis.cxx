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
/// \file k892analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>, Sawan Sawan <sawan.sawan@cern.ch>

#include <TLorentzVector.h>
#include "TF1.h"
#include "TRandom3.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct K892analysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  ///// Configurables
  /// Histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning of the occupancy axis"};
  // ConfigurableAxis binsCent{"binsCent", {200, 0.0f, 200.0f}, "Binning of the centrality axis"};
  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};
  Configurable<int> cPIDBins{"cPIDBins", 65, "PID binning"};
  Configurable<float> cPIDQALimit{"cPIDQALimit", 6.5, "PID QA limit"};
  Configurable<int> cDCABins{"cDCABins", 150, "DCA binning"};
  Configurable<bool> invmass1D{"invmass1D", false, "Invariant mass 1D"};
  Configurable<bool> studyAntiparticle{"studyAntiparticle", false, "Study anti-particles separately"};
  Configurable<bool> fillPidPlots{"fillPidPlots", false, "Make TPC and TOF PID plots"};
  // Configurable<bool> applyOccupancyCut{"applyOccupancyCut", false, "Apply occupancy cut"};
  // Configurable<int> occupancyCut{"occupancyCut", 1000, "Mimimum Occupancy cut"};

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};  // Combined
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};               // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};               // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                           // By pass TOF PID selection
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAeventPlots{"additionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> additionalMEPlots{"additionalMEPlots", false, "Additional Mixed event plots"};

  Configurable<bool> tofAtHighPt{"tofAtHighPt", false, "Use TOF at high pT"};
  Configurable<bool> additionalEvsel{"additionalEvsel", true, "Additional event selcection"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};

  // Rotational background
  Configurable<bool> isCalcRotBkg{"isCalcRotBkg", true, "Calculate rotational background"};
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<int> cNofRotations{"cNofRotations", 3, "Number of random rotations in the rotational background"};

  // MC Event selection
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  // Event selection cuts - Alex (Temporary, need to fix!)
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  // cuts on mother
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", false, "Enamble additional cuts on mother"};
  Configurable<double> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
  Configurable<double> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum Minv of mother cut"};
  Configurable<bool> cfgCutsOnDaughters{"cfgCutsOnDaughters", false, "Enamble additional cuts on daughters"};
  Configurable<int> cetaphiBins{"cetaphiBins", 400, "number of eta and phi bins"};
  Configurable<double> cMaxDeltaEtaCut{"cMaxDeltaEtaCut", 0.7, "Maximum deltaEta between daughters"};
  Configurable<double> cMaxDeltaPhiCut{"cMaxDeltaPhiCut", 1.5, "Maximum deltaPhi between daughters"};
  TRandom* rn = new TRandom();

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};
    // AxisSpec occupancyAxis = {occupancyBins, "Occupancy [-40,100]"};

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

    // Mass QA (quick check)
    if (invmass1D) {
      histos.add("k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTH1F, {invMassAxis});
      histos.add("k892invmassLS", "Invariant mass of K(892)0 like sign", kTH1F, {invMassAxis});
      histos.add("k892invmassME", "Invariant mass of K(892)0 mixed event", kTH1F, {invMassAxis});
      if (studyAntiparticle) {
        histos.add("k892invmassDSAnti", "Invariant mass of Anti-K(892)0 differnt sign", kTH1F, {invMassAxis});
        histos.add("k892invmassLSAnti", "Invariant mass of Anti-K(892)0 like sign", kTH1F, {invMassAxis});
      }
    }

    if (additionalMEPlots) {
      histos.add("k892invmassME_DS", "Invariant mass of K(892)0 mixed event DS", kTH1F, {invMassAxis});
      histos.add("k892invmassME_DSAnti", "Invariant mass of K(892)0 mixed event DSAnti", kTH1F, {invMassAxis});
    }

    if (additionalQAplots) {
      // TPC ncluster distirbutions
      histos.add("TPCncluster/TPCnclusterpi", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("TPCncluster/TPCnclusterka", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("TPCncluster/TPCnclusterPhipi", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
      histos.add("TPCncluster/TPCnclusterPhika", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

      // Multiplicity correlation calibrations
      histos.add("MultCalib/centglopi_before", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      // histos.add("MultCalib/centgloka_before", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      histos.add("MultCalib/GloPVpi_before", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
      // histos.add("MultCalib/GloPVka_before", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});

      histos.add("MultCalib/centglopi_after", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      // histos.add("MultCalib/centgloka_after", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      histos.add("MultCalib/GloPVpi_after", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
      // histos.add("MultCalib/GloPVka_after", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
    }

    // DCA QA
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAbefore/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    // pT QA
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAbefore/trkpT_ka", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_ka", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    // PID QA before cuts
    if (fillPidPlots) {
      histos.add("QAbefore/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_Nsigma_ka_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_Nsigma_ka_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    }

    // eta phi QA
    if (cfgCutsOnDaughters) {
      histos.add("QAbefore/deltaEta", "deltaEta of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAbefore/deltaPhi", "deltaPhi of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAbefore/Eta", "Eta of  tracks candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
      histos.add("QAbefore/Phi", "Phi of tracks candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});

      histos.add("QAafter/deltaEta", "deltaEta of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAafter/deltaPhi", "deltaPhi of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAafter/EtaPi", "Eta of  pion candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
      histos.add("QAafter/PhiPi", "Phi of  pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});
      histos.add("QAafter/EtaKa", "Eta of kaon  candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
      histos.add("QAafter/PhiKa", "Phi of kaon  candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});

      histos.add("QAafter/deltaEtaafter", "deltaEta of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAafter/deltaPhiafter", "deltaPhi of kaon and pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 3.15}});
      histos.add("QAafter/EtaPiafter", "Eta of  pion candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
      histos.add("QAafter/PhiPiafter", "Phi of  pion candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});
      histos.add("QAafter/EtaKaafter", "Eta of kaon  candidates", HistType::kTH1F, {{cetaphiBins, -1.6, 1.6}});
      histos.add("QAafter/PhiKaafter", "Phi of kaon  candidates", HistType::kTH1F, {{cetaphiBins, 0.0, 6.30}});
    }

    // 3d histogram
    histos.add("h3k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLS", "Invariant mass of K(892)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassME", "Invariant mass of K(892)0 mixed event", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLSAnti", "Invariant mass of Anti-K(892)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis});

    if (studyAntiparticle) {
      histos.add("h3k892invmassDSAnti", "Invariant mass of Anti-K(892)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
    }

    if (isCalcRotBkg) {
      histos.add("h3K892InvMassRotation", "Invariant mass of K(892)0 rotation", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
    }

    if (additionalMEPlots) {
      histos.add("h3k892invmassME_DS", "Invariant mass of K(892)0 mixed event DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassME_DSAnti", "Invariant mass of K(892)0 mixed event DSAnti", kTH3F, {centAxis, ptAxis, invMassAxis});
    }

    if (doprocessMCLight) {
      // MC QA
      histos.add("QAMCTrue/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMCTrue/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
      if (fillPidPlots) {
        histos.add("QAMCTrue/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TOF_Nsigma_ka_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      }
      histos.add("h3Reck892invmass", "Invariant mass of Reconstructed MC K(892)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
      histos.add("h3Reck892invmassAnti", "Invariant mass of Reconstructed MC Anti-K(892)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis});
      histos.add("k892Gen", "pT distribution of True MC K(892)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis});
      histos.add("k892GenAnti", "pT distribution of True MC Anti-K(892)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis});
      histos.add("k892Rec", "pT distribution of Reconstructed MC K(892)0", kTH2F, {ptAxis, centAxis});
      histos.add("k892RecAnti", "pT distribution of Reconstructed MC Anti-K(892)0", kTH2F, {ptAxis, centAxis});
      histos.add("k892Recinvmass", "Inv mass distribution of Reconstructed MC Phi", kTH1F, {invMassAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
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

  double massKa = MassKaonCharged;
  double massPi = MassPionCharged;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    // if (collision.alias_bit(kTVXinTRD)) {
    //   // TRD triggered
    //   // return 0;
    // }
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

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;

    return true;
  }
  // PID selection tools
  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    if (tofAtHighPt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
        tpcPIDPassed = true;
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
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
    }
    return false;
  }
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
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
        tpcPIDPassed = true;
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedKaon > 0) && (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa() < nsigmaCutCombinedKaon * nsigmaCutCombinedKaon)) {
          tofPIDPassed = true;
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

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    // auto multNTracksPV = collision.multNTracksPV();
    auto multiplicity = collision.cent();
    if (additionalEvsel && !eventSelected(collision, multiplicity)) {
      return;
    }
    // auto occupancyNo = collision.trackOccupancyInTimeRange();
    // if (applyOccupancyCut && occupancyNo < occupancyCut) {
    //   return;
    // }

    if (additionalQAplots) {
      histos.fill(HIST("MultCalib/centglopi_before"), multiplicity, dTracks1.size());            // centrality vs global tracks before the multiplicity calibration cuts
      histos.fill(HIST("MultCalib/GloPVpi_before"), dTracks1.size(), collision.multNTracksPV()); // global tracks vs PV tracks before the multiplicity calibration cuts
    }

    if (additionalQAeventPlots) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }

    if (additionalQAplots) {
      histos.fill(HIST("MultCalib/centglopi_after"), multiplicity, dTracks1.size());            // centrality vs global tracks after the multiplicity calibration cuts
      histos.fill(HIST("MultCalib/GloPVpi_after"), dTracks1.size(), collision.multNTracksPV()); // global tracks vs PV tracks after the multiplicity calibration cuts
    }

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance, ldaughterRot, lresonanceRot;
    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {

      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (additionalQAeventPlots) {
        if constexpr (!IsMix) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
        }
      }

      //// Initialize variables
      // Trk1: Pion, Trk2: Kaon
      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();
      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      auto deltaEta = std::abs(trk1.eta() - trk2.eta());
      auto deltaPhi = std::abs(trk1.phi() - trk2.phi());
      deltaPhi = (deltaPhi > constants::math::PI) ? (constants::math::TwoPI - deltaPhi) : deltaPhi;

      if constexpr (!IsMix) {
        //// QA plots before the selection
        //  --- PID QA Pion
        if (fillPidPlots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
          //  --- PID QA Kaon
          histos.fill(HIST("QAbefore/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
          }
        }
        histos.fill(HIST("QAbefore/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QAbefore/trkpT_ka"), trk2ptKa);
        histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAxy_ka"), trk2.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QAbefore/trkDCAz_ka"), trk2.dcaZ());

        if (cfgCutsOnDaughters) {
          histos.fill(HIST("QAbefore/Eta"), trk1.eta());
          histos.fill(HIST("QAbefore/Phi"), trk1.phi());
          histos.fill(HIST("QAbefore/deltaEta"), deltaEta);
          histos.fill(HIST("QAbefore/deltaPhi"), deltaPhi);
        }
      }

      // Apply the selection
      if (cUseOnlyTOFTrackPi && !isTrk1hasTOF)
        continue;
      if (cUseOnlyTOFTrackKa && !isTrk2hasTOF)
        continue;
      if (!selectionPIDPion(trk1) || !selectionPIDKaon(trk2))
        continue;

      if (additionalQAplots) {
        // TPCncluster distributions
        histos.fill(HIST("TPCncluster/TPCnclusterpi"), trk1.tpcNClsFound());
        histos.fill(HIST("TPCncluster/TPCnclusterka"), trk2.tpcNClsFound());
        histos.fill(HIST("TPCncluster/TPCnclusterPhipi"), trk1.tpcNClsFound(), trk1.phi());
        histos.fill(HIST("TPCncluster/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());
      }

      if constexpr (!IsMix) {
        //// QA plots after the selection
        //  --- PID QA Pion
        if (fillPidPlots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAafter/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
          //  --- PID QA Kaon
          histos.fill(HIST("QAafter/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QAafter/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
          }
        }
        histos.fill(HIST("QAafter/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QAafter/trkpT_ka"), trk2ptKa);
        histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QAafter/trkDCAxy_ka"), trk2.dcaXY());
        histos.fill(HIST("QAafter/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QAafter/trkDCAz_ka"), trk2.dcaZ());

        if (cfgCutsOnDaughters) {
          histos.fill(HIST("QAafter/EtaPi"), trk1.eta());
          histos.fill(HIST("QAafter/PhiPi"), trk1.phi());
          histos.fill(HIST("QAafter/EtaKa"), trk2.eta());
          histos.fill(HIST("QAafter/PhiKa"), trk2.phi());

          histos.fill(HIST("QAafter/deltaEta"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhi"), deltaPhi);
        }
      }

      //// Resonance reconstruction
      lDecayDaughter1.SetPtEtaPhiM(trk1.pt(), trk1.eta(), trk1.phi(), massPi);
      lDecayDaughter2.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (std::abs(lResonance.Rapidity()) >= 0.5)
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
          histos.fill(HIST("QAafter/EtaPiafter"), trk1.eta());
          histos.fill(HIST("QAafter/PhiPiafter"), trk1.phi());
          histos.fill(HIST("QAafter/EtaKaafter"), trk2.eta());
          histos.fill(HIST("QAafter/PhiKaafter"), trk2.phi());
          histos.fill(HIST("QAafter/deltaEtaafter"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhiafter"), deltaPhi);
        }
      }

      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          if (isCalcRotBkg) {
            for (int i = 0; i < cNofRotations; i++) {
              float theta2 = rn->Uniform(constants::math::PI - constants::math::PI / rotationalCut, constants::math::PI + constants::math::PI / rotationalCut);
              ldaughterRot.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi() + theta2, massKa); // for rotated background
              lresonanceRot = lDecayDaughter1 + ldaughterRot;
              histos.fill(HIST("h3K892InvMassRotation"), multiplicity, lresonanceRot.Pt(), lresonanceRot.M());
            }
          }
          if (studyAntiparticle) {
            if (trk1.sign() < 0) {
              if (invmass1D)
                histos.fill(HIST("k892invmassDS"), lResonance.M());
              histos.fill(HIST("h3k892invmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (trk1.sign() > 0) {
              if (invmass1D)
                histos.fill(HIST("k892invmassDSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassDSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          } else {
            if (invmass1D)
              histos.fill(HIST("k892invmassDS"), lResonance.M());
            histos.fill(HIST("h3k892invmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        } else {
          if (invmass1D)
            histos.fill(HIST("k892invmassME"), lResonance.M());
          histos.fill(HIST("h3k892invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (additionalMEPlots) {
            if (trk1.sign() < 0) {
              if (invmass1D)
                histos.fill(HIST("k892invmassME_DS"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (trk1.sign() > 0) {
              if (invmass1D)
                histos.fill(HIST("k892invmassME_DSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          }
        }

        // MC
        if constexpr (IsMC) {
          if (std::abs(trk1.pdgCode()) != 211 || std::abs(trk2.pdgCode()) != 321)
            continue;
          if (trk1.motherId() != trk2.motherId()) // Same mother
            continue;
          if (std::abs(trk1.motherPDG()) != 313)
            continue;

          // Track selection check.
          histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), trk1.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAxy_ka"), trk2.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAz_pi"), trk1.dcaZ());
          histos.fill(HIST("QAMCTrue/trkDCAz_ka"), trk2.dcaZ());
          if (fillPidPlots) {
            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
            if (isTrk1hasTOF) {
              histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            }
            histos.fill(HIST("QAMCTrue/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
            if (isTrk2hasTOF) {
              histos.fill(HIST("QAMCTrue/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
            }
          }

          // MC histograms
          if (trk1.motherPDG() < 0) {
            histos.fill(HIST("k892Rec"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("k892Recinvmass"), lResonance.M());
            histos.fill(HIST("h3Reck892invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("k892RecAnti"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("k892Recinvmass"), lResonance.M());
            histos.fill(HIST("h3Reck892invmassAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      } else if (trk1.sign() * trk2.sign() > 0) {
        if constexpr (!IsMix) {
          if (trk1.sign() < 0) {
            if (invmass1D)
              histos.fill(HIST("k892invmassLS"), lResonance.M());
            histos.fill(HIST("h3k892invmassLS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (trk1.sign() > 0) {
            if (invmass1D)
              histos.fill(HIST("k892invmassLSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassLSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processDataLight(aod::ResoCollision const& collision,
                        aod::ResoTracks const& resotracks)
  {
    // LOG(info) << "new collision, zvtx: " << collision.posZ();
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K892analysis, processDataLight, "Process Event for data", false);

  void processMCLight(ResoMCCols::iterator const& collision,
                      soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    if (!collision.isInAfterAllCuts() || (std::abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;
    fillHistograms<true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K892analysis, processMCLight, "Process Event for MC (Reconstructed)", false);

  void processMCTrue(ResoMCCols::iterator const& collision, aod::ResoMCParents const& resoParents)
  {
    auto multiplicity = collision.cent();
    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != 313 || std::abs(part.y()) >= 0.5)
        continue;
      bool pass1 = std::abs(part.daughterPDG1()) == 211 || std::abs(part.daughterPDG2()) == 211;
      bool pass2 = std::abs(part.daughterPDG1()) == 321 || std::abs(part.daughterPDG2()) == 321;

      if (!pass1 || !pass2)
        continue;

      if (collision.isVtxIn10()) // INEL10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("k892Gen"), 0, part.pt(), multiplicity);
        else
          histos.fill(HIST("k892GenAnti"), 0, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isInSel8()) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("k892Gen"), 1, part.pt(), multiplicity);
        else
          histos.fill(HIST("k892GenAnti"), 1, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isTriggerTVX()) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("k892Gen"), 2, part.pt(), multiplicity);
        else
          histos.fill(HIST("k892GenAnti"), 2, part.pt(), multiplicity);
      }
      if (collision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("k892Gen"), 3, part.pt(), multiplicity);
        else
          histos.fill(HIST("k892GenAnti"), 3, part.pt(), multiplicity);
      }
    }
  }
  PROCESS_SWITCH(K892analysis, processMCTrue, "Process Event for MC (Generated)", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMELight(o2::aod::ResoCollisions const& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K892analysis, processMELight, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K892analysis>(cfgc)};
}
