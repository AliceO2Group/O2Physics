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

/// \author Youpeng Su (yousu@cern.ch)

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include <TLorentzVector.h>
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include <TTree.h>
#include <TFile.h>
#include <TMatrixD.h>
#include "PWGLF/DataModel/lambdaJetpolarization.h"

using std::cout;
using std::endl;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myAnalysis {
  Produces<o2::aod::MyTable> myTable;
  Produces<o2::aod::MyTableJet> myTableJet;
  Produces<o2::aod::MyCollision> outputCollisions;
  Produces<o2::aod::MyCollisionV0s> outputCollisionsV0;
  Produces<o2::aod::MyTableLeadingJet> myTableLeadingJet;

  HistogramRegistry registry{"registry"};
  Configurable<float> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.05, "DCA Pos To PV"};
  SliceCache cache;
  HistogramRegistry JEhistos{"JEhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> cfgeventSelections{"cfgeventSelections", "sel8", "choose event selection"};
  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};
  Configurable<bool> cfgDataHists{"cfgDataHists", true, "Enables DataHists"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  // Others configure
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};
  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> cfgjetR{"cfgjetR", 0.4, "jet resolution parameter"};
  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.8, "set track max Eta"};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<int> cfgnTPCPID{"cfgnTPCPID", 4, "nTPC PID"};
  Configurable<int> cfgnTOFPID{"cfgnTOFPID", 4, "nTOF PID"};
  // V0 track selection////////////////////////////////////////////////////////////////
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireTPC{"requireTPC", true, "require TPC hit"};
  Configurable<bool> requirepassedSingleTrackSelection{"requirepassedSingleTrackSelection", false, "requirepassedSingleTrackSelection"};
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min track"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max track"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.2f, "Minimum V0 Radius"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -5.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +5.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -5.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +5.0f, "Maximum nsigma TOF"};
  Configurable<float> yMin{"yMin", -0.5f, "minimum y"};
  Configurable<float> yMax{"yMax", +0.5f, "maximum y"};
  Configurable<float> v0rejLambda{"v0rejLambda", 0.01, "V0 rej K0s"};
  Configurable<float> CtauLambda{"ctauLambda", 30, "C tau Lambda (cm)"};
  Configurable<bool> ifpasslambda{"passedLambdaSelection", 1, "passedLambdaSelection"};
  Configurable<bool> ifpassantilambda{"passedANtiLambdaSelection", 1, "passedAntiLambdaSelection"};
  // Event Selection/////////////////////////////////
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};
  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "TVX trigger"};
  Configurable<bool> iscutzvertex{"iscutzvertex", 1, "Accepted z-vertex range (cm)"};
  Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", 1, "TF border cut"};
  Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", 1, "ITS ROF border cut"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", 1, "Is Vertex TOF matched"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", 0, "isGoodZvtxFT0vsPV"};
  /////////////////////////V0 QA analysis///////////////////////////////
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<bool> doArmenterosCut{"doArmenterosCut", 1, "do Armenteros Cut"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  // CONFIG DONE
  /////////////////////////////////////////  //INIT////////////////////////////////////////////////////////////////////
  // int eventSelection = -1;
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec axisPt{200, 0, +200, "#pt"};
    const AxisSpec MinvAxis = {500, 0.1, 1.25};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec MultAxis = {100, 0, 100};
    const AxisSpec dRAxis = {100, 0, 100};

    const AxisSpec axisPx{200, -10, 10, "#px"};
    const AxisSpec axisPy{200, -10, 10, "#py"};
    const AxisSpec axisPz{200, -10, 10, "#pz"};
    const AxisSpec massAxis{200, 0.9f, 1.2f, "mass"};
    const AxisSpec eventAxis{1000000, 0.5f, 1000000.5f, "event"};

    if (cfgDataHists) {

      JEhistos.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", kTH1F, {{200, 0., 200.}});
      JEhistos.add("h_track_eta", "track #eta;#eta_{track};entries", kTH1F, {{100, -1.f, 1.f}});
      JEhistos.add("h_track_phi", "track #varphi;#varphi_{track};entries", kTH1F, {{80, -1.f, 7.f}});
      JEhistos.add("nJetsPerEvent", "nJetsPerEvent", kTH1F, {{10, 0.0, 10.0}});
      JEhistos.add("FJetaHistogram", "FJetaHistogram", kTH1F, {axisEta});
      JEhistos.add("FJphiHistogram", "FJphiHistogram", kTH1F, {axisPhi});
      JEhistos.add("FJptHistogram", "FJptHistogram", kTH1F, {axisPt});
      JEhistos.add("hDCArToPv", "DCArToPv", kTH1F, {{300, 0.0, 3.0}});
      JEhistos.add("hDCAzToPv", "DCAzToPv", kTH1F, {{300, 0.0, 3.0}});
      JEhistos.add("rawpT", "rawpT", kTH1F, {{1000, 0.0, 10.0}});
      JEhistos.add("rawDpT", "rawDpT", kTH2F, {{1000, 0.0, 10.0}, {300, -1.5, 1.5}});
      JEhistos.add("hIsPrim", "hIsPrim", kTH1F, {{2, -0.5, +1.5}});
      JEhistos.add("hIsGood", "hIsGood", kTH1F, {{2, -0.5, +1.5}});
      JEhistos.add("hIsPrimCont", "hIsPrimCont", kTH1F, {{2, -0.5, +1.5}});
      JEhistos.add("hFindableTPCClusters", "hFindableTPCClusters", kTH1F, {{200, 0, 200}});
      JEhistos.add("hFindableTPCRows", "hFindableTPCRows", kTH1F, {{200, 0, 200}});
      JEhistos.add("hClustersVsRows", "hClustersVsRows", kTH1F, {{200, 0, 2}});
      JEhistos.add("hTPCChi2", "hTPCChi2", kTH1F, {{200, 0, 100}});
      JEhistos.add("hITSChi2", "hITSChi2", kTH1F, {{200, 0, 100}});
      JEhistos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
      JEhistos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
      JEhistos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
      JEhistos.add("V0Counts", "V0Counts", kTH1F, {{10, 0, 10}});
      JEhistos.add("hUSS_1D", "hUSS_1D", kTH1F, {MinvAxis});
      JEhistos.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
      JEhistos.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
      JEhistos.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
      JEhistos.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 0.9f, 1.2f}}});
      JEhistos.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 0.9f, 1.2f}}});
      JEhistos.add("V0Radius", "V0Radius", {HistType::kTH1D, {{100, 0.0f, 20.0f}}});
      JEhistos.add("CosPA", "CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
      JEhistos.add("V0DCANegToPV", "V0DCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      JEhistos.add("V0DCAPosToPV", "V0DCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      JEhistos.add("V0DCAV0Daughters", "V0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
      JEhistos.add("TPCNSigmaPosPi", "TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      JEhistos.add("TPCNSigmaNegPi", "TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      JEhistos.add("TPCNSigmaPosPr", "TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      JEhistos.add("TPCNSigmaNegPr", "TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      JEhistos.add("hNEvents", "hNEvents", {HistType::kTH1I, {{10, 0.f, 10.f}}});
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isGoodZvtxFT0vsPV");
      JEhistos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "Applied selected");
      registry.add("hNEventsJet", "hNEventsJet", {HistType::kTH1I, {{4, 0.f, 4.f}}});
      registry.get<TH1>(HIST("hNEventsJet"))->GetXaxis()->SetBinLabel(1, "all");
      registry.get<TH1>(HIST("hNEventsJet"))->GetXaxis()->SetBinLabel(2, "zvertex");
      registry.get<TH1>(HIST("hNEventsJet"))->GetXaxis()->SetBinLabel(3, "JCollisionSel::sel8");
      JEhistos.add("v0Lambdapx", "v0Lambdapx", kTH1F, {axisPx});
      JEhistos.add("v0Lambdapy", "v0Lambdapy", kTH1F, {axisPy});
      JEhistos.add("v0Lambdapz", "v0Lambdapz", kTH1F, {axisPz});
      JEhistos.add("v0AntiLambdapx", "v0AntiLambdapx", kTH1F, {axisPx});
      JEhistos.add("v0AntiLambdapy", "v0AntiLambdapy", kTH1F, {axisPy});
      JEhistos.add("v0AntiLambdapz", "v0AntiLambdapz", kTH1F, {axisPz});
      JEhistos.add("jetpx", "jetpx", kTH1F, {axisPx});
      JEhistos.add("jetpy", "jetpy", kTH1F, {axisPy});
      JEhistos.add("jetpz", "jetpz", kTH1F, {axisPz});
      JEhistos.add("hV0Lambda", "V0Lambda",
                   {HistType::kTHnSparseF, {eventAxis, axisPx, axisPy, axisPz, massAxis, axisPx, axisPy, axisPz}});
      JEhistos.add("EventIndexselection", "EventIndexselection", {HistType::kTH1F, {{1000000, 0.5f, 1000000.5f}}});
    }
    JEhistos.add("hKaonplusCounts", "hKaonplusCounts", {HistType::kTH1F, {{1, -0.5, 0.5f}}});
    JEhistos.add("hKaonminusCounts", "hKaonminusCounts", {HistType::kTH1F, {{1, -0.5, 0.5f}}});

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(cfgeventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  } // end of init

  double massPi = o2::constants::physics::MassPiMinus;
  double massPr = o2::constants::physics::MassProton;
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr, aod::TrackSelection>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr>;
  using JCollisions = soa::Join<aod::JCollisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFV0As>;
  using V0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFV0As>;

  Filter jetCuts = aod::jet::pt > cfgjetPtMin&& aod::jet::r == nround(cfgjetR.node() * 100.0f);
  template <typename TrackType>
  bool TrackSelection(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cfgtrkMinPt)
      return false;

    if (std::abs(track.eta()) > cfgtrkMaxEta)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;

    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows)
      return false;

    if (track.tpcCrossedRowsOverFindableCls() > cfgnRowsOverFindable)
      return false;

    if (track.tpcChi2NCl() > cfgnTPCChi2)
      return false;

    if (track.itsChi2NCl() > cfgnITSChi2)
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;

    return true;
  };

  template <typename T>
  bool trackPIDPion(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cfgnTOFPID) {
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

  template <typename T>
  bool trackPIDProton(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPr()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPr()) < cfgnTOFPID) {
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

  // init Selection
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedInitLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    if (v0.v0radius() < minimumV0Radius || v0.v0cosPA() < v0cospa ||
        TMath::Abs(ptrack.eta()) > etaMax ||
        TMath::Abs(ntrack.eta()) > etaMax) {
      return false;
    }
    if (v0.dcaV0daughters() > dcav0dau) {
      return false;
    }

    if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv) {
      return false;
    }

    if (TMath::Abs(v0.dcapostopv()) < dcapostopv) {
      return false;
    }
    return true;
  }

  // Lambda Selections
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {

    if (ptrack.tpcNSigmaKa() > nsigmaTPCmin && ptrack.tpcNSigmaKa() < nsigmaTPCmax) {
      JEhistos.fill(HIST("hKaonplusCounts"), 1);
    }

    if (ntrack.tpcNSigmaKa() > nsigmaTPCmin && ntrack.tpcNSigmaKa() < nsigmaTPCmax) {
      JEhistos.fill(HIST("hKaonminusCounts"), 1);
    }

    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0cosPA() < v0cospa ||
        TMath::Abs(ptrack.eta()) > etaMax ||
        TMath::Abs(ntrack.eta()) > etaMax) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv)
      return false;
    if (TMath::Abs(v0.dcapostopv()) < dcapostopv)
      return false;
    if (v0.dcaV0daughters() > dcav0dau)
      return false;

    // PID Selections (TPC)
    if (requireTPC) {
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;
    }

    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);

    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax) {
      return false;
    }

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) > 0.075) {
      return false;
    }

    return true;
  }
  // AntiLambda Selections
  template <typename AntiLambda, typename TrackPos, typename TrackNeg>
  bool passedAntiLambdaSelection(const AntiLambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0cosPA() < v0cospa ||
        TMath::Abs(ptrack.eta()) > etaMax ||
        TMath::Abs(ntrack.eta()) > etaMax) {
      return false;
    }

    if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv) //
      return false;
    if (TMath::Abs(v0.dcapostopv()) < dcapostopv) //
      return false;
    if (v0.dcaV0daughters() > dcav0dau) //
      return false;
    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }
    if (requireTPC) {
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
    }

    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax) {
      return false;
    }

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > 0.075) {
      return false;
    }
    return true;
  }
  ///////Event selection
  template <typename TCollision>
  bool AcceptEvent(TCollision const& collision)
  {
    if (sel8 && !collision.sel8()) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 1.5);

    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 2.5);

    if (iscutzvertex && TMath::Abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 3.5);

    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    JEhistos.fill(HIST("hNEvents"), 4.5);

    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 5.5);
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 6.5);
    if (isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    JEhistos.fill(HIST("hNEvents"), 7.5);

    return true;
  }
  // Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcaV0DaughtersMax;

  int nEventsJet = 0;
  void processJetTracks(aod::JetCollision const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedjets, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, TrackCandidates const&)
  {

    registry.fill(HIST("hNEventsJet"), 0.5);
    if (fabs(collision.posZ()) > cfgVtxCut) {

      return;
    }
    registry.fill(HIST("hNEventsJet"), 1.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("hNEventsJet"), 2.5);
    outputCollisions(collision.posZ());

    for (auto const& track : tracks) {

      auto originalTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr>>();
      JEhistos.fill(HIST("hDCArToPv"), originalTrack.dcaXY());
      JEhistos.fill(HIST("hDCAzToPv"), originalTrack.dcaZ());
      JEhistos.fill(HIST("rawpT"), originalTrack.pt());
      JEhistos.fill(HIST("rawDpT"), track.pt(), track.pt() - originalTrack.pt());
      JEhistos.fill(HIST("hIsPrim"), originalTrack.isPrimaryTrack());
      JEhistos.fill(HIST("hIsGood"), originalTrack.isGlobalTrackWoDCA());
      JEhistos.fill(HIST("hIsPrimCont"), originalTrack.isPVContributor());
      JEhistos.fill(HIST("hFindableTPCClusters"), originalTrack.tpcNClsFindable());
      JEhistos.fill(HIST("hFindableTPCRows"), originalTrack.tpcNClsCrossedRows());
      JEhistos.fill(HIST("hClustersVsRows"), originalTrack.tpcCrossedRowsOverFindableCls());
      JEhistos.fill(HIST("hTPCChi2"), originalTrack.tpcChi2NCl());
      JEhistos.fill(HIST("hITSChi2"), originalTrack.itsChi2NCl());
      JEhistos.fill(HIST("h_track_pt"), track.pt());
      JEhistos.fill(HIST("h_track_eta"), track.eta());
      JEhistos.fill(HIST("h_track_phi"), track.phi());

      if (track.pt() < cfgtrkMinPt && std::abs(track.eta()) > cfgtrkMaxEta) {
        continue;
      }
      JEhistos.fill(HIST("ptHistogram"), track.pt());
      JEhistos.fill(HIST("etaHistogram"), track.eta());
      JEhistos.fill(HIST("phiHistogram"), track.phi());
    }
    int nJets = 0;
    int lastindex = 0;
    int collisionId = 0;
    float maxJetpx = 0;
    float maxJetpy = 0;
    float maxJetpz = 0;
    float maxJetpT = 0;
    float maxJetPt = -999;
    nEventsJet++;
    for (const auto& chargedjet : chargedjets) {
      JEhistos.fill(HIST("FJetaHistogram"), chargedjet.eta());
      JEhistos.fill(HIST("FJphiHistogram"), chargedjet.phi());
      JEhistos.fill(HIST("FJptHistogram"), chargedjet.pt());

      JEhistos.fill(HIST("jetpx"), chargedjet.px());
      JEhistos.fill(HIST("jetpy"), chargedjet.py());
      JEhistos.fill(HIST("jetpz"), chargedjet.pz());

      myTableJet(outputCollisions.lastIndex(), chargedjet.collisionId(), chargedjet.px(), chargedjet.py(), chargedjet.pz(), chargedjet.pt());

      nJets++;
      if (chargedjet.pt() > maxJetPt) {
        maxJetpx = chargedjet.px();
        maxJetpy = chargedjet.py();
        maxJetpz = chargedjet.pz();
        maxJetpT = chargedjet.pt();
        collisionId = chargedjet.collisionId();
        lastindex = outputCollisions.lastIndex();
      }
    }
    if (maxJetpT > 0) {
      myTableLeadingJet(lastindex, collisionId, maxJetpx, maxJetpy, maxJetpz, maxJetpT);
    }
    JEhistos.fill(HIST("nJetsPerEvent"), nJets);
  }
  PROCESS_SWITCH(myAnalysis, processJetTracks, "process JE Framework", true);
  int nEventsV0 = 0;

  void processV0(V0Collisions::iterator const& collision, aod::V0Datas const& V0s, TrackCandidates const&)
  {
    nEventsV0++;
    JEhistos.fill(HIST("hNEvents"), 0.5);
    if (!AcceptEvent(collision)) {
      return;
    }
    JEhistos.fill(HIST("hNEvents"), 8.5);
    int V0NumbersPerEvent = 0;
    int V0LambdaNumbers = 0;
    outputCollisionsV0(collision.posX());
    for (auto& v0 : V0s) {
      float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;

      const auto& pos = v0.posTrack_as<TrackCandidates>();
      const auto& neg = v0.negTrack_as<TrackCandidates>();
      V0NumbersPerEvent = V0NumbersPerEvent + 1;
      if (passedLambdaSelection(v0, pos, neg) && ctauLambda < CtauLambda && ifpasslambda) {
        JEhistos.fill(HIST("hPt"), v0.pt());
        JEhistos.fill(HIST("V0Radius"), v0.v0radius());
        JEhistos.fill(HIST("CosPA"), v0.v0cosPA());
        JEhistos.fill(HIST("V0DCANegToPV"), v0.dcanegtopv());
        JEhistos.fill(HIST("V0DCAPosToPV"), v0.dcapostopv());
        JEhistos.fill(HIST("V0DCAV0Daughters"), v0.dcaV0daughters());
      } else if (passedInitLambdaSelection(v0, pos, neg)) {
        JEhistos.fill(HIST("hPt"), v0.pt());
        JEhistos.fill(HIST("V0Radius"), v0.v0radius());
        JEhistos.fill(HIST("CosPA"), v0.v0cosPA());
        JEhistos.fill(HIST("V0DCANegToPV"), v0.dcanegtopv());
        JEhistos.fill(HIST("V0DCAPosToPV"), v0.dcapostopv());
        JEhistos.fill(HIST("V0DCAV0Daughters"), v0.dcaV0daughters());
      }
      if (passedLambdaSelection(v0, pos, neg) && ctauAntiLambda < CtauLambda && ifpasslambda) {

        V0LambdaNumbers = V0LambdaNumbers + 1;
        JEhistos.fill(HIST("hMassVsPtLambda"), v0.pt(), v0.mLambda());
        JEhistos.fill(HIST("hMassLambda"), v0.mLambda());
        JEhistos.fill(HIST("TPCNSigmaPosPr"), pos.tpcNSigmaPr());
        JEhistos.fill(HIST("TPCNSigmaNegPi"), neg.tpcNSigmaPi());

        JEhistos.fill(HIST("v0Lambdapx"), v0.px());
        JEhistos.fill(HIST("v0Lambdapy"), v0.py());
        JEhistos.fill(HIST("v0Lambdapz"), v0.pz());
        myTable(outputCollisionsV0.lastIndex(), v0.collisionId(), v0.px(), v0.py(), v0.pz(), v0.pt(), v0.mLambda(), pos.px(), pos.py(), pos.pz());
        JEhistos.fill(HIST("hV0Lambda"), nEventsV0, v0.px(), v0.py(), v0.pz(), v0.mLambda(), pos.px(), pos.py(), pos.pz());
      } else if (passedInitLambdaSelection(v0, pos, neg)) {
        V0LambdaNumbers = V0LambdaNumbers + 1;
        JEhistos.fill(HIST("hMassVsPtLambda"), v0.pt(), v0.mLambda());
        JEhistos.fill(HIST("hMassLambda"), v0.mLambda());
        JEhistos.fill(HIST("TPCNSigmaPosPr"), pos.tpcNSigmaPr());
        JEhistos.fill(HIST("TPCNSigmaNegPi"), neg.tpcNSigmaPi());

        JEhistos.fill(HIST("v0Lambdapx"), v0.px());
        JEhistos.fill(HIST("v0Lambdapy"), v0.py());
        JEhistos.fill(HIST("v0Lambdapz"), v0.pz());
        myTable(outputCollisionsV0.lastIndex(), v0.collisionId(), v0.px(), v0.py(), v0.pz(), v0.pt(), v0.mLambda(), pos.px(), pos.py(), pos.pz());
        JEhistos.fill(HIST("hV0Lambda"), nEventsV0, v0.px(), v0.py(), v0.pz(), v0.mLambda(), pos.px(), pos.py(), pos.pz());
      }
      if (passedAntiLambdaSelection(v0, pos, neg) && ifpassantilambda) {

        JEhistos.fill(HIST("hMassVsPtAntiLambda"), v0.pt(), v0.mAntiLambda());
        JEhistos.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
        JEhistos.fill(HIST("TPCNSigmaPosPi"), pos.tpcNSigmaPi());
        JEhistos.fill(HIST("TPCNSigmaNegPr"), neg.tpcNSigmaPr());

        JEhistos.fill(HIST("v0AntiLambdapx"), v0.px());
        JEhistos.fill(HIST("v0AntiLambdapy"), v0.py());
        JEhistos.fill(HIST("v0AntiLambdapz"), v0.pz());
      } else if (passedInitLambdaSelection(v0, pos, neg)) {
        JEhistos.fill(HIST("hMassVsPtAntiLambda"), v0.pt(), v0.mAntiLambda());
        JEhistos.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
        JEhistos.fill(HIST("TPCNSigmaPosPi"), pos.tpcNSigmaPi());
        JEhistos.fill(HIST("TPCNSigmaNegPr"), neg.tpcNSigmaPr());

        JEhistos.fill(HIST("v0AntiLambdapx"), v0.px());
        JEhistos.fill(HIST("v0AntiLambdapy"), v0.py());
        JEhistos.fill(HIST("v0AntiLambdapz"), v0.pz());
      }
    }
    JEhistos.fill(HIST("V0Counts"), V0NumbersPerEvent);
  }
  PROCESS_SWITCH(myAnalysis, processV0, "processV0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myAnalysis>(cfgc),
  };
}
