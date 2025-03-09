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

/// \file phiInJets.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance hadrochemistry analysis.
///
///
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <algorithm>
#include <iostream>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phiInJets {
  SliceCache cache;
  HistogramRegistry JEhistos{"JEhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> cfgeventSelections{"cfgeventSelections", "sel8", "choose event selection"};
  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};

  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.9, "set track max Eta"};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgnTPCPID{"cfgnTPCPID", 4, "nTPC PID"};
  Configurable<int> cfgnTOFPID{"cfgnTOFPID", 4, "nTOF PID"};
  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> cfgjetR{"cfgjetR", 0.4, "jet resolution parameter"};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};
  Configurable<bool> cfgBR{"cfgBR", false, "Forces Gen. Charged BR Only"};
  Configurable<bool> cfgSimPID{"cfgSimPID", false, "Enforces PID on the Gen. Rec level"};
  Configurable<bool> cfgSingleJet{"cfgSingleJet", false, "Enforces strict phi-jet correspondance"};
  Configurable<bool> cfgIsKstar{"cfgIsKstar", false, "Swaps Phi for Kstar analysis"};
  Configurable<bool> cfgDataHists{"cfgDataHists", false, "Enables DataHists"};
  Configurable<bool> cfgMCRecHists{"cfgMCRecHists", false, "Enables MCRecHists"};
  Configurable<bool> cfgMCGenHists{"cfgMCGenHists", false, "Enables MCGenHists"};
  Configurable<bool> cfgMCGenMATCHEDHists{"cfgMCGenMATCHEDHists", false, "Enables MCGenMATCHEDHists"};
  Configurable<bool> cfgMCRecMATCHEDHists{"cfgMCRecMATCHEDHists", false, "Enables MCRecMATCHEDHists"};

  // CONFIG DONE
  /////////////////////////////////////////  //INIT

  std::vector<int> eventSelectionBits;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec axisPt{200, 0, +200, "#pt"};
    const AxisSpec MinvAxis = {500, 0.75, 1.25};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec MultAxis = {100, 0, 100};
    const AxisSpec dRAxis = {100, 0, 100};

    // Now we define histograms. Due to saving space, we only reserve histogram memory for histograms we will be filling.

    if (cfgDataHists) {
      JEhistos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
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

      JEhistos.add("hNResoPerEvent", "hNResoPerEvent", kTH1F, {{10, 0, 10}});
      JEhistos.add("hNResoPerEventWJet", "hNResoPerEventWJet", kTH1F, {{10, 0, 10}});
      JEhistos.add("hNResoPerEventInJet", "hNResoPerEventInJet", kTH1F, {{10, 0, 10}});

      JEhistos.add("FJetaHistogram", "FJetaHistogram", kTH1F, {axisEta});
      JEhistos.add("FJphiHistogram", "FJphiHistogram", kTH1F, {axisPhi});
      JEhistos.add("FJptHistogram", "FJptHistogram", kTH1F, {axisPt});
      JEhistos.add("nJetsPerEvent", "nJetsPerEvent", kTH1F, {{10, 0.0, 10.0}});

      JEhistos.add("hUSS", "hUSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hUSS_1D", "hUSS_1D", kTH1F, {MinvAxis});
      JEhistos.add("hUSS_1D_2_3", "hUSS_1D_2_3", kTH1F, {MinvAxis});
      JEhistos.add("hLSS", "hLSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hLSS_1D", "hLSS_1D", kTH1F, {MinvAxis});
      JEhistos.add("hLSS_1D_2_3", "hLSS_1D_2_3", kTH1F, {MinvAxis});
      JEhistos.add("hUSS_INSIDE", "hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hUSS_INSIDE_1D", "hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hUSS_INSIDE_1D_2_3", "hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
      JEhistos.add("hLSS_INSIDE", "hLSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hLSS_INSIDE_1D", "hLSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hLSS_INSIDE_1D_2_3", "hLSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    }

    if (cfgMCRecHists) {
      JEhistos.add("nEvents_MCRec", "nEvents_MCRec", kTH1F, {{4, 0.0, 4.0}});

      JEhistos.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{4000, 0., 200.}}});
      JEhistos.add("h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      JEhistos.add("h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});

      JEhistos.add("ptJEHistogramPion", "ptJEHistogramPion", kTH1F, {PtAxis});
      JEhistos.add("ptJEHistogramKaon", "ptJEHistogramKaon", kTH1F, {PtAxis});
      JEhistos.add("ptJEHistogramProton", "ptJEHistogramProton", kTH1F, {PtAxis});
      JEhistos.add("ptJEHistogramPhi", "ptJEHistogramPhi", kTH1F, {PtAxis});
      JEhistos.add("ptJEHistogramPhi_JetTrigger", "ptJEHistogramPhi_JetTrigger", kTH1F, {PtAxis});
      JEhistos.add("minvJEHistogramPhi", "minvJEHistogramPhi", kTH1F, {MinvAxis});
      JEhistos.add("hNRealPhiVPhiCand", "hNRealPhiVPhiCand", kTH2F, {{10, 0, 10}, {10, 0, 10}});
      JEhistos.add("hNRealPhiWithJetVPhiCand", "hNRealPhiWithJetVPhiCand", kTH2F, {{10, 0, 10}, {10, 0, 10}});
      JEhistos.add("hNRealPhiInJetVPhiCand", "hNRealPhiInJetVPhiCand", kTH2F, {{10, 0, 10}, {10, 0, 10}});

      JEhistos.add("hMCRec_nonmatch_hUSS_KtoKangle_v_pt", "hMCRec_nonmatch_hUSS_KtoKangle_v_pt", kTH2F, {axisEta, PtAxis});
      JEhistos.add("hMCRec_nonmatch_hUSS_Kangle_v_pt", "hMCRec_nonmatch_hUSS_Kangle_v_pt", kTH2F, {axisEta, PtAxis});
      JEhistos.add("hMCRec_nonmatch_hUSS_INSIDE_pt_v_eta", "hMCRec_nonmatch_hUSS_INSIDE_pt_v_eta", kTH2F, {PtAxis, axisEta});
      JEhistos.add("JetVsPhi_REC", "JetVsPhi_REC", kTH2F, {{4000, 0., 200.}, {200, 0, 20.0}});

      JEhistos.add("hMCRecTrue_hUSS", "hMCRecTrue_hUSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRecTrue_hLSS", "hMCRecTrue_hLSS", kTH3F, {dRAxis, PtAxis, MinvAxis});

      JEhistos.add("hMCRec_hUSS", "hMCRec_hUSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_hUSS_1D", "hMCRec_hUSS_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_hUSS_1D_2_3", "hMCRec_hUSS_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRec_hLSS", "hMCRec_hLSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_hLSS_1D", "hMCRec_hLSS_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_hLSS_1D_2_3", "hMCRec_hLSS_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRec_hUSS_INSIDE", "hMCRec_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_hUSS_INSIDE_1D", "hMCRec_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_hUSS_INSIDE_1D_2_3", "hMCRec_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRec_hLSS_INSIDE", "hMCRec_hLSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_hLSS_INSIDE_1D", "hMCRec_hLSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_hLSS_INSIDE_1D_2_3", "hMCRec_hLSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRecTrue_hUSS_INSIDE", "hMCRecTrue_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRecTrue_hUSS_INSIDE_1D", "hMCRecTrue_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRecTrue_hUSS_INSIDE_1D_2_3", "hMCRecTrue_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRecTrue_hLSS_INSIDE", "hMCRecTrue_hLSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRecTrue_hLSS_INSIDE_1D", "hMCRecTrue_hLSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRecTrue_hLSS_INSIDE_1D_2_3", "hMCRecTrue_hLSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

      JEhistos.add("hMCRec_nonmatch_hUSS_INSIDE", "hMCRec_nonmatch_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_nonmatch_hUSS_INSIDE_1D", "hMCRec_nonmatch_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_nonmatch_hUSS_INSIDE_1D_2_3", "hMCRec_nonmatch_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    }

    if (cfgMCGenHists) {
      JEhistos.add("nEvents_MCGen", "nEvents_MCGen", kTH1F, {{4, 0.0, 4.0}});

      JEhistos.add("h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{4000, 0., 200.}}});
      JEhistos.add("h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      JEhistos.add("h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});

      JEhistos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {PtAxis});
      JEhistos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {PtAxis});
      JEhistos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {PtAxis});
      JEhistos.add("ptGeneratedPhi", "ptGeneratedPhi", kTH1F, {PtAxis});
      JEhistos.add("ptGeneratedPhi_ALLBR", "ptGeneratedPhi_ALLBR", kTH1F, {PtAxis});

      JEhistos.add("ptGeneratedPhi_JetTrigger", "ptGeneratedPhi_JetTrigger", kTH1F, {PtAxis});
      JEhistos.add("mGeneratedPhi", "mGeneratedPhi", kTH1F, {MinvAxis});

      JEhistos.add("hMCTrue_nonmatch_hUSS_Kangle_v_pt", "hMCTrue_nonmatch_hUSS_Kangle_v_pt", kTH2F, {axisEta, PtAxis});
      JEhistos.add("hMCTrue_nonmatch_hUSS_INSIDE_pt_v_eta", "hMCTrue_nonmatch_hUSS_INSIDE_pt_v_eta", kTH2F, {PtAxis, axisEta});
      JEhistos.add("JetVsPhi_GEN", "JetVsPhi_GEN", kTH2F, {{4000, 0., 200.}, {200, 0, 20.0}});

      JEhistos.add("hMCTrue_nonmatch_hUSS_INSIDE", "hMCTrue_nonmatch_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCTrue_nonmatch_hUSS_INSIDE_1D", "hMCTrue_nonmatch_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCTrue_nonmatch_hUSS_INSIDE_1D_2_3", "hMCTrue_nonmatch_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    }

    if (cfgMCGenMATCHEDHists) {
      JEhistos.add("nEvents_MCGen_MATCHED", "nEvents_MCGen_MATCHED", kTH1F, {{4, 0.0, 4.0}});

      JEhistos.add("h_matched_GEN_jet_pt", "matched_GEN level jet pT;#it{p}_{T,jet part} (GeV/#it{c});Delta", {HistType::kTH2F, {{200, 0., 200.}, {400, -20., 20.}}});
      JEhistos.add("h_matched_GEN_jet_eta", "matched_GEN level jet #eta;#eta_{jet part};Delta", {HistType::kTH2F, {{100, -1.0, 1.0}, {400, -20., 20.}}});
      JEhistos.add("h_matched_GEN_jet_phi", "matched_GEN level jet #phi;#phi_{jet part};Delta", {HistType::kTH2F, {{80, -1.0, 7.}, {400, -20., 20.}}});
      JEhistos.add("2DGenToRec", "2DGenToRec", kTH2F, {PtAxis, axisPt});
      JEhistos.add("2DGenToRec_constrained", "2DGenToRec_constrained", kTH2F, {PtAxis, axisPt});

      JEhistos.add("RespGen_Matrix_MATCHED", "RespGen_Matrix_MATCHED", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt});             // REC(Phi,Jet), GEN(Phi,Jet)
      JEhistos.add("RespGen_Matrix_MATCHED_rand0", "RespGen_Matrix_MATCHED_rand0", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt}); // REC(Phi,Jet), GEN(Phi,Jet)
      JEhistos.add("RespGen_Matrix_MATCHED_rand1", "RespGen_Matrix_MATCHED_rand1", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt}); // REC(Phi,Jet), GEN(Phi,Jet)

      JEhistos.add("hMCTrue_hUSS_INSIDE", "hMCTrue_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCTrue_hUSS_INSIDE_1D", "hMCTrue_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCTrue_hUSS_INSIDE_1D_2_3", "hMCTrue_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    }

    if (cfgMCRecMATCHEDHists) {
      JEhistos.add("nEvents_MCRec_MATCHED", "nEvents_MCRec_MATCHED", kTH1F, {{4, 0.0, 4.0}});

      JEhistos.add("h_matched_REC_jet_pt", "matched_REC level jet pT;#it{p}_{T,jet part} (GeV/#it{c});Delta", {HistType::kTH2F, {{200, 0., 200.}, {400, -20., 20.}}});
      JEhistos.add("h_matched_REC_jet_eta", "matched_REC level jet #eta;#eta_{jet part};Delta", {HistType::kTH2F, {{100, -1.0, 1.0}, {400, -20., 20.}}});
      JEhistos.add("h_matched_REC_jet_phi", "matched_REC level jet #phi;#phi_{jet part};Delta", {HistType::kTH2F, {{80, -1.0, 7.}, {400, -20., 20.}}});

      JEhistos.add("Resp_Matrix_MATCHED", "Resp_Matrix_MATCHED", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt});             // REC(Phi,Jet), GEN(Phi,Jet)
      JEhistos.add("Resp_Matrix_MATCHED_rand0", "Resp_Matrix_MATCHED_rand0", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt}); // REC(Phi,Jet), GEN(Phi,Jet)
      JEhistos.add("Resp_Matrix_MATCHED_rand1", "Resp_Matrix_MATCHED_rand1", HistType::kTHnSparseD, {PtAxis, axisPt, PtAxis, axisPt}); // REC(Phi,Jet), GEN(Phi,Jet)

      JEhistos.add("2DRecToGen", "2DRecToGen", kTH2F, {PtAxis, axisPt});
      JEhistos.add("2DRecToGen_constrained", "2DRecToGen_constrained", kTH2F, {PtAxis, axisPt});

      JEhistos.add("hMCRec_hUSS_INSIDE", "hMCRec_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
      JEhistos.add("hMCRec_hUSS_INSIDE_1D", "hMCRec_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
      JEhistos.add("hMCRec_hUSS_INSIDE_1D_2_3", "hMCRec_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    }
    // EVENT SELECTION
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(cfgeventSelections));

  } // end of init

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>;
  Filter jetCuts = aod::jet::pt > cfgjetPtMin&& aod::jet::r == nround(cfgjetR.node() * 100.0f);

  // Function for track quality cuts
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename TrackType>
  bool trackSelection(const TrackType track)
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
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <typename T>
  bool trackPID(const T& candidate, bool FT)
  {
    bool pid = false;

    if (!cfgIsKstar) {
      pid = trackPIDKaon<T>(candidate);

    } else {
      if (!FT)
        pid = trackPIDPion<T>(candidate);
      else
        pid = trackPIDKaon<T>(candidate);
    }
    return pid;
  }

  template <typename T>
  bool trackPIDKaon(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfgnTOFPID) {
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

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <typename JetType>
  double DistinguishJets(const JetType& jets, const TLorentzVector lResonance)
  {
    if (cDebugLevel > 0)
      std::cout << "oof, multiple jets fit to the same phi. Time to find the best phi-jet link" << std::endl;

    double best_R = 0;
    double best_jetpt = 0;
    for (auto const& jet : jets) {
      double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
      double etadiff = jet.eta() - lResonance.Eta();
      double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
      if (R < cfgjetR && best_R == 0) {
        best_R = R;
        best_jetpt = jet.pt();
      } else if (R < best_R) {
        best_R = R;
        best_jetpt = jet.pt();
      }
    } // jetloop
    return best_jetpt;
  }

  template <typename Jet_pt, typename Jet_phi, typename Jet_eta>
  double DistinguishJetsMC(const Jet_pt& jet_pt, const Jet_phi& jet_phi, const Jet_eta& jet_eta, const TLorentzVector lResonance)
  {
    if (cDebugLevel > 0)
      std::cout << "oof, multiple jets fit to the same phi. Time to find the best phi-jet link" << std::endl;

    double best_R = 0;
    double best_jetpt = 0;
    for (std::vector<double>::size_type i = 0; i < jet_pt.size(); i++) {
      double phidiff = TVector2::Phi_mpi_pi(jet_phi[i] - lResonance.Phi());
      double etadiff = jet_eta[i] - lResonance.Eta();
      double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
      if (R < cfgjetR && best_R == 0) {
        best_R = R;
        best_jetpt = jet_pt[i];
      } else if (R < best_R) {
        best_R = R;
        best_jetpt = jet_pt[i];
      }
    } // jetloop
    return best_jetpt;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <bool IsMC, bool IsMix, typename TracksType, typename JetType>
  int minvReconstruction(double mult, const TracksType& trk1, const TracksType& trk2, const JetType& jets)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!trackSelection(trk1) || !trackSelection(trk2))
      return -1;

    if (!trackPID(trk1, true) || !trackPID(trk2, false))
      return -1;

    if (!cfgIsKstar) {
      if (trk1.globalIndex() >= trk2.globalIndex())
        return -1; // For Phi, we only need to iterate each pair once
    } else {
      if (trk1.globalIndex() == trk2.globalIndex())
        return -1; // For Kstar, we need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
    }

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    if (!cfgIsKstar)
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
    else
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);

    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgtrkMaxEta)
      return -1;

    /////////////////////////////////////////////////////////////////////////////
    // Fill Global Event Minv
    if (trk1.sign() * trk2.sign() < 0) {
      JEhistos.fill(HIST("hUSS_1D"), lResonance.M());
      if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
        JEhistos.fill(HIST("hUSS_1D_2_3"), lResonance.M());
      JEhistos.fill(HIST("hUSS"), mult, lResonance.Pt(), lResonance.M());

    } else if (trk1.sign() * trk2.sign() > 0) {

      JEhistos.fill(HIST("hLSS_1D"), lResonance.M());
      if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
        JEhistos.fill(HIST("hLSS_1D_2_3"), lResonance.M());
      JEhistos.fill(HIST("hLSS"), mult, lResonance.Pt(), lResonance.M());
    }
    /////////////////////////////////////////////////////////////////////////////

    bool jetFlag = false;
    int goodjets = 0;
    double jetpt = 0;
    for (auto const& jet : jets) {
      double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
      double etadiff = jet.eta() - lResonance.Eta();
      double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
      if (R < cfgjetR) {
        jetFlag = true;
        jetpt = jet.pt();
        goodjets++;
      }
    }

    if (cfgSingleJet)
      if (goodjets > 1)
        jetpt = DistinguishJets<JetType>(jets, lResonance);

    /////////////////////////////////////////////////////////////////////////////
    // Fill inside Jet
    if (jetFlag) {
      if (trk1.sign() * trk2.sign() < 0) {
        JEhistos.fill(HIST("hUSS_INSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hUSS_INSIDE_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hUSS_INSIDE"), jetpt, lResonance.Pt(), lResonance.M());

      } else if (trk1.sign() * trk2.sign() > 0) {

        JEhistos.fill(HIST("hLSS_INSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hLSS_INSIDE_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hLSS_INSIDE"), jetpt, lResonance.Pt(), lResonance.M());
      }
    } // jetflag
    /////////////////////////////////////////////////////////////////////////////

    if (!cfgIsKstar) {
      if (lResonance.M() > 1.005 && lResonance.M() < 1.035) {
        if (jetFlag)
          return 3;
        if (goodjets > 0)
          return 2;
        return 1;
      } else {
        return -1;
      }
    } else {
      if (lResonance.M() > 0.85 && lResonance.M() < 0.95) {
        if (jetFlag)
          return 3;
        if (goodjets > 0)
          return 2;
        return 1;
      } else {
        return -1;
      }
    }
  } // MinvReconstruction

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  int nEvents = 0;
  void processJetTracks(aod::JetCollision const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedjets, soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, TrackCandidates const&)
  {
    if (cDebugLevel > 0) {
      nEvents++;
      if ((nEvents + 1) % 10000 == 0) {
        std::cout << "Processed Data Events: " << nEvents << std::endl;
      }
    }
    JEhistos.fill(HIST("nEvents"), 0.5);

    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;

    int nReso = 0;
    int nResoWTrig = 0;
    int nResoInTrig = 0;
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>>();
      auto trk2 = track2.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>>();
      int Reso = minvReconstruction<false, false>(1.0, trk1, trk2, chargedjets);
      if (Reso > 0)
        nReso++;
      if (Reso > 1)
        nResoWTrig++;
      if (Reso > 2)
        nResoInTrig++;
    }
    // Here we have to fill the number of reso candidates per event
    JEhistos.fill(HIST("hNResoPerEvent"), nReso);
    JEhistos.fill(HIST("hNResoPerEventWJet"), nResoWTrig);
    JEhistos.fill(HIST("hNResoPerEventInJet"), nResoInTrig);

    int nJets = 0;
    for (auto chargedjet : chargedjets) {
      JEhistos.fill(HIST("FJetaHistogram"), chargedjet.eta());
      JEhistos.fill(HIST("FJphiHistogram"), chargedjet.phi());
      JEhistos.fill(HIST("FJptHistogram"), chargedjet.pt());
      nJets++;
    }

    JEhistos.fill(HIST("nJetsPerEvent"), nJets);

    JEhistos.fill(HIST("nEvents"), 1.5);
    //    return;

    for (auto& trackC : tracks) {
      auto originalTrack = trackC.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>>();
      JEhistos.fill(HIST("hDCArToPv"), originalTrack.dcaXY());
      JEhistos.fill(HIST("hDCAzToPv"), originalTrack.dcaZ());
      JEhistos.fill(HIST("rawpT"), originalTrack.pt());
      JEhistos.fill(HIST("rawDpT"), trackC.pt(), trackC.pt() - originalTrack.pt());
      JEhistos.fill(HIST("hIsPrim"), originalTrack.isPrimaryTrack());
      JEhistos.fill(HIST("hIsGood"), originalTrack.isGlobalTrackWoDCA());
      JEhistos.fill(HIST("hIsPrimCont"), originalTrack.isPVContributor());
      JEhistos.fill(HIST("hFindableTPCClusters"), originalTrack.tpcNClsFindable());
      JEhistos.fill(HIST("hFindableTPCRows"), originalTrack.tpcNClsCrossedRows());
      JEhistos.fill(HIST("hClustersVsRows"), originalTrack.tpcCrossedRowsOverFindableCls());
      JEhistos.fill(HIST("hTPCChi2"), originalTrack.tpcChi2NCl());
      JEhistos.fill(HIST("hITSChi2"), originalTrack.itsChi2NCl());

      if (!trackSelection(originalTrack))
        continue;

      JEhistos.fill(HIST("etaHistogram"), trackC.eta());
      JEhistos.fill(HIST("phiHistogram"), trackC.phi());
    } // JTrack Loop

  }; // Process Switch
  PROCESS_SWITCH(phiInJets, processJetTracks, "process JE Framework", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////MC STUFF////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>;
  using myCompleteJetTracks = soa::Join<aod::JetTracks, aod::JTrackPIs, aod::McTrackLabels>;
  int nJEEvents = 0;
  int nprocessRecEvents = 0;
  void processRec(o2::aod::JetCollision const& collision, myCompleteJetTracks const& tracks, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, aod::McParticles const&, myCompleteTracks const& /*originalTracks*/)
  {
    if (cDebugLevel > 0) {
      nprocessRecEvents++;
      if ((nprocessRecEvents + 1) % 10000 == 0) {
        double histmem = JEhistos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "processRec: " << nprocessRecEvents << std::endl;
      }
    }
    //=================
    // # of Events
    //=================
    JEhistos.fill(HIST("nEvents_MCRec"), 0.5);
    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (fabs(track.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    JEhistos.fill(HIST("nEvents_MCRec"), 1.5);

    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};

    bool hasJets = kFALSE;
    for (auto& mcdjet : mcdjets) {
      hasJets = kTRUE;
      mcd_pt.push_back(mcdjet.pt());
      mcd_eta.push_back(mcdjet.eta());
      mcd_phi.push_back(mcdjet.phi());
      JEhistos.fill(HIST("h_jet_pt"), mcdjet.pt());
      JEhistos.fill(HIST("h_jet_eta"), mcdjet.eta());
      JEhistos.fill(HIST("h_jet_phi"), mcdjet.phi());
    }
    if (hasJets)
      JEhistos.fill(HIST("nEvents_MCRec"), 2.5);

    double PhiCand = 0;
    double RealPhiCand = 0;
    double RealPhiCandWithJet = 0;
    double RealPhiCandInJet = 0;
    //============
    // Track Eff
    for (const auto& track : tracks) {
      auto originalTrack = track.track_as<myCompleteTracks>();
      if (!trackSelection(originalTrack))
        continue;
      if (cfgSimPID)
        if (!trackPID(originalTrack, true))
          continue;

      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();

        if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.eta()) <= cfgtrkMaxEta) {
          if (abs(mcParticle.pdgCode()) == 211)
            JEhistos.fill(HIST("ptJEHistogramPion"), mcParticle.pt());
          if (abs(mcParticle.pdgCode()) == 321)
            JEhistos.fill(HIST("ptJEHistogramKaon"), mcParticle.pt());
          if (abs(mcParticle.pdgCode()) == 2212)
            JEhistos.fill(HIST("ptJEHistogramProton"), mcParticle.pt());
        }
      }
      for (const auto& track2 : tracks) {
        auto originalTrack2 = track2.track_as<myCompleteTracks>();
        if (!trackSelection(originalTrack2))
          continue;
        if (cfgSimPID)
          if (!trackPID(originalTrack2, false))
            continue;
        if (!cfgIsKstar) {
          if (originalTrack.globalIndex() >= originalTrack2.globalIndex())
            continue;
        } else {
          if (originalTrack.globalIndex() == originalTrack2.globalIndex())
            continue;
        }
        if (fabs(originalTrack.eta()) > cfgtrkMaxEta || fabs(originalTrack2.eta()) > cfgtrkMaxEta)
          continue;

        TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
        lDecayDaughter1.SetXYZM(originalTrack.px(), originalTrack.py(), originalTrack.pz(), massKa);
        if (!cfgIsKstar)
          lDecayDaughter2.SetXYZM(originalTrack2.px(), originalTrack2.py(), originalTrack2.pz(), massKa);
        else
          lDecayDaughter2.SetXYZM(originalTrack2.px(), originalTrack2.py(), originalTrack2.pz(), massPi);
        lResonance = lDecayDaughter1 + lDecayDaughter2;

        if (fabs(lResonance.Eta()) > cfgtrkMaxEta)
          continue;

        if (lResonance.M() > 1.005 && lResonance.M() < 1.035)
          PhiCand++;

        //==================
        // 1.MB REC Closure
        //==================
        if (originalTrack.sign() * originalTrack2.sign() < 0) {
          JEhistos.fill(HIST("hMCRec_hUSS"), 1.0, lResonance.Pt(), lResonance.M());
        } else if (originalTrack.sign() * originalTrack2.sign() > 0) {
          JEhistos.fill(HIST("hMCRec_hLSS"), 1.0, lResonance.Pt(), lResonance.M());
        }
        //============================================
        // 2.Check if particle is inside a jet or not
        //============================================
        bool jetFlag = false;
        int goodjets = 0;
        double jetpt = 0;

        for (std::size_t i = 0; i < mcd_pt.size(); i++) {
          if (i == 0) {
            if (lResonance.M() > 1.005 && lResonance.M() < 1.035) {
              RealPhiCandWithJet++;
            }
          }
          double phidiff = TVector2::Phi_mpi_pi(mcd_phi[i] - lResonance.Phi());
          double etadiff = mcd_eta[i] - lResonance.Eta();
          double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));

          double phidiff_K1 = TVector2::Phi_mpi_pi(mcd_phi[i] - lDecayDaughter1.Phi());
          double etadiff_K1 = mcd_eta[i] - lDecayDaughter1.Eta();
          double R_K1 = TMath::Sqrt((etadiff_K1 * etadiff_K1) + (phidiff_K1 * phidiff_K1));

          double phidiff_K2 = TVector2::Phi_mpi_pi(mcd_phi[i] - lDecayDaughter2.Phi());
          double etadiff_K2 = mcd_eta[i] - lDecayDaughter2.Eta();
          double R_K2 = TMath::Sqrt((etadiff_K2 * etadiff_K2) + (phidiff_K2 * phidiff_K2));
          if (R < cfgjetR) {
            JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_Kangle_v_pt"), R_K1, lResonance.Pt());
            JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_Kangle_v_pt"), R_K2, lResonance.Pt());
          }
          if (R < cfgjetR) {
            jetFlag = true;
            jetpt = mcd_pt[i];
            goodjets++;
          }
        } // R check for jets

        //======================
        // 3.INSIDE REC Closure
        //======================
        if (jetFlag) {
          if (originalTrack.sign() * originalTrack2.sign() < 0) {
            JEhistos.fill(HIST("hMCRec_hUSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());
          } else if (originalTrack.sign() * originalTrack2.sign() > 0) {
            JEhistos.fill(HIST("hMCRec_hLSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());
          }
        }

        // check PID
        if (track.has_mcParticle() && track2.has_mcParticle()) {
          auto part1 = track.mcParticle();
          auto part2 = track2.mcParticle();
          if (fabs(part1.pdgCode()) != 321)
            continue; // Not Kaon
          if (!cfgIsKstar) {
            if (fabs(part2.pdgCode()) != 321)
              continue; // Not Kaon
          } else {
            if (fabs(part2.pdgCode()) != 211)
              continue; // Not Kaon
          }

          if (!part1.has_mothers())
            continue; // Not decaying Kaon
          if (!part2.has_mothers())
            continue; // Not decaying Kaon

          std::vector<int> mothers1{};
          std::vector<int> mothers1PDG{};
          for (auto& part1_mom : part1.mothers_as<aod::McParticles>()) {
            mothers1.push_back(part1_mom.globalIndex());
            mothers1PDG.push_back(part1_mom.pdgCode());
          }

          std::vector<int> mothers2{};
          std::vector<int> mothers2PDG{};
          for (auto& part2_mom : part2.mothers_as<aod::McParticles>()) {
            mothers2.push_back(part2_mom.globalIndex());
            mothers2PDG.push_back(part2_mom.pdgCode());
          }

          if (!cfgIsKstar) {
            if (mothers1PDG[0] != 333)
              continue; // mother not phi
            if (mothers2PDG[0] != 333)
              continue; // mother not phi
          } else {
            if (mothers1PDG[0] != 313)
              continue; // mother not phi
            if (mothers2PDG[0] != 313)
              continue; // mother not phi
          }
          if (mothers1[0] != mothers2[0])
            continue; // Kaons not from the same phi

          double phidiff_Kaons = TVector2::Phi_mpi_pi(lDecayDaughter2.Phi() - lDecayDaughter1.Phi());
          double etadiff_Kaons = lDecayDaughter2.Eta() - lDecayDaughter1.Eta();
          double R_Kaons = TMath::Sqrt((etadiff_Kaons * etadiff_Kaons) + (phidiff_Kaons * phidiff_Kaons));
          JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_KtoKangle_v_pt"), R_Kaons, lResonance.Pt());
          JEhistos.fill(HIST("ptJEHistogramPhi"), lResonance.Pt());

          //=====================
          // 4.MB True Closure
          //=====================
          if (originalTrack.sign() * originalTrack2.sign() < 0) {
            JEhistos.fill(HIST("hMCRecTrue_hUSS"), 1.0, lResonance.Pt(), lResonance.M());
          } else if (originalTrack.sign() * originalTrack2.sign() > 0) {
            JEhistos.fill(HIST("hMCRecTrue_hLSS"), 1.0, lResonance.Pt(), lResonance.M());
          }

          //===========================
          // 5.INSIDE REC True Closure
          //===========================
          if (jetFlag) {
            if (originalTrack.sign() * originalTrack2.sign() < 0) {
              JEhistos.fill(HIST("hMCRecTrue_hUSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());
            } else if (originalTrack.sign() * originalTrack2.sign() > 0) {
              JEhistos.fill(HIST("hMCRecTrue_hLSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());
            }
          }

          if (lResonance.M() > 1.005 && lResonance.M() < 1.035)
            RealPhiCand++;

          // Now we do jets
          if (cfgSingleJet)
            if (goodjets > 1)
              jetpt = DistinguishJetsMC(mcd_pt, mcd_phi, mcd_eta, lResonance);

          if (jetFlag) {
            if (lResonance.M() > 1.005 && lResonance.M() < 1.035)
              RealPhiCandInJet++;
            JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_INSIDE_pt_v_eta"), lResonance.Pt(), lResonance.Eta());
            JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_INSIDE_1D"), lResonance.M());
            if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
              JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_INSIDE_1D_2_3"), lResonance.M());
            JEhistos.fill(HIST("hMCRec_nonmatch_hUSS_INSIDE"), jetpt, lResonance.Pt(), lResonance.M());
          }

          if (hasJets) {
            JEhistos.fill(HIST("ptJEHistogramPhi_JetTrigger"), lResonance.Pt());
            auto triggerjet = std::min_element(mcd_pt.begin(), mcd_pt.end());
            double triggerjet_pt = *triggerjet;
            JEhistos.fill(HIST("JetVsPhi_REC"), triggerjet_pt, lResonance.Pt());
          }
          JEhistos.fill(HIST("minvJEHistogramPhi"), lResonance.M());
        } // mcpart check
      } // tracks2
    } // tracks1
    JEhistos.fill(HIST("hNRealPhiVPhiCand"), PhiCand, RealPhiCand);
    JEhistos.fill(HIST("hNRealPhiWithJetVPhiCand"), PhiCand, RealPhiCandWithJet);
    JEhistos.fill(HIST("hNRealPhiInJetVPhiCand"), PhiCand, RealPhiCandInJet);

    // Jet Eff
  }
  PROCESS_SWITCH(phiInJets, processRec, "pikp detector level MC JE", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int nprocessSimEvents = 0;
  //  Preslice<aod::JCollisions> slice = o2::aod::JCollision::collisionId;
  void processSim(o2::aod::JetMcCollision const& collision, soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls, aod::JetParticles const& mcParticles, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
    if (cDebugLevel > 0) {
      nprocessSimEvents++;
      if ((nprocessSimEvents + 1) % 10000 == 0) {
        double histmem = JEhistos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "processSim: " << nprocessSimEvents << std::endl;
      }
    }

    JEhistos.fill(HIST("nEvents_MCGen"), 0.5);

    if (recocolls.size() <= 0) // not reconstructed
      return;

    for (auto& recocoll : recocolls) { // poorly reconstructed
      if (!jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits))
        return;
    }

    if (fabs(collision.posZ()) > cfgVtxCut) // bad vertex
      return;
    bool INELgt0 = false;
    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }

    if (!INELgt0) // not INEL
      return;
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};
    JEhistos.fill(HIST("nEvents_MCGen"), 1.5);

    bool hasJets = kFALSE;
    // Check jets
    for (auto& mcpjet : mcpjets) {
      hasJets = kTRUE;
      mcp_pt.push_back(mcpjet.pt());
      mcp_eta.push_back(mcpjet.eta());
      mcp_phi.push_back(mcpjet.phi());

      JEhistos.fill(HIST("h_part_jet_pt"), mcpjet.pt());
      JEhistos.fill(HIST("h_part_jet_eta"), mcpjet.eta());
      JEhistos.fill(HIST("h_part_jet_phi"), mcpjet.phi());
    }
    if (hasJets)
      JEhistos.fill(HIST("nEvents_MCGen"), 2.5);

    // Check pikp and phi
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.eta()) <= cfgtrkMaxEta) { // watch out for context!!!
        if (abs(mcParticle.pdgCode()) == 211)
          JEhistos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
        if (abs(mcParticle.pdgCode()) == 321)
          JEhistos.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
        if (abs(mcParticle.pdgCode()) == 2212)
          JEhistos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
      }
      if (fabs(mcParticle.eta()) <= cfgtrkMaxEta) { // watch out for context!!!
        TLorentzVector lResonance;
        lResonance.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());

        int GenPID = 0;
        if (!cfgIsKstar)
          GenPID = 333;
        else
          GenPID = 313;

        if (abs(mcParticle.pdgCode()) == GenPID) {
          JEhistos.fill(HIST("ptGeneratedPhi_ALLBR"), mcParticle.pt());

          bool skip = false;

          // First we check for Forced BR
          // if we check for Phi
          if (!cfgIsKstar) {
            if (mcParticle.has_daughters())
              for (auto& dgth : mcParticle.daughters_as<aod::JetParticles>())
                if (fabs(dgth.pdgCode()) != 321)
                  skip = true;
          } else {
            if (mcParticle.has_daughters())
              for (auto& dgth : mcParticle.daughters_as<aod::JetParticles>())
                if (fabs(dgth.pdgCode()) != 321 || fabs(dgth.pdgCode()) != 211)
                  skip = true;
          }

          if (skip && cfgBR)
            continue;
          JEhistos.fill(HIST("ptGeneratedPhi"), mcParticle.pt());
          JEhistos.fill(HIST("mGeneratedPhi"), lResonance.M());

          ////////////////////////////Implementation of phi finding

          int goodjets = 0;
          double jetpt = 0;
          TLorentzVector lResonance;
          lResonance.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
          bool jetFlag = false;
          for (std::vector<double>::size_type i = 0; i < mcp_pt.size(); i++) {
            double phidiff = TVector2::Phi_mpi_pi(mcp_phi[i] - lResonance.Phi());
            double etadiff = mcp_eta[i] - lResonance.Eta();
            double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
            if (mcParticle.has_daughters()) {
              for (auto& dgth : mcParticle.daughters_as<aod::JetParticles>()) {
                double phidiff_K = TVector2::Phi_mpi_pi(mcp_phi[i] - dgth.phi());
                double etadiff_K = mcp_eta[i] - dgth.eta();
                double R_K = TMath::Sqrt((etadiff_K * etadiff_K) + (phidiff_K * phidiff_K));
                if (R < cfgjetR)
                  JEhistos.fill(HIST("hMCTrue_nonmatch_hUSS_Kangle_v_pt"), R_K, lResonance.Pt());
              }
            }
            if (R < cfgjetR) {
              jetFlag = true;
              jetpt = mcp_pt[i];
              goodjets++;
            }
          }
          if (cfgSingleJet)
            if (goodjets > 1)
              jetpt = DistinguishJetsMC(mcp_pt, mcp_phi, mcp_eta, lResonance);

          if (jetFlag) {
            JEhistos.fill(HIST("hMCTrue_nonmatch_hUSS_INSIDE_pt_v_eta"), lResonance.Pt(), lResonance.Eta());
            JEhistos.fill(HIST("hMCTrue_nonmatch_hUSS_INSIDE_1D"), lResonance.M());
            if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
              JEhistos.fill(HIST("hMCTrue_nonmatch_hUSS_INSIDE_1D_2_3"), lResonance.M());
            JEhistos.fill(HIST("hMCTrue_nonmatch_hUSS_INSIDE"), jetpt, lResonance.Pt(), lResonance.M());
          }

          ////////////////////////////Phi found
          if (hasJets) {
            JEhistos.fill(HIST("ptGeneratedPhi_JetTrigger"), mcParticle.pt());
            auto triggerjet = std::min_element(mcp_pt.begin(), mcp_pt.end());
            double triggerjet_pt = *triggerjet;
            JEhistos.fill(HIST("JetVsPhi_GEN"), triggerjet_pt, mcParticle.pt());
          } // check for jets

        } // check for phi
      } // check for rapidity
    } // loop over particles
  } // process switch
  PROCESS_SWITCH(phiInJets, processSim, "pikp particle level MC", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  using JetMCDTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;

  int nprocessSimJEEvents = 0;
  void processMatchedGen(aod::JetMcCollision const& collision,
                         soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls,
                         JetMCDTable const& /*mcdjets*/,
                         JetMCPTable const& mcpjets,
                         myCompleteJetTracks const& tracks,
                         myCompleteTracks const&,
                         aod::JetParticles const& mcParticles,
                         aod::McParticles const&)

  {
    if (cDebugLevel > 0) {
      nprocessSimJEEvents++;
      if ((nprocessSimJEEvents + 1) % 10000 == 0) {
        double histmem = JEhistos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "processMatchedGen Events: " << nprocessSimJEEvents << std::endl;
      }
    }

    JEhistos.fill(HIST("nEvents_MCGen_MATCHED"), 0.5);

    if (fabs(collision.posZ()) > cfgVtxCut)
      return;

    if (recocolls.size() <= 0) // not reconstructed
      return;
    for (auto& recocoll : recocolls) { // poorly reconstructed
      if (!jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits))
        return;
    }

    bool INELgt0 = false;
    for (const auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    JEhistos.fill(HIST("nEvents_MCGen_MATCHED"), 1.5);

    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};

    bool hasJets = kFALSE;

    for (auto& mcpjet : mcpjets) {
      if (!mcpjet.has_matchedJetGeo())
        continue;

      for (auto& mcdjet : mcpjet.template matchedJetGeo_as<JetMCDTable>()) {
        if (!mcdjet.has_matchedJetGeo())
          continue;
        hasJets = kTRUE;
        JEhistos.fill(HIST("h_matched_GEN_jet_pt"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        JEhistos.fill(HIST("h_matched_GEN_jet_phi"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        JEhistos.fill(HIST("h_matched_GEN_jet_eta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        mcd_pt.push_back(mcdjet.pt());
        mcd_eta.push_back(mcdjet.eta());
        mcd_phi.push_back(mcdjet.phi());
        mcp_pt.push_back(mcpjet.pt());
        mcp_eta.push_back(mcpjet.eta());
        mcp_phi.push_back(mcpjet.phi());
      } // mcpjets
    } // mcdjets

    if (hasJets)
      JEhistos.fill(HIST("nEvents_MCGen_MATCHED"), 2.5);

    // First we do GEN part
    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.eta()) > cfgtrkMaxEta)
        continue;

      int GenPID = 0;

      if (!cfgIsKstar)
        GenPID = 333;
      else
        GenPID = 313;

      if (fabs(mcParticle.pdgCode()) == GenPID) {
        bool skip = false;
        double phi_dgth_px[2] = {0};
        double phi_dgth_py[2] = {0};
        double phi_dgth_pz[2] = {0};
        double TEMP_phi_dgth_px[2] = {0};
        double TEMP_phi_dgth_py[2] = {0};
        double TEMP_phi_dgth_pz[2] = {0};

        bool good_daughter[2] = {false};
        int dgth_index = 0;
        // First we check for Forced BR
        // if we check for Phi
        if (!cfgIsKstar) {
          if (mcParticle.has_daughters()) {
            for (auto& dgth : mcParticle.daughters_as<aod::JetParticles>()) {
              if (fabs(dgth.pdgCode()) != 321) {
                skip = true;
                break;
              }
              for (const auto& track : tracks) {
                auto trk = track.track_as<myCompleteTracks>();
                if (!trackSelection(trk))
                  continue;
                if (fabs(trk.eta()) > cfgtrkMaxEta)
                  continue;
                if (cfgSimPID) {
                  if (!trackPID(trk, true))
                    continue;
                }
                if (!trk.has_mcParticle())
                  continue;
                auto part = trk.mcParticle();
                if (part.globalIndex() == dgth.globalIndex()) {
                  TEMP_phi_dgth_px[dgth_index] = track.px();
                  TEMP_phi_dgth_py[dgth_index] = track.py();
                  TEMP_phi_dgth_pz[dgth_index] = track.pz();
                  good_daughter[dgth_index] = true;
                  dgth_index++;
                  if (dgth_index == 2) {
                    phi_dgth_px[0] = TEMP_phi_dgth_px[0];
                    phi_dgth_py[0] = TEMP_phi_dgth_py[0];
                    phi_dgth_pz[0] = TEMP_phi_dgth_pz[0];
                    phi_dgth_px[1] = TEMP_phi_dgth_px[1];
                    phi_dgth_py[1] = TEMP_phi_dgth_py[1];
                    phi_dgth_pz[1] = TEMP_phi_dgth_pz[1];
                    break;
                  }
                } // index check
              } // track loop
            } // mc daughter loop
          } // check if particle has daughters
        } else { // check for kstar
          if (mcParticle.has_daughters())
            for (auto& dgth : mcParticle.daughters_as<aod::JetParticles>())
              if (fabs(dgth.pdgCode()) != 321 || fabs(dgth.pdgCode()) != 211)
                skip = true;
        }

        if (skip && cfgBR)
          continue;

        int goodjets = 0;
        double jetpt_mcd = 0;
        double jetpt_mcp = 0;
        TLorentzVector lResonance;
        TLorentzVector lResonance_REC;
        TLorentzVector lDecayDaughter1_REC;
        TLorentzVector lDecayDaughter2_REC;
        lResonance.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
        lDecayDaughter1_REC.SetXYZM(phi_dgth_px[0], phi_dgth_py[0], phi_dgth_pz[0], massKa);
        lDecayDaughter2_REC.SetXYZM(phi_dgth_px[1], phi_dgth_py[1], phi_dgth_pz[1], massKa);
        lResonance_REC = lDecayDaughter1_REC + lDecayDaughter2_REC;

        bool jetFlag = false;
        for (std::vector<double>::size_type i = 0; i < mcp_pt.size(); i++) {
          double phidiff = TVector2::Phi_mpi_pi(mcp_phi[i] - lResonance.Phi());
          double etadiff = mcp_eta[i] - lResonance.Eta();
          double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
          if (R < cfgjetR) {
            jetFlag = true;
            jetpt_mcp = mcp_pt[i];
            jetpt_mcd = mcd_pt[i];
            goodjets++;
          }
        }
        if (cfgSingleJet) {
          if (goodjets > 1) {
            jetpt_mcp = DistinguishJetsMC(mcp_pt, mcp_phi, mcp_eta, lResonance);
            jetpt_mcd = DistinguishJetsMC(mcd_pt, mcd_phi, mcd_eta, lResonance_REC);
          }
        }

        if (jetFlag) {
          if (cDebugLevel > 0) {
            std::cout << "******************************************" << std::endl;
            std::cout << "GEN TO REC LEVEL: " << std::endl;
            std::cout << "Rec. Phi Pt: " << lResonance_REC.Pt() << std::endl;
            std::cout << "Rec. Phi Eta: " << lResonance_REC.Eta() << std::endl;
            std::cout << "Rec. Jet Pt: " << jetpt_mcd << std::endl;
            std::cout << "Gen. Phi Pt: " << lResonance.Pt() << std::endl;
            std::cout << "Gen. Phi Eta: " << lResonance.Eta() << std::endl;
            std::cout << "Gen. Jet Pt: " << jetpt_mcp << std::endl;
            std::cout << "******************************************" << std::endl;
          }

          JEhistos.fill(HIST("RespGen_Matrix_MATCHED"), lResonance_REC.Pt(), jetpt_mcd, lResonance.Pt(), jetpt_mcp);
          unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
          int dice = rand_r(&seed) % 2;
          if (dice > 0)
            JEhistos.fill(HIST("RespGen_Matrix_MATCHED_rand0"), lResonance_REC.Pt(), jetpt_mcd, lResonance.Pt(), jetpt_mcp);
          else
            JEhistos.fill(HIST("RespGen_Matrix_MATCHED_rand1"), lResonance_REC.Pt(), jetpt_mcd, lResonance.Pt(), jetpt_mcp);

          JEhistos.fill(HIST("2DGenToRec"), lResonance.Pt(), jetpt_mcp);
          // //check constrained eff
          if (good_daughter[0] && good_daughter[1] && lResonance_REC.Pt() > 0 && lResonance_REC.Pt() < 20.0 && jetpt_mcd > 8 && jetpt_mcd < 200)
            JEhistos.fill(HIST("2DGenToRec_constrained"), lResonance.Pt(), jetpt_mcp);

          JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE"), jetpt_mcp, lResonance.Pt(), lResonance.M());
        }
      } // chech for phi
    } // MC Particles
  } // main fcn
  PROCESS_SWITCH(phiInJets, processMatchedGen, "phi matched level MC", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int nprocessRecJEEvents = 0;
  //  void processMatchedRec(o2::aod::JCollision const& collision, myCompleteJetTracks const& tracks, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, aod::McParticles const&, myCompleteTracks const& originalTracks)
  void processMatchedRec(aod::JetCollision const& collision,
                         JetMCDTable const& mcdjets,
                         JetMCPTable const&,
                         myCompleteJetTracks const& tracks,
                         myCompleteTracks const&,
                         aod::McParticles const&)
  {

    if (cDebugLevel > 0) {
      nprocessRecJEEvents++;
      if ((nprocessRecJEEvents + 1) % 10000 == 0) {
        double histmem = JEhistos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "processMatched Rec Events: " << nprocessRecJEEvents << std::endl;
      }
    }
    JEhistos.fill(HIST("nEvents_MCRec_MATCHED"), 0.5);

    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (fabs(track.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    JEhistos.fill(HIST("nEvents_MCRec_MATCHED"), 1.5);

    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};

    bool hasJets = kFALSE;
    for (auto& mcdjet : mcdjets) {
      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {
        if (!mcpjet.has_matchedJetGeo())
          continue;
        hasJets = kTRUE;
        JEhistos.fill(HIST("h_matched_REC_jet_pt"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        JEhistos.fill(HIST("h_matched_REC_jet_phi"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        JEhistos.fill(HIST("h_matched_REC_jet_eta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        mcd_pt.push_back(mcdjet.pt());
        mcd_eta.push_back(mcdjet.eta());
        mcd_phi.push_back(mcdjet.phi());
        mcp_pt.push_back(mcpjet.pt());
        mcp_eta.push_back(mcpjet.eta());
        mcp_phi.push_back(mcpjet.phi());
      } // mcpjets
    } // mcdjets
    // Now we do REC part
    if (hasJets)
      JEhistos.fill(HIST("nEvents_MCRec_MATCHED"), 2.5);

    for (const auto& track1 : tracks) {
      auto trk1 = track1.track_as<myCompleteTracks>();
      for (const auto& track2 : tracks) {
        auto trk2 = track2.track_as<myCompleteTracks>();
        if (!cfgIsKstar) {
          if (trk1.globalIndex() >= trk2.globalIndex())
            continue;
        } else {
          if (trk1.globalIndex() == trk2.globalIndex())
            continue;
        }
        if (fabs(trk1.eta()) > cfgtrkMaxEta || fabs(trk2.eta()) > cfgtrkMaxEta)
          continue;
        if ((trk1.sign() * trk2.sign()) > 0)
          continue; // Not K+K-
        if (trackSelection(trk1) && trackSelection(trk2)) {
          if (cfgSimPID) {
            if (!trackPID(trk1, true))
              continue;
            if (!trackPID(trk2, false))
              continue;
          }
          if (track1.has_mcParticle() && track2.has_mcParticle()) {
            auto part1 = track1.mcParticle();
            auto part2 = track2.mcParticle();
            if (fabs(part1.pdgCode()) != 321)
              continue; // Not Kaon
            if (!cfgIsKstar) {
              if (fabs(part2.pdgCode()) != 321)
                continue; // Not Kaon
            } else {
              if (fabs(part2.pdgCode()) != 211)
                continue; // Not Kaon
            }
            if (!part1.has_mothers())
              continue; // Not decaying Kaon
            if (!part2.has_mothers())
              continue; // Not decaying Kaon

            std::vector<int> mothers1{};
            std::vector<int> mothers1PDG{};
            std::vector<float> mothers1Pt{};
            for (auto& part1_mom : part1.mothers_as<aod::McParticles>()) {
              mothers1.push_back(part1_mom.globalIndex());
              mothers1PDG.push_back(part1_mom.pdgCode());
              mothers1Pt.push_back(part1_mom.pt());
            }

            std::vector<int> mothers2{};
            std::vector<int> mothers2PDG{};
            std::vector<float> mothers2Pt{};
            for (auto& part2_mom : part2.mothers_as<aod::McParticles>()) {
              mothers2.push_back(part2_mom.globalIndex());
              mothers2PDG.push_back(part2_mom.pdgCode());
              mothers2Pt.push_back(part2_mom.pt());
            }
            if (!cfgIsKstar) {
              if (mothers1PDG[0] != 333)
                continue; // mother not phi
              if (mothers2PDG[0] != 333)
                continue; // mother not phi
            } else {
              if (mothers1PDG[0] != 313)
                continue; // mother not phi
              if (mothers2PDG[0] != 313)
                continue; // mother not phi
            }
            if (mothers1[0] != mothers2[0])
              continue; // Kaons not from the same phi

            TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
            lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
            if (!cfgIsKstar)
              lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
            else
              lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);

            lResonance = lDecayDaughter1 + lDecayDaughter2;
            if (fabs(lResonance.Eta()) > cfgtrkMaxEta)
              continue;

            bool jetFlag = false;
            int goodjets = 0;
            double jetpt_mcp = 0;
            double jetpt_mcd = 0;
            for (std::vector<double>::size_type i = 0; i < mcd_pt.size(); i++) {
              double phidiff = TVector2::Phi_mpi_pi(mcd_phi[i] - lResonance.Phi());
              double etadiff = mcd_eta[i] - lResonance.Eta();
              double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
              if (R < cfgjetR) {
                jetpt_mcd = mcd_pt[i];
                jetpt_mcp = mcp_pt[i];
                goodjets++;
                jetFlag = true;
              }
            }

            if (cfgSingleJet) {
              if (goodjets > 1) {
                jetpt_mcd = DistinguishJetsMC(mcd_pt, mcd_phi, mcd_eta, lResonance);
                jetpt_mcp = DistinguishJetsMC(mcp_pt, mcp_phi, mcp_eta, lResonance);
              }
            }
            if (jetFlag) { // Fill Resp. Matrix
              if (cDebugLevel > 0) {
                std::cout << "******************************************" << std::endl;
                std::cout << "REC TO GEN LEVEL: " << std::endl;
                std::cout << "Rec. Phi Pt: " << lResonance.Pt() << std::endl;
                std::cout << "Rec. Jet Pt: " << jetpt_mcd << std::endl;
                std::cout << "Gen. Phi Pt: " << mothers1Pt[0] << std::endl;
                std::cout << "Gen. Jet Pt: " << jetpt_mcp << std::endl;
                std::cout << "******************************************" << std::endl;
              }
              JEhistos.fill(HIST("Resp_Matrix_MATCHED"), lResonance.Pt(), jetpt_mcd, mothers1Pt[0], jetpt_mcp);
              unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
              int dice = rand_r(&seed) % 2;
              if (dice > 0)
                JEhistos.fill(HIST("Resp_Matrix_MATCHED_rand0"), lResonance.Pt(), jetpt_mcd, mothers1Pt[0], jetpt_mcp);
              else
                JEhistos.fill(HIST("Resp_Matrix_MATCHED_rand1"), lResonance.Pt(), jetpt_mcd, mothers1Pt[0], jetpt_mcp);

              JEhistos.fill(HIST("2DRecToGen"), lResonance.Pt(), jetpt_mcd);
              // check constrained eff
              if (mothers1Pt[0] > 0 && mothers1Pt[0] < 20.0 && jetpt_mcp > 8 && jetpt_mcp < 200)
                JEhistos.fill(HIST("2DRecToGen_constrained"), lResonance.Pt(), jetpt_mcd);
            }

            // Fill 3D Invariant mass distributions
            if (jetFlag) {
              JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D"), lResonance.M());
              if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
                JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D_2_3"), lResonance.M());
              JEhistos.fill(HIST("hMCRec_hUSS_INSIDE"), jetpt_mcd, lResonance.Pt(), lResonance.M());
            }
            // else if (!jetFlag && mcd_pt.size() > 0) {
            //    JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_TRIG_1D"), lResonance.M());

            //    if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            //      JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_TRIG_1D_2_3"), lResonance.M());

            //    JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_TRIG"), jetpt_mcd, lResonance.Pt(), lResonance.M());

            //  } else if (!jetFlag) {
            //    JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D"), lResonance.M());

            //    if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            //      JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D_2_3"), lResonance.M());

            //    JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE"), jetpt_mcd, lResonance.Pt(), lResonance.M());

            //         }  //! jetflag

          } // pass track cut
        } // has mc particle

      } // tracks 2
    } // tracks 1
    // }   // tracks
  } // main fcn
  PROCESS_SWITCH(phiInJets, processMatchedRec, "phi matched Rec level MC", true);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiInJets>(cfgc)};
};
