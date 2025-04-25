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
// This is a task that reads sigma0 tables (from sigma0builder) to perform analysis.
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Sigma0 analysis task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "CCDB/BasicCCDBManager.h"
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using std::array;
using V0MCSigmas = soa::Join<aod::Sigma0Cores, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras, aod::SigmaMCCores>;
using V0Sigmas = soa::Join<aod::Sigma0Cores, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras>;

static const std::vector<std::string> PhotonSels = {"NoSel", "V0Type", "DaupT", "DCADauToPV",
                                                    "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                    "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                    "PsiPair", "Phi", "Mass"};

static const std::vector<std::string> LambdaSels = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                    "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                    "TPCTOFPID", "DCADauToPV", "Mass"};

static const std::vector<std::string> DirList = {"BeforeSel", "AfterSel"};

struct sigmaanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> fillQAhistos{"fillQAhistos", false, "if true, fill QA histograms"};
  Configurable<bool> fillBkgQAhistos{"fillBkgQAhistos", false, "if true, fill MC QA histograms for Bkg study. Only works with MC."};
  Configurable<bool> fillpTResoQAhistos{"fillpTResoQAhistos", false, "if true, fill MC QA histograms for pT resolution study. Only works with MC."};

  // Interaction rate selection:
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<float> minIR{"minIR", -1, "Min Interaction Rate (kHz). Leave -1 if no selection desired."};
  Configurable<float> maxIR{"maxIR", -1, "Max Interaction Rate (kHz). Leave -1 if no selection desired."};

  // Analysis strategy:
  Configurable<bool> fUseMLSel{"fUseMLSel", false, "Flag to use ML selection. If False, the standard selection is applied."};
  Configurable<bool> fselLambdaTPCPID{"fselLambdaTPCPID", true, "Flag to select lambda-like candidates using TPC NSigma."};
  Configurable<bool> fselLambdaTOFPID{"fselLambdaTOFPID", false, "Flag to select lambda-like candidates using TOF NSigma."};
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "Flag to process only signal candidates. Use only with processMonteCarlo!"};
  Configurable<bool> doPhotonLambdaSelQA{"doPhotonLambdaSelQA", false, "Flag to fill photon and lambda QA histos!"};

  // For ML Selection
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For Standard Selection:
  //// Lambda standard criteria::
  Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> ALambdaMinDCANegToPv{"ALambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> ALambdaMinDCAPosToPv{"ALambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
  Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 40, "Max V0 radius (cm)"};
  Configurable<float> LambdaMinQt{"LambdaMinQt", 0.01, "Min lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMaxQt{"LambdaMaxQt", 0.17, "Max lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMinAlpha{"LambdaMinAlpha", 0.25, "Min lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMaxAlpha{"LambdaMaxAlpha", 1.0, "Max lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMinv0cospa{"LambdaMinv0cospa", 0.95, "Min V0 CosPA"};
  Configurable<float> LambdaMaxLifeTime{"LambdaMaxLifeTime", 30, "Max lifetime"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.015, "Mass window around expected (in GeV/c2)"};
  Configurable<float> LambdaMaxRap{"LambdaMaxRap", 0.8, "Max lambda rapidity"};
  Configurable<float> LambdaMaxDauEta{"LambdaMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
  Configurable<float> LambdaMaxTPCNSigmas{"LambdaMaxTPCNSigmas", 1e+9, "Max TPC NSigmas for daughters"};
  Configurable<float> LambdaPrMaxTOFNSigmas{"LambdaPrMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
  Configurable<float> LambdaPiMaxTOFNSigmas{"LambdaPiMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
  Configurable<int> LambdaMinTPCCrossedRows{"LambdaMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
  Configurable<int> LambdaMinITSclusters{"LambdaMinITSclusters", 1, "minimum ITS clusters"};
  Configurable<bool> LambdaRejectPosITSafterburner{"LambdaRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
  Configurable<bool> LambdaRejectNegITSafterburner{"LambdaRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};

  //// Photon standard criteria:
  Configurable<int> Photonv0TypeSel{"Photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
  Configurable<float> PhotonDauMinPt{"PhotonDauMinPt", 0.0, "Min daughter pT (GeV/c)"};
  Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<int> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
  Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
  Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
  Configurable<float> PhotonMinPt{"PhotonMinPt", 0.0, "Min photon pT (GeV/c)"};
  Configurable<float> PhotonMaxPt{"PhotonMaxPt", 50.0, "Max photon pT (GeV/c)"};
  Configurable<float> PhotonMaxRap{"PhotonMaxRap", 0.5, "Max photon rapidity"};
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
  Configurable<float> PhotonMaxZ{"PhotonMaxZ", 240, "Max photon conversion point z value (cm)"};
  Configurable<float> PhotonMaxQt{"PhotonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
  Configurable<float> PhotonMaxAlpha{"PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
  Configurable<float> PhotonMinV0cospa{"PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
  Configurable<float> PhotonPsiPairMax{"PhotonPsiPairMax", 1e+9, "maximum psi angle of the track pair"};
  Configurable<float> PhotonMaxDauEta{"PhotonMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonLineCutZ0{"PhotonLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
  Configurable<float> PhotonPhiMin1{"PhotonPhiMin1", -1, "Phi min value to reject photons, region 1 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMax1{"PhotonPhiMax1", -1, "Phi max value to reject photons, region 1 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMin2{"PhotonPhiMin2", -1, "Phi max value to reject photons, region 2 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMax2{"PhotonPhiMax2", -1, "Phi min value to reject photons, region 2 (leave negative if no selection desired)"};

  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.5, "Max sigma0 rapidity"};

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisInvPt{"axisInvPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {400, -50.0, 50.0}, ""};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {150, 0, 1500}, "Binning for the interaction rate (kHz)"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality, PID and other axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisNCls{"axisNCls", {8, -0.5, 7.5}, "NCls"};
  ConfigurableAxis axisChi2PerNcl{"axisChi2PerNcl", {80, -40, 40}, "Chi2 Per Ncl"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {120, -30, 30}, "TPC NSigma"};
  ConfigurableAxis axisTOFNSigma{"axisTOFNSigma", {120, -30, 30}, "TOF NSigma"};
  ConfigurableAxis axisLifetime{"axisLifetime", {100, 0, 100}, "Chi2 Per Ncl"};

  // topological variable QA axes
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.5f, 1.0f}, "Cosine of pointing angle"};
  ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  ConfigurableAxis axisPsiPair{"axisPsiPair", {250, -5.0f, 5.0f}, "Psipair for photons"};
  ConfigurableAxis axisPhi{"axisPhi", {200, 0, 2 * o2::constants::math::PI}, "Phi for photons"};
  ConfigurableAxis axisZ{"axisZ", {120, -120.0f, 120.0f}, "V0 Z position (cm)"};

  ConfigurableAxis axisCandSel{"axisCandSel", {20, 0.5f, +20.5f}, "Candidate Selection"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};

  void init(InitContext const&)
  {

    for (const auto& histodir : DirList) {

      histos.add(histodir + "/Photon/hTrackCode", "hTrackCode", kTH1F, {{11, 0.5f, 11.5f}});
      histos.add(histodir + "/Photon/hV0Type", "hV0Type", kTH1F, {{8, 0.5f, 8.5f}});
      histos.add(histodir + "/Photon/hNegpT", "hNegpT", kTH1F, {axisPt});
      histos.add(histodir + "/Photon/hPospT", "hPospT", kTH1F, {axisPt});
      histos.add(histodir + "/Photon/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Photon/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Photon/hDCADau", "hDCADau", kTH1F, {axisDCAdau});
      histos.add(histodir + "/Photon/hPosTPCCR", "hPosTPCCR", kTH1F, {axisTPCrows});
      histos.add(histodir + "/Photon/hNegTPCCR", "hNegTPCCR", kTH1F, {axisTPCrows});
      histos.add(histodir + "/Photon/h2dPosTPCNSigmaEl", "h2dPosTPCNSigmaEl", kTH2F, {axisPt, axisTPCNSigma});
      histos.add(histodir + "/Photon/h2dNegTPCNSigmaEl", "h2dNegTPCNSigmaEl", kTH2F, {axisPt, axisTPCNSigma});
      histos.add(histodir + "/Photon/h2dPosTPCNSigmaPi", "h2dPosTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
      histos.add(histodir + "/Photon/h2dNegTPCNSigmaPi", "h2dNegTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
      histos.add(histodir + "/Photon/hpT", "hpT", kTH1F, {axisPt});
      histos.add(histodir + "/Photon/hY", "hY", kTH1F, {axisRapidity});
      histos.add(histodir + "/Photon/hPosEta", "hPosEta", kTH1F, {axisRapidity});
      histos.add(histodir + "/Photon/hNegEta", "hNegEta", kTH1F, {axisRapidity});
      histos.add(histodir + "/Photon/hRadius", "hRadius", kTH1F, {axisRadius});
      histos.add(histodir + "/Photon/hZ", "hZ", kTH1F, {axisZ});
      histos.add(histodir + "/Photon/h2dRZCut", "h2dRZCut", kTH2F, {axisZ, axisRadius});
      histos.add(histodir + "/Photon/h2dRZPlane", "h2dRZPlane", kTH2F, {axisZ, axisRadius});
      histos.add(histodir + "/Photon/hCosPA", "hCosPA", kTH1F, {axisCosPA});
      histos.add(histodir + "/Photon/hPsiPair", "hPsiPair", kTH1F, {axisPsiPair});
      histos.add(histodir + "/Photon/hPhi", "hPhi", kTH1F, {axisPhi});
      histos.add(histodir + "/Photon/h3dMass", "h3dMass", kTH3F, {axisCentrality, axisPt, axisPhotonMass});
      histos.add(histodir + "/Photon/hMass", "hMass", kTH1F, {axisPhotonMass});

      histos.add(histodir + "/Lambda/hTrackCode", "hTrackCode", kTH1F, {{11, 0.5f, 11.5f}});
      histos.add(histodir + "/Lambda/hRadius", "hRadius", kTH1F, {axisRadius});
      histos.add(histodir + "/Lambda/hDCADau", "hDCADau", kTH1F, {axisDCAdau});
      histos.add(histodir + "/Lambda/hCosPA", "hCosPA", kTH1F, {axisCosPA});
      histos.add(histodir + "/Lambda/hY", "hY", kTH1F, {axisRapidity});
      histos.add(histodir + "/Lambda/hPosEta", "hPosEta", kTH1F, {axisRapidity});
      histos.add(histodir + "/Lambda/hNegEta", "hNegEta", kTH1F, {axisRapidity});
      histos.add(histodir + "/Lambda/hPosTPCCR", "hPosTPCCR", kTH1F, {axisTPCrows});
      histos.add(histodir + "/Lambda/hNegTPCCR", "hNegTPCCR", kTH1F, {axisTPCrows});
      histos.add(histodir + "/Lambda/hPosITSCls", "hPosITSCls", kTH1F, {axisNCls});
      histos.add(histodir + "/Lambda/hNegITSCls", "hNegITSCls", kTH1F, {axisNCls});
      histos.add(histodir + "/Lambda/hPosChi2PerNc", "hPosChi2PerNc", kTH1F, {axisChi2PerNcl});
      histos.add(histodir + "/Lambda/hNegChi2PerNc", "hNegChi2PerNc", kTH1F, {axisChi2PerNcl});
      histos.add(histodir + "/Lambda/hLifeTime", "hLifeTime", kTH1F, {axisLifetime});
      histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_LambdaPr", "h2dTPCvsTOFNSigma_LambdaPr", kTH2F, {axisTPCNSigma, axisTOFNSigma});
      histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_LambdaPi", "h2dTPCvsTOFNSigma_LambdaPi", kTH2F, {axisTPCNSigma, axisTOFNSigma});
      histos.add(histodir + "/Lambda/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Lambda/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Lambda/hLambdapT", "hLambdapT", kTH1F, {axisPt});
      histos.add(histodir + "/Lambda/hLambdaMass", "hLambdaMass", kTH1F, {axisLambdaMass});
      histos.add(histodir + "/Lambda/h3dLambdaMass", "h3dLambdaMass", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
      histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_ALambdaPr", "h2dTPCvsTOFNSigma_ALambdaPr", kTH2F, {axisTPCNSigma, axisTOFNSigma});
      histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_ALambdaPi", "h2dTPCvsTOFNSigma_ALambdaPi", kTH2F, {axisTPCNSigma, axisTOFNSigma});
      histos.add(histodir + "/Lambda/hALambdaDCANegToPV", "hALambdaDCANegToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Lambda/hALambdaDCAPosToPV", "hALambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
      histos.add(histodir + "/Lambda/hALambdapT", "hALambdapT", kTH1F, {axisPt});
      histos.add(histodir + "/Lambda/hAntiLambdaMass", "hAntiLambdaMass", kTH1F, {axisLambdaMass});
      histos.add(histodir + "/Lambda/h3dAntiLambdaMass", "h3dAntiLambdaMass", kTH3F, {axisCentrality, axisPt, axisLambdaMass});

      histos.add(histodir + "/h2dArmenteros", "h2dArmenteros", kTH2F, {axisAPAlpha, axisAPQt});

      histos.add(histodir + "/Sigma0/hMass", "hMass", kTH1F, {axisSigmaMass});
      histos.add(histodir + "/Sigma0/hPt", "hPt", kTH1F, {axisPt});
      histos.add(histodir + "/Sigma0/hY", "hY", kTH1F, {axisRapidity});
      histos.add(histodir + "/Sigma0/h3dMass", "h3dMass", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add(histodir + "/Sigma0/h3dPhotonRadiusVsMassSigma", "h3dPhotonRadiusVsMassSigma", kTH3F, {axisCentrality, axisRadius, axisSigmaMass});
      histos.add(histodir + "/Sigma0/h2dpTVsOPAngle", "h2dpTVsOPAngle", kTH2F, {axisPt, {140, 0.0f, +7.0f}});

      histos.add(histodir + "/ASigma0/hMass", "hMass", kTH1F, {axisSigmaMass});
      histos.add(histodir + "/ASigma0/hPt", "hPt", kTH1F, {axisPt});
      histos.add(histodir + "/ASigma0/hY", "hY", kTH1F, {axisRapidity});
      histos.add(histodir + "/ASigma0/h3dMass", "h3dMass", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add(histodir + "/ASigma0/h3dPhotonRadiusVsMassSigma", "h3dPhotonRadiusVsMassSigma", kTH3F, {axisCentrality, axisRadius, axisSigmaMass});
      histos.add(histodir + "/ASigma0/h2dpTVsOPAngle", "h2dpTVsOPAngle", kTH2F, {axisPt, {140, 0.0f, +7.0f}});

      // Process MC
      if (doprocessMonteCarlo) {
        histos.add(histodir + "/MC/Photon/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1F, {{2, 0.0f, 2.0f}});
        histos.add(histodir + "/MC/Photon/hPt", "hPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Photon/hMCPt", "hMCPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Photon/h2dPosTPCNSigmaEl", "h2dPosTPCNSigmaEl", kTH2F, {axisPt, axisTPCNSigma});
        histos.add(histodir + "/MC/Photon/h2dNegTPCNSigmaEl", "h2dNegTPCNSigmaEl", kTH2F, {axisPt, axisTPCNSigma});
        histos.add(histodir + "/MC/Photon/h2dPosTPCNSigmaPi", "h2dPosTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
        histos.add(histodir + "/MC/Photon/h2dNegTPCNSigmaPi", "h2dNegTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
        histos.add(histodir + "/MC/Photon/h2dIRVsPt", "h2dIRVsPt", kTH2F, {axisIRBinning, axisPt});
        histos.add(histodir + "/MC/Photon/h3dPAVsIRVsPt", "h3dPAVsIRVsPt", kTH3F, {axisPA, axisIRBinning, axisPt});
        histos.add(histodir + "/MC/Photon/h2dIRVsPt_BadCollAssig", "h2dIRVsPt_BadCollAssig", kTH2F, {axisIRBinning, axisPt});
        histos.add(histodir + "/MC/Photon/h3dPAVsIRVsPt_BadCollAssig", "h3dPAVsIRVsPt_BadCollAssig", kTH3F, {axisPA, axisIRBinning, axisPt});

        histos.add(histodir + "/MC/Lambda/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1F, {{2, 0.0f, 2.0f}});
        histos.add(histodir + "/MC/Lambda/hPt", "hPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Lambda/hMCPt", "hMCPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Lambda/h3dTPCvsTOFNSigma_Pr", "h3dTPCvsTOFNSigma_Pr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
        histos.add(histodir + "/MC/Lambda/h3dTPCvsTOFNSigma_Pi", "h3dTPCvsTOFNSigma_Pi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});

        histos.add(histodir + "/MC/ALambda/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1F, {{2, 0.0f, 2.0f}});
        histos.add(histodir + "/MC/ALambda/hPt", "hPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/ALambda/hMCPt", "hMCPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/ALambda/h3dTPCvsTOFNSigma_Pr", "h3dTPCvsTOFNSigma_Pr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
        histos.add(histodir + "/MC/ALambda/h3dTPCvsTOFNSigma_Pi", "h3dTPCvsTOFNSigma_Pi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});

        histos.add(histodir + "/MC/h2dArmenteros", "h2dArmenteros", kTH2F, {axisAPAlpha, axisAPQt});

        histos.add(histodir + "/MC/Sigma0/hPt", "hPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Sigma0/hMCPt", "hMCPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/Sigma0/h2dMCPtVsLambdaMCPt", "h2dMCPtVsLambdaMCPt", kTH2F, {axisPt, axisPt});
        histos.add(histodir + "/MC/Sigma0/h2dMCPtVsGammaMCPt", "h2dMCPtVsGammaMCPt", kTH2F, {axisPt, axisPt});
        histos.add(histodir + "/MC/Sigma0/hMass", "hMass", kTH1F, {axisSigmaMass});
        histos.add(histodir + "/MC/Sigma0/h3dMass", "h3dMass", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

        histos.add(histodir + "/MC/ASigma0/hPt", "hPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/ASigma0/hMCPt", "hMCPt", kTH1F, {axisPt});
        histos.add(histodir + "/MC/ASigma0/h2dMCPtVsLambdaMCPt", "h2dMCPtVsLambdaMCPt", kTH2F, {axisPt, axisPt});
        histos.add(histodir + "/MC/ASigma0/h2dMCPtVsPhotonMCPt", "h2dMCPtVsPhotonMCPt", kTH2F, {axisPt, axisPt});
        histos.add(histodir + "/MC/ASigma0/hMass", "hMass", kTH1F, {axisSigmaMass});
        histos.add(histodir + "/MC/ASigma0/h3dMass", "h3dMass", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

        // 1/pT Resolution:
        if (fillpTResoQAhistos && histodir == "BeforeSel") {
          histos.add(histodir + "/MC/pTReso/h3dGammaPtResoVsTPCCR", "h3dGammaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h3dGammaPtResoVsTPCCR", "h3dGammaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h2dGammaPtResolution", "h2dGammaPtResolution", kTH2F, {axisInvPt, axisDeltaPt});
          histos.add(histodir + "/MC/pTReso/h2dLambdaPtResolution", "h2dLambdaPtResolution", kTH2F, {axisInvPt, axisDeltaPt});
          histos.add(histodir + "/MC/pTReso/h3dLambdaPtResoVsTPCCR", "h3dLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h3dLambdaPtResoVsTPCCR", "h3dLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h2dAntiLambdaPtResolution", "h2dAntiLambdaPtResolution", kTH2F, {axisInvPt, axisDeltaPt});
          histos.add(histodir + "/MC/pTReso/h3dAntiLambdaPtResoVsTPCCR", "h3dAntiLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h3dAntiLambdaPtResoVsTPCCR", "h3dAntiLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, axisTPCrows});
          histos.add(histodir + "/MC/pTReso/h2dSigma0PtResolution", "h2dSigma0PtResolution", kTH2F, {axisInvPt, axisDeltaPt});
          histos.add(histodir + "/MC/pTReso/h2dAntiSigma0PtResolution", "h2dAntiSigma0PtResolution", kTH2F, {axisInvPt, axisDeltaPt});
        }

        // For background decomposition study
        if (fillBkgQAhistos) {
          histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_All", "h2dPtVsMassSigma_All", kTH2F, {axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_TrueDaughters", "h2dPtVsMassSigma_TrueDaughters", kTH2F, {axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/BkgStudy/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
          histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda", "h2dPtVsMassSigma_TrueGammaFakeLambda", kTH2F, {axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda", "h2dPtVsMassSigma_FakeGammaTrueLambda", kTH2F, {axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_FakeDaughters", "h2dPtVsMassSigma_FakeDaughters", kTH2F, {axisPt, axisSigmaMass});
        }
      }
    }

    // Selections
    histos.add("Selection/Photon/hCandidateSel", "hCandidateSel", kTH1F, {axisCandSel});
    histos.add("Selection/Lambda/hCandidateSel", "hCandidateSel", kTH1F, {axisCandSel});

    for (size_t i = 0; i < PhotonSels.size(); ++i) {
      const auto& sel = PhotonSels[i];

      histos.add(Form("Selection/Photon/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2F, {axisPt, axisPhotonMass});
      histos.get<TH1>(HIST("Selection/Photon/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
      histos.add(Form("Selection/Sigma0/h2dPhoton%s", sel.c_str()), ("h2dPhoton" + sel).c_str(), kTH2F, {axisPt, axisSigmaMass});
    }

    for (size_t i = 0; i < LambdaSels.size(); ++i) {
      const auto& sel = LambdaSels[i];

      histos.add(Form("Selection/Lambda/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2F, {axisPt, axisLambdaMass});
      histos.get<TH1>(HIST("Selection/Lambda/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
      histos.add(Form("Selection/Sigma0/h2dLambda%s", sel.c_str()), ("h2dLambda" + sel).c_str(), kTH2F, {axisPt, axisSigmaMass});
    }
  }

  //__________________________________________
  template <bool isGamma, typename TV0Object>
  int retrieveV0TrackCode(TV0Object const& sigma)
  {

    int TrkCode = 10; // 1: TPC-only, 2: TPC+Something, 3: ITS-Only, 4: ITS+TPC + Something, 10: anything else

    if (isGamma) {
      if (sigma.photonPosTrackCode() == 1 && sigma.photonNegTrackCode() == 1)
        TrkCode = 1;
      if ((sigma.photonPosTrackCode() != 1 && sigma.photonNegTrackCode() == 1) || (sigma.photonPosTrackCode() == 1 && sigma.photonNegTrackCode() != 1))
        TrkCode = 2;
      if (sigma.photonPosTrackCode() == 3 && sigma.photonNegTrackCode() == 3)
        TrkCode = 3;
      if (sigma.photonPosTrackCode() == 2 || sigma.photonNegTrackCode() == 2)
        TrkCode = 4;
    } else {
      if (sigma.lambdaPosTrackCode() == 1 && sigma.lambdaNegTrackCode() == 1)
        TrkCode = 1;
      if ((sigma.lambdaPosTrackCode() != 1 && sigma.lambdaNegTrackCode() == 1) || (sigma.lambdaPosTrackCode() == 1 && sigma.lambdaNegTrackCode() != 1))
        TrkCode = 2;
      if (sigma.lambdaPosTrackCode() == 3 && sigma.lambdaNegTrackCode() == 3)
        TrkCode = 3;
      if (sigma.lambdaPosTrackCode() == 2 || sigma.lambdaNegTrackCode() == 2)
        TrkCode = 4;
    }

    return TrkCode;
  }

  template <typename TV0Object>
  void getpTResolution(TV0Object const& sigma)
  {

    //_______________________________________
    // Gamma MC association
    if (sigma.photonCandPDGCode() == 22) {
      if (sigma.photonMCPt() > 0) {
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.photonNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.photonPosTPCCrossedRows());      // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h2dGammaPtResolution"), 1.f / sigma.photonMCPt(), 1.f / sigma.photonPt() - 1.f / sigma.photonMCPt());                                        // pT resolution
      }
    }

    //_______________________________________
    // Lambda MC association
    if (sigma.lambdaCandPDGCode() == 3122) {
      if (sigma.lambdaMCPt() > 0) {
        histos.fill(HIST("BeforeSel/MC/pTReso/h2dLambdaPtResolution"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt());                                        // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
      }
    }

    //_______________________________________
    // AntiLambda MC association
    if (sigma.lambdaCandPDGCode() == -3122) {
      if (sigma.lambdaMCPt() > 0) {
        histos.fill(HIST("BeforeSel/MC/pTReso/h2dAntiLambdaPtResolution"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt());                                        // pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/pTReso/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
      }
    }

    //_______________________________________
    // Sigma and AntiSigma MC association
    if (sigma.isSigma()) {
      if (sigma.sigmaMCPt() > 0)
        histos.fill(HIST("BeforeSel/MC/pTReso/h2dSigma0PtResolution"), 1.f / sigma.sigmaMCPt(), 1.f / sigma.sigmapT() - 1.f / sigma.sigmaMCPt()); // pT resolution
    }
    if (sigma.isAntiSigma()) {
      if (sigma.sigmaMCPt() > 0)
        histos.fill(HIST("BeforeSel/MC/pTReso/h2dAntiSigma0PtResolution"), 1.f / sigma.sigmaMCPt(), 1.f / sigma.sigmapT() - 1.f / sigma.sigmaMCPt()); // pT resolution
    }
  }

  // To save histograms for background analysis
  template <int mode, typename TV0Object>
  void runBkgAnalysis(TV0Object const& sigma)
  {
    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    bool fIsSigma = sigma.isSigma();
    bool fIsAntiSigma = sigma.isAntiSigma();
    int PhotonPDGCode = sigma.photonCandPDGCode();
    int PhotonPDGCodeMother = sigma.photonCandPDGCodeMother();
    int LambdaPDGCode = sigma.lambdaCandPDGCode();
    int LambdaPDGCodeMother = sigma.lambdaCandPDGCodeMother();
    float sigmapT = sigma.sigmapT();
    float sigmaMass = sigma.sigmaMass();

    histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_All"), sigmapT, sigmaMass);

    //_______________________________________
    // Real Gamma x Real Lambda - but not from the same sigma0/antisigma0!
    if ((PhotonPDGCode == 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122)) && (!fIsSigma && !fIsAntiSigma)) {
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_TrueDaughters"), sigmapT, sigmaMass);
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dTrueDaughtersMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }

    //_______________________________________
    // Real Gamma x fake Lambda
    if ((PhotonPDGCode == 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda"), sigmapT, sigmaMass);

    //_______________________________________
    // Fake Gamma x Real Lambda
    if ((PhotonPDGCode != 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122)))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda"), sigmapT, sigmaMass);

    //_______________________________________
    // Fake Gamma x Fake Lambda
    if ((PhotonPDGCode != 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_FakeDaughters"), sigmapT, sigmaMass);
  }

  template <int mode, typename TV0Object>
  void fillQAHistos(TV0Object const& sigma)
  {

    // Check whether it is before or after selections
    // static std::string main_dir;
    // main_dir = IsBeforeSel ? "BeforeSel" : "AfterSel";
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    // Get V0trackCode
    int GammaTrkCode = retrieveV0TrackCode<true>(sigma);
    int LambdaTrkCode = retrieveV0TrackCode<false>(sigma);

    float photonRZLineCut = TMath::Abs(sigma.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-PhotonMaxDauEta))) - PhotonLineCutZ0;
    //_______________________________________
    // Photon
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hTrackCode"), GammaTrkCode);
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hV0Type"), sigma.photonV0Type());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegpT"), sigma.photonNegPt());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPospT"), sigma.photonPosPt());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCANegToPV"), sigma.photonDCANegPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCAPosToPV"), sigma.photonDCAPosPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCADau"), sigma.photonDCADau());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosTPCCR"), sigma.photonPosTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegTPCCR"), sigma.photonNegTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dPosTPCNSigmaEl"), sigma.photonPosPt(), sigma.photonPosTPCNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dNegTPCNSigmaEl"), sigma.photonNegPt(), sigma.photonNegTPCNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dPosTPCNSigmaPi"), sigma.photonPosPt(), sigma.photonPosTPCNSigmaPi());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dNegTPCNSigmaPi"), sigma.photonNegPt(), sigma.photonNegTPCNSigmaPi());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hpT"), sigma.photonPt());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hY"), sigma.photonY());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosEta"), sigma.photonPosEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegEta"), sigma.photonNegEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hRadius"), sigma.photonRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hZ"), sigma.photonZconv());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZCut"), sigma.photonRadius(), photonRZLineCut);
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZPlane"), sigma.photonZconv(), sigma.photonRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hCosPA"), sigma.photonCosPA());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPsiPair"), sigma.photonPsiPair());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPhi"), sigma.photonPhi());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h3dMass"), sigma.sigmaCentrality(), sigma.photonPt(), sigma.photonMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hMass"), sigma.photonMass());

    //_______________________________________
    // Lambdas
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hTrackCode"), LambdaTrkCode);
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hRadius"), sigma.lambdaRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hDCADau"), sigma.lambdaDCADau());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hCosPA"), sigma.lambdaCosPA());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hY"), sigma.lambdaY());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosEta"), sigma.lambdaPosEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegEta"), sigma.lambdaNegEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosTPCCR"), sigma.lambdaPosTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegTPCCR"), sigma.lambdaNegTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosITSCls"), sigma.lambdaPosITSCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegITSCls"), sigma.lambdaNegITSCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosChi2PerNc"), sigma.lambdaPosChi2PerNcl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegChi2PerNc"), sigma.lambdaNegChi2PerNcl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLifeTime"), sigma.lambdaLifeTime());

    //_______________________________________
    // Sigmas and Lambdas
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

    if (sigma.lambdaAlpha() > 0) {
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_LambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_LambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaDCANegToPV"), sigma.lambdaDCANegPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaDCAPosToPV"), sigma.lambdaDCAPosPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdapT"), sigma.lambdaPt());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaMass"), sigma.lambdaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h3dLambdaMass"), sigma.sigmaCentrality(), sigma.lambdaPt(), sigma.lambdaMass());

      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hMass"), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hPt"), sigma.sigmapT());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hY"), sigma.sigmaRapidity());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dMass"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dPhotonRadiusVsMassSigma"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h2dpTVsOPAngle"), sigma.sigmapT(), sigma.sigmaOPAngle());
    } else {
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_ALambdaPr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_ALambdaPi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdaDCANegToPV"), sigma.lambdaDCANegPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdaDCAPosToPV"), sigma.lambdaDCAPosPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdapT"), sigma.lambdaPt());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hAntiLambdaMass"), sigma.antilambdaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h3dAntiLambdaMass"), sigma.sigmaCentrality(), sigma.lambdaPt(), sigma.antilambdaMass());

      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hMass"), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hPt"), sigma.sigmapT());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hY"), sigma.sigmaRapidity());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dMass"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dPhotonRadiusVsMassSigma"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h2dpTVsOPAngle"), sigma.sigmapT(), sigma.sigmaOPAngle());
    }

    //_______________________________________
    // MC specific
    if (doprocessMonteCarlo) {
      if constexpr (requires { sigma.lambdaCandPDGCode(); sigma.photonCandPDGCode(); }) {

        //_______________________________________
        // Gamma MC association
        if (sigma.photonCandPDGCode() == 22) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hV0ToCollAssoc"), sigma.photonIsCorrectlyAssoc());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPt"), sigma.photonPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hMCPt"), sigma.photonMCPt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPosTPCNSigmaEl"), sigma.photonPosPt(), sigma.photonPosTPCNSigmaEl());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dNegTPCNSigmaEl"), sigma.photonNegPt(), sigma.photonNegTPCNSigmaEl());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPosTPCNSigmaPi"), sigma.photonPosPt(), sigma.photonPosTPCNSigmaPi());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dNegTPCNSigmaPi"), sigma.photonNegPt(), sigma.photonNegTPCNSigmaPi());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dIRVsPt"), sigma.sigmaIR(), sigma.photonMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h3dPAVsIRVsPt"), TMath::ACos(sigma.photonCosPA()), sigma.sigmaIR(), sigma.photonMCPt());

          if (!sigma.photonIsCorrectlyAssoc()) {
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dIRVsPt_BadCollAssig"), sigma.sigmaIR(), sigma.photonMCPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h3dPAVsIRVsPt_BadCollAssig"), TMath::ACos(sigma.photonCosPA()), sigma.sigmaIR(), sigma.photonMCPt());
          }
        }

        //_______________________________________
        // Lambda MC association
        if (sigma.lambdaCandPDGCode() == 3122) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hV0ToCollAssoc"), sigma.lambdaIsCorrectlyAssoc());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hPt"), sigma.lambdaPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hMCPt"), sigma.lambdaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/h3dTPCvsTOFNSigma_Pr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/h3dTPCvsTOFNSigma_Pi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());
        }

        //_______________________________________
        // AntiLambda MC association
        if (sigma.lambdaCandPDGCode() == -3122) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hV0ToCollAssoc"), sigma.lambdaIsCorrectlyAssoc());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hPt"), sigma.lambdaPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hMCPt"), sigma.lambdaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/h3dTPCvsTOFNSigma_Pr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/h3dTPCvsTOFNSigma_Pi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma(), sigma.lambdaPt());
        }

        //_______________________________________
        // Sigma0 MC association
        if (sigma.isSigma()) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hPt"), sigma.sigmapT());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hMCPt"), sigma.sigmaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h2dMCPtVsLambdaMCPt"), sigma.sigmaMCPt(), sigma.lambdaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h2dMCPtVsGammaMCPt"), sigma.sigmaMCPt(), sigma.photonMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hMass"), sigma.sigmaMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h3dMass"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
        }

        //_______________________________________
        // AntiSigma0 MC association
        if (sigma.isAntiSigma()) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hPt"), sigma.sigmapT());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hMCPt"), sigma.sigmaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h2dMCPtVsLambdaMCPt"), sigma.sigmaMCPt(), sigma.lambdaMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h2dMCPtVsPhotonMCPt"), sigma.sigmaMCPt(), sigma.photonMCPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hMass"), sigma.sigmaMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h3dMass"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
        }

        // For background studies:
        if (fillBkgQAhistos)
          runBkgAnalysis<mode>(sigma);

        //_______________________________________
        // pT resolution histos
        if ((mode == 0) && fillpTResoQAhistos)
          getpTResolution(sigma);
      }
    }
  }

  template <int selection_index, typename TV0Object>
  void fillSelHistos(TV0Object const& sigma, int PDGRequired)
  {

    static constexpr std::string_view PhotonSelsLocal[] = {"NoSel", "V0Type", "DaupT", "DCADauToPV",
                                                           "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                           "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                           "PsiPair", "Phi", "Mass"};

    static constexpr std::string_view LambdaSelsLocal[] = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                           "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                           "TPCTOFPID", "DCADauToPV", "Mass"};

    if (PDGRequired == 22) {
      if constexpr (selection_index >= 0 && selection_index < (int)std::size(PhotonSelsLocal)) {
        histos.fill(HIST("Selection/Photon/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/Photon/h2d") + HIST(PhotonSelsLocal[selection_index]), sigma.photonPt(), sigma.photonMass());
        histos.fill(HIST("Selection/Sigma0/h2dPhoton") + HIST(PhotonSelsLocal[selection_index]), sigma.sigmapT(), sigma.sigmaMass());
      }
    }

    if (PDGRequired == 3122) {
      if constexpr (selection_index >= 0 && selection_index < (int)std::size(LambdaSelsLocal)) {
        histos.fill(HIST("Selection/Lambda/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/Lambda/h2d") + HIST(LambdaSelsLocal[selection_index]), sigma.lambdaPt(), sigma.lambdaMass());
        histos.fill(HIST("Selection/Sigma0/h2dLambda") + HIST(LambdaSelsLocal[selection_index]), sigma.sigmapT(), sigma.sigmaMass());
      }
    }
  }

  // Apply specific selections for photons
  template <typename TV0Object>
  bool selectPhoton(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, 22);
    if (cand.photonV0Type() != Photonv0TypeSel && Photonv0TypeSel > -1)
      return false;

    fillSelHistos<1>(cand, 22);
    if ((cand.photonPosPt() < PhotonDauMinPt) || (cand.photonNegPt() < PhotonDauMinPt))
      return false;

    fillSelHistos<2>(cand, 22);
    if ((TMath::Abs(cand.photonDCAPosPV()) < PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < PhotonMinDCADauToPv))
      return false;

    fillSelHistos<3>(cand, 22);
    if (TMath::Abs(cand.photonDCADau()) > PhotonMaxDCAV0Dau)
      return false;

    fillSelHistos<4>(cand, 22);
    if ((cand.photonPosTPCCrossedRows() < PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < PhotonMinTPCCrossedRows))
      return false;

    fillSelHistos<5>(cand, 22);
    if (((cand.photonPosTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
      return false;

    if (((cand.photonNegTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
      return false;

    fillSelHistos<6>(cand, 22);
    if ((cand.photonPt() < PhotonMinPt) || (cand.photonPt() > PhotonMaxPt))
      return false;

    fillSelHistos<7>(cand, 22);
    if ((TMath::Abs(cand.photonY()) > PhotonMaxRap) || (TMath::Abs(cand.photonPosEta()) > PhotonMaxDauEta) || (TMath::Abs(cand.photonNegEta()) > PhotonMaxDauEta))
      return false;

    fillSelHistos<8>(cand, 22);
    if ((cand.photonRadius() < PhotonMinRadius) || (cand.photonRadius() > PhotonMaxRadius))
      return false;

    fillSelHistos<9>(cand, 22);
    float photonRZLineCut = TMath::Abs(cand.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-PhotonMaxDauEta))) - PhotonLineCutZ0;
    if ((TMath::Abs(cand.photonRadius()) < photonRZLineCut) || (TMath::Abs(cand.photonZconv()) > PhotonMaxZ))
      return false;

    fillSelHistos<10>(cand, 22);
    if (cand.photonQt() > PhotonMaxQt)
      return false;

    if (TMath::Abs(cand.photonAlpha()) > PhotonMaxAlpha)
      return false;

    fillSelHistos<11>(cand, 22);
    if (cand.photonCosPA() < PhotonMinV0cospa)
      return false;

    fillSelHistos<12>(cand, 22);
    if (TMath::Abs(cand.photonPsiPair()) > PhotonPsiPairMax)
      return false;

    fillSelHistos<13>(cand, 22);
    if ((((cand.photonPhi() > PhotonPhiMin1) && (cand.photonPhi() < PhotonPhiMax1)) || ((cand.photonPhi() > PhotonPhiMin2) && (cand.photonPhi() < PhotonPhiMax2))) && ((PhotonPhiMin1 != -1) && (PhotonPhiMax1 != -1) && (PhotonPhiMin2 != -1) && (PhotonPhiMax2 != -1)))
      return false;

    fillSelHistos<14>(cand, 22);
    if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
      return false;

    fillSelHistos<15>(cand, 22);
    return true;
  }

  // Apply specific selections for lambdas
  template <typename TV0Object>
  bool selectLambda(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, 3122);
    if ((cand.lambdaRadius() < LambdaMinv0radius) || (cand.lambdaRadius() > LambdaMaxv0radius))
      return false;

    fillSelHistos<1>(cand, 3122);
    if (TMath::Abs(cand.lambdaDCADau()) > LambdaMaxDCAV0Dau)
      return false;

    fillSelHistos<2>(cand, 3122);
    if ((cand.lambdaQt() < LambdaMinQt) || (cand.lambdaQt() > LambdaMaxQt))
      return false;

    if ((TMath::Abs(cand.lambdaAlpha()) < LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > LambdaMaxAlpha))
      return false;

    fillSelHistos<3>(cand, 3122);
    if (cand.lambdaCosPA() < LambdaMinv0cospa)
      return false;

    fillSelHistos<4>(cand, 3122);
    if ((TMath::Abs(cand.lambdaY()) > LambdaMaxRap) || (TMath::Abs(cand.lambdaPosEta()) > LambdaMaxDauEta) || (TMath::Abs(cand.lambdaNegEta()) > LambdaMaxDauEta))
      return false;

    fillSelHistos<5>(cand, 3122);
    if ((cand.lambdaPosTPCCrossedRows() < LambdaMinTPCCrossedRows) || (cand.lambdaNegTPCCrossedRows() < LambdaMinTPCCrossedRows))
      return false;

    fillSelHistos<6>(cand, 3122);
    // check minimum number of ITS clusters + reject ITS afterburner tracks if requested
    bool posIsFromAfterburner = cand.lambdaPosChi2PerNcl() < 0;
    bool negIsFromAfterburner = cand.lambdaNegChi2PerNcl() < 0;
    if (cand.lambdaPosITSCls() < LambdaMinITSclusters && (!LambdaRejectPosITSafterburner || posIsFromAfterburner))
      return false;
    if (cand.lambdaNegITSCls() < LambdaMinITSclusters && (!LambdaRejectNegITSafterburner || negIsFromAfterburner))
      return false;

    fillSelHistos<7>(cand, 3122);
    if (cand.lambdaLifeTime() > LambdaMaxLifeTime)
      return false;

    // Separating lambda and antilambda selections:
    fillSelHistos<8>(cand, 3122);
    if (cand.lambdaAlpha() > 0) { // Lambda selection
      // TPC Selection
      if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPrTPCNSigma()) > LambdaMaxTPCNSigmas))
        return false;
      if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPiTPCNSigma()) > LambdaMaxTPCNSigmas))
        return false;

      // TOF Selection
      if (fselLambdaTOFPID && (TMath::Abs(cand.lambdaPrTOFNSigma()) > LambdaPrMaxTOFNSigmas))
        return false;
      if (fselLambdaTOFPID && (TMath::Abs(cand.lambdaPiTOFNSigma()) > LambdaPiMaxTOFNSigmas))
        return false;

      // DCA Selection
      fillSelHistos<9>(cand, 3122);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < LambdaMinDCANegToPv))
        return false;

      // Mass Selection
      fillSelHistos<10>(cand, 3122);
      if (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) > LambdaWindow)
        return false;

      fillSelHistos<11>(cand, 3122);

    } else { // AntiLambda selection

      // TPC Selection
      if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPiTPCNSigma()) > LambdaMaxTPCNSigmas))
        return false;
      if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPrTPCNSigma()) > LambdaMaxTPCNSigmas))
        return false;

      // TOF Selection
      if (fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPrTOFNSigma()) > LambdaPrMaxTOFNSigmas))
        return false;
      if (fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPiTOFNSigma()) > LambdaPiMaxTOFNSigmas))
        return false;

      // DCA Selection
      fillSelHistos<9>(cand, 3122);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < ALambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < ALambdaMinDCANegToPv))
        return false;

      // Mass Selection
      fillSelHistos<10>(cand, 3122);
      if (TMath::Abs(cand.antilambdaMass() - o2::constants::physics::MassLambda0) > LambdaWindow)
        return false;

      fillSelHistos<11>(cand, 3122);
    }

    return true;
  }

  // Apply selections in sigma0 candidates
  template <typename TV0Object>
  bool processSigmaCandidate(TV0Object const& cand)
  {
    // Optionally Select on Interaction Rate
    if (fGetIR && (maxIR != -1) && (minIR != -1) && ((cand.sigmaIR() <= minIR) || (cand.sigmaIR() >= maxIR))) {
      return false;
    }

    // Do ML analysis
    if (fUseMLSel) {
      if ((cand.gammaBDTScore() == -1) || (cand.lambdaBDTScore() == -1) || (cand.antilambdaBDTScore() == -1)) {
        LOGF(fatal, "ML Score is not available! Please, enable gamma and lambda selection with ML in sigmabuilder!");
      }
      // Photon selection:
      if (cand.gammaBDTScore() <= Gamma_MLThreshold)
        return false;

      // Lambda selection:
      if (cand.lambdaBDTScore() <= Lambda_MLThreshold)
        return false;

      // AntiLambda selection:
      if (cand.antilambdaBDTScore() <= AntiLambda_MLThreshold)
        return false;

    }

    // Go for standard analysis
    else {

      // Photon specific selections
      if (!selectPhoton(cand))
        return false;

      // Lambda specific selections
      if (!selectLambda(cand))
        return false;

      // Sigma0 specific selections
      if (TMath::Abs(cand.sigmaRapidity()) > SigmaMaxRap)
        return false;
    }

    return true;
  }

  void processMonteCarlo(V0MCSigmas const& sigmas)
  {
    for (auto& sigma : sigmas) { // selecting Sigma0-like candidates
      if (doMCAssociation && !(sigma.isSigma() || sigma.isAntiSigma())) {
        continue;
      }

      // Fill histos before any selection
      fillQAHistos<0>(sigma);

      // Select sigma0 candidates
      if (!processSigmaCandidate(sigma))
        continue;

      // Fill histos after all selections
      fillQAHistos<1>(sigma);
    }
  }

  void processRealData(V0Sigmas const& sigmas)
  {
    for (auto& sigma : sigmas) { // selecting Sigma0-like candidates

      // Fill histos before any selection
      fillQAHistos<0>(sigma);

      // Select sigma0 candidates
      if (!processSigmaCandidate(sigma))
        continue;

      // Fill histos after all selections
      fillQAHistos<1>(sigma);
    }
  }

  PROCESS_SWITCH(sigmaanalysis, processMonteCarlo, "Do Monte-Carlo-based analysis", false);
  PROCESS_SWITCH(sigmaanalysis, processRealData, "Do real data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigmaanalysis>(cfgc)};
}
