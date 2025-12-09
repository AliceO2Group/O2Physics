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

/// \file jetSpectraEseTask.cxx
/// \brief jet spectra analysis framework with ESE (19/08/2024)
///
/// \author Joachim C. K. B. Hansen, Lund University

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/Qvectors.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetSpectraEseTask {
  ConfigurableAxis binJetPt{"binJetPt", {250, -50., 200.}, ""};
  ConfigurableAxis bindPhi{"bindPhi", {180, -o2::constants::math::PI, o2::constants::math::PI}, ""};
  ConfigurableAxis binESE{"binESE", {100, 0, 100}, ""};
  ConfigurableAxis binCos{"binCos", {100, -1.05, 1.05}, ""};
  ConfigurableAxis binOccupancy{"binOccupancy", {5000, 0, 25000}, ""};
  ConfigurableAxis binQVec{"binQVec", {500, -3, 3}, ""};
  ConfigurableAxis binCentrality{"binCentrality", {101, -1, 100}, ""};
  ConfigurableAxis binPhi{"binPhi", {60, -1.0, 7.0}, ""};
  ConfigurableAxis binEta{"binEta", {80, -0.9, 0.9}, ""};
  ConfigurableAxis binFit0{"binFit0", {100, 0, 50}, ""};
  ConfigurableAxis binFit13{"binFit13", {100, 0, 1.4}, ""};
  ConfigurableAxis binFit24{"binFit24", {100, 0, 10}, ""};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.2, "jet resolution parameter"};
  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};
  Configurable<float> randomConeLeadJetDeltaR{"randomConeLeadJetDeltaR", -99.0, "min distance between leading jet axis and random cone (RC) axis; if negative, min distance is set to automatic value of R_leadJet+R_RC "};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "vertex z cut"};
  Configurable<std::vector<float>> centRange{"centRange", {0, 90}, "centrality region of interest"};
  Configurable<bool> cfgSelCentrality{"cfgSelCentrality", true, "Flag for centrality selection"};
  // Configurable<double> leadingTrackPtCut{"leadingTrackPtCut", 5.0, "leading jet pT cut"};
  Configurable<double> jetAreaFractionMin{"jetAreaFractionMin", -99, "used to make a cut on the jet areas"};
  Configurable<bool> cfgCentVariant{"cfgCentVariant", false, "Flag for centrality variant 1"};
  Configurable<bool> cfgisPbPb{"cfgisPbPb", false, "Flag for using MC centrality in PbPb"};
  Configurable<bool> cfgbkgSubMC{"cfgbkgSubMC", true, "Flag for MC background subtraction"};
  Configurable<bool> cfgUseMCEventWeights{"cfgUseMCEventWeights", false, "Flag for using MC event weights"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<bool> checkLeadConstituentMinPtForMcpJets{"checkLeadConstituentMinPtForMcpJets", false, "flag to choose whether particle level jets should have their lead track pt above leadingConstituentPtMin to be accepted; off by default, as leadingConstituentPtMin cut is only applied on MCD jets for the Pb-Pb analysis using pp MC anchored to Pb-Pb for the response matrix"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 200.0, "maximum pT acceptance for tracks"};

  Configurable<float> jetEtaMin{"jetEtaMin", -0.7, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "maximum jet pseudorapidity"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8FullPbPb", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<bool> cfgEvSelOccupancy{"cfgEvSelOccupancy", true, "Flag for occupancy cut"};

  Configurable<std::vector<int>> cfgCutOccupancy{"cfgCutOccupancy", {0, 1000}, "Occupancy cut"};
  Configurable<std::vector<float>> cfgOccupancyPtCut{"cfgOccupancyPtCut", {0, 100}, "pT cut"};

  Configurable<bool> cfgrhoPhi{"cfgrhoPhi", true, "Flag for rho(phi)"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number // look in Qvector table for this number"};
  Configurable<int> cfgnCorrLevel{"cfgnCorrLevel", 3, "QVector step: 0 = no corr, 1 = rect, 2 = twist, 3 = full"};

  Configurable<std::string> cfgEPRefA{"cfgEPRefA", "FT0A", "EP reference A"};
  Configurable<std::string> cfgEPRefB{"cfgEPRefB", "TPCpos", "EP reference B"};
  Configurable<std::string> cfgEPRefC{"cfgEPRefC", "TPCneg", "EP reference C"};

  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T,jet}"};
  AxisSpec dPhiAxis = {bindPhi, "#Delta#phi"};
  AxisSpec eseAxis = {binESE, "#it{q}_{2}"};
  AxisSpec cosAxis = {binCos, ""};
  AxisSpec occAxis = {binOccupancy, "Occupancy"};
  AxisSpec qvecAxis = {binQVec, "Q-vector"};
  AxisSpec centAxis = {binCentrality, "Centrality"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec fitAxis0 = {binFit0, "Fit0"};
  AxisSpec fitAxis13 = {binFit13, "Fit13"};
  AxisSpec fitAxis24 = {binFit24, "Fit24"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  std::vector<int> eventSelectionBits;
  int trackSelection{-1};

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && nabs(aod::jet::eta) < 0.9f - jetR;
  Filter colFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using ChargedMCDJets = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>>;
  Preslice<ChargedMCDJets> mcdjetsPerJCollision = o2::aod::jet::collisionId;

  enum class DetID { FT0C,
                     FT0A,
                     FT0M,
                     FV0A,
                     TPCpos,
                     TPCneg,
                     TPCall };

  struct EventPlane {
    float psi2;
    float psi3;
  };

  struct EventPlaneFiller {
    bool psi;
    bool hist;
  };

  static constexpr EventPlaneFiller PsiFillerEP = {true, true};
  static constexpr EventPlaneFiller PsiFillerEse = {true, false};
  static constexpr EventPlaneFiller PsiFillerOcc = {false, false};

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    LOGF(info, "jetSpectraEse::init()");

    if (doprocessESEDataCharged) {
      LOGF(info, "JetSpectraEseTask::init() - ESE Data Process");
      registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("hCentralitySel", ";Centrality;entries", {HistType::kTH1F, {{centAxis}}});
      registry.add("hCentralityAnalyzed", ";Centrality;entries", {HistType::kTH1F, {{centAxis}}});
      registry.add("hTrackCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("hJetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
      registry.add("hJetPt_bkgsub", "jet pT background sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
      registry.add("hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{etaAxis}}});
      registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{phiAxis}}});
      registry.add("hRho", ";#rho;entries", {HistType::kTH2F, {{centAxis}, {100, 0, 200.}}});
      registry.add("hRhoPhi", ";#rho;entries", {HistType::kTH2F, {{centAxis}, {100, 0, 200.}}});
      registry.add("hJetArea", ";area_{jet};entries", {HistType::kTH1F, {{80, 0, 10.}}});
      registry.add("hCentJetPtdPhiq2", "", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {dPhiAxis}, {eseAxis}}});
      registry.add("hCentJetPtdPhiq2RhoPhi", "", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {dPhiAxis}, {eseAxis}}});

      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(1, "Input event");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(2, "Event selection");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(3, "Occupancy cut");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(4, "ESE available");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(5, "Lead track");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(6, "Jet loop");
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(7, "Centrality analyzed");

      registry.add("hPhiPtsum", "jet sumpt;sum p_{T};entries", {HistType::kTH1F, {{40, 0., o2::constants::math::TwoPI}}});
      registry.add("hfitPar0", "", {HistType::kTH2F, {{centAxis}, {fitAxis0}}});
      registry.add("hfitPar1", "", {HistType::kTH2F, {{centAxis}, {fitAxis13}}});
      registry.add("hfitPar2", "", {HistType::kTH2F, {{centAxis}, {fitAxis24}}});
      registry.add("hfitPar3", "", {HistType::kTH2F, {{centAxis}, {fitAxis13}}});
      registry.add("hfitPar4", "", {HistType::kTH2F, {{centAxis}, {fitAxis24}}});
      registry.add("hPValueCentCDF", "p-value cDF vs centrality; centrality; p-value", {HistType::kTH2F, {{centAxis}, {40, 0, 1}}});
      registry.add("hCentChi2Ndf", "Chi2 vs centrality; centrality; #tilde{#chi^{2}}", {HistType::kTH2F, {{centAxis}, {100, 0, 5}}});
      registry.add("hCentPhi", "centrality vs #rho(#varphi); centrality;  #rho(#varphi) ", {HistType::kTH2F, {{centAxis}, {210, -10.0, 200.0}}});
      registry.add("hdPhiRhoPhi", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho(#varphi) ", {HistType::kTH2F, {{40, -o2::constants::math::PI, o2::constants::math::PI}, {210, -10.0, 200.0}}});
    }
    if (doprocessESEEPData) {
      LOGF(info, "JetSpectraEseTask::init() - Event Plane Process");
      registry.add("hPsi2FT0C", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{centAxis}, {150, -2.5, 2.5}}});
      registry.addClone("hPsi2FT0C", "hPsi2FT0A");
      registry.addClone("hPsi2FT0C", "hPsi2FV0A");
      registry.addClone("hPsi2FT0C", "hPsi2TPCpos");
      registry.addClone("hPsi2FT0C", "hPsi2TPCneg");
      registry.add("hCosPsi2AmC", ";Centrality;cos(2(#Psi_{2}^{A}-#Psi_{2}^{B}));#it{q}_{2}", {HistType::kTH3F, {{centAxis}, {cosAxis}, {eseAxis}}});
      registry.addClone("hCosPsi2AmC", "hCosPsi2AmB");
      registry.addClone("hCosPsi2AmC", "hCosPsi2BmC");
      registry.add("hQvecUncorV2", ";Centrality;Q_x;Q_y", {HistType::kTH3F, {{centAxis}, {qvecAxis}, {qvecAxis}}});
      registry.addClone("hQvecUncorV2", "hQvecRectrV2");
      registry.addClone("hQvecUncorV2", "hQvecTwistV2");
      registry.addClone("hQvecUncorV2", "hQvecFinalV2");
      registry.addClone("hPsi2FT0C", "hEPUncorV2");
      registry.addClone("hPsi2FT0C", "hEPRectrV2");
      registry.addClone("hPsi2FT0C", "hEPTwistV2");
    }
    if (doprocessESEBackground) {
      LOGF(info, "JetSpectraEseTask::init() - Background Process");
      registry.add("hCentRhoRandomCone", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("hCentRhoRandomConeRandomTrackDir", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("hCentRhoRandomConewoLeadingJet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTHnSparseF, {{centAxis}, {800, -400.0, 400.0}, {dPhiAxis}, {eseAxis}}});
      registry.add("hCentRhoRandomConeRndTrackDirwoOneLeadingJet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTHnSparseF, {{centAxis}, {800, -400.0, 400.0}, {dPhiAxis}, {eseAxis}}});
      registry.add("hCentRhoRandomConeRndTrackDirwoTwoLeadingJet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTHnSparseF, {{centAxis}, {800, -400.0, 400.0}, {dPhiAxis}, {eseAxis}}});
    }
    if (doprocessMCParticleLevel) {
      LOGF(info, "JetSpectraEseTask::init() - MC Particle level");
      registry.add("mcp/hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("mcp/hCentralitySel", ";centr;entries", {HistType::kTH1F, {{centAxis}}});
      // registry.add("mcp/hJetSparse", ";Centrality;#it{p}_{T,jet part}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}});
      registry.add("mcp/hJetSparse", ";Centrality;#it{p}_{T,jet part}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}});
      // registry.add("mcp/hJetPt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
      // registry.add("mcp/hJetEta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{etaAxis}}});
      // registry.add("mcp/hJetPhi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{phiAxis}}});

      registry.get<TH1>(HIST("mcp/hEventCounter"))->GetXaxis()->SetBinLabel(1, "Input event");
      registry.get<TH1>(HIST("mcp/hEventCounter"))->GetXaxis()->SetBinLabel(2, "Collision size < 1");
      registry.get<TH1>(HIST("mcp/hEventCounter"))->GetXaxis()->SetBinLabel(3, "MCD size != 1");
      registry.get<TH1>(HIST("mcp/hEventCounter"))->GetXaxis()->SetBinLabel(4, "Occupancy cut");
    }
    if (doprocessMCDetectorLevel) {
      LOGF(info, "JetSpectraEseTask::init() - MC Detector level");
      registry.add("mcd/hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("mcd/hCentralitySel", ";centr;entries", {HistType::kTH1F, {{centAxis}}});
      // registry.add("mcd/hJetPt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
      // registry.add("mcd/hJetSparse", ";Centrality;#it{p}_{T,jet det}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}});
      registry.add("mcd/hJetSparse", ";Centrality;#it{p}_{T,jet det}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}});
      // registry.add("mcd/hJetEta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{etaAxis}}});
      // registry.add("mcd/hJetPhi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{phiAxis}}});

      registry.get<TH1>(HIST("mcd/hEventCounter"))->GetXaxis()->SetBinLabel(1, "Input event");
      registry.get<TH1>(HIST("mcd/hEventCounter"))->GetXaxis()->SetBinLabel(2, "Collision size < 1");
      registry.get<TH1>(HIST("mcd/hEventCounter"))->GetXaxis()->SetBinLabel(3, "Occupancy cut");
    }
    if (doprocessMCChargedMatched) {
      LOGF(info, "JetSpectraEseTask::init() - MC Charged Matched");

      registry.add("mcm/hJetSparse", ";Centrality;#it{p}_{T,jet det}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}}); /* detector level */
      // registry.add("mcm/hPartSparseMatch", ";Centrality;#it{p}_{T,jet part}; #eta; #phi", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {etaAxis}, {phiAxis}}});
      registry.addClone("mcm/hJetSparse", "mcm/hDetSparseMatch");
      registry.addClone("mcm/hJetSparse", "mcm/hPartSparseMatch");

      registry.add("mcm/hMCEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("mcm/hMCDMatchedEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("mcm/hCentralityAnalyzed", ";Centrality;entries", {HistType::kTH1F, {{centAxis}}});

      registry.add("mcm/hMatchedJetsPtDelta", "#it{p}_{T,jet part}; det - part", {HistType::kTH2F, {{jetPtAxis}, {100, -20., 20.0}}});
      registry.add("mcm/hGenMatchedJetsPtDeltadPt", "#it{p}_{T,jet part}; det - part / part", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {200, -20., 20.0}}});
      registry.add("mcm/hRecoMatchedJetsPtDeltadPt", "#it{p}_{T,jet det}; det - part / det", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {200, -20., 20.0}}});
      registry.add("mcm/hMatchedJetsEtaDelta", "#eta_{jet part}; det - part", {HistType::kTH2F, {{etaAxis}, {200, -0.8, 0.8}}});
      registry.add("mcm/hMatchedJetsPhiDelta", "#phi_{jet part}; det - part", {HistType::kTH2F, {{phiAxis}, {200, -10.0, 10.}}});

      registry.add("mcm/hRespMcDMcPMatch", ";Centrality,#it{p}_{T, jet det}; #it{p}_{T, jet part}", {HistType::kTHnSparseF, {{centAxis}, {jetPtAxis}, {jetPtAxis}}});

      registry.get<TH1>(HIST("mcm/hMCEventCounter"))->GetXaxis()->SetBinLabel(1, "Input event");
      registry.get<TH1>(HIST("mcm/hMCEventCounter"))->GetXaxis()->SetBinLabel(2, "Collision size < 1");
      registry.get<TH1>(HIST("mcm/hMCEventCounter"))->GetXaxis()->SetBinLabel(3, "Vertex cut");
      registry.get<TH1>(HIST("mcm/hMCEventCounter"))->GetXaxis()->SetBinLabel(4, "After analysis");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(1, "Input event");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(2, "Vertex cut");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(3, "Event selection");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(4, "Occupancy cut");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(5, "Centrality cut1:cut2");
      registry.get<TH1>(HIST("mcm/hMCDMatchedEventCounter"))->GetXaxis()->SetBinLabel(6, "After analysis");
    }
    if (doprocessESEOccupancy) {
      LOGF(info, "JetSpectraEseTask::init() - Occupancy QA");
      registry.add("hEventCounterOcc", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("hTrackPt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTHnSparseF, {{centAxis}, {100, 0, 100}, {eseAxis}, {occAxis}}});
      registry.add("hTrackEta", "track #eta;#eta_{track};entries", {HistType::kTH3F, {{centAxis}, {etaAxis}, {occAxis}}});
      registry.add("hTrackPhi", "track #phi;#phi_{track};entries", {HistType::kTH3F, {{centAxis}, {phiAxis}, {occAxis}}});
      registry.add("hOccupancy", "Occupancy;Occupancy;entries", {HistType::kTH1F, {{occAxis}}});
      registry.add("hPsiOccupancy", "Occupancy;#Psi_{2};entries", {HistType::kTH3F, {{centAxis}, {150, -2.5, 2.5}, {occAxis}}});
    }
  }

  void processESEDataCharged(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors, aod::QPercentileFT0Cs>::iterator const& collision,
                             soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets,
                             aod::JetTracks const& tracks)
  {
    float counter{0.5f};
    registry.fill(HIST("hEventCounter"), counter++);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    if (cfgEvSelOccupancy && !isOccupancyAccepted(collision))
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    auto centrality = cfgCentVariant ? collision.centFT0CVariant1() : collision.centFT0M();
    if (cfgSelCentrality && !isCentralitySelected(centrality))
      return;
    registry.fill(HIST("hCentralitySel"), centrality);

    const auto psi{procEP<PsiFillerEse>(collision)};
    const auto qPerc{collision.qPERCFT0C()};
    if (qPerc[0] < 0)
      return;
    registry.fill(HIST("hEventCounter"), counter++);
    std::unique_ptr<TF1> rhoFit{nullptr};
    if (cfgrhoPhi) {
      rhoFit = fitRho<true>(collision, psi, tracks, jets);
      if (!rhoFit)
        return;
    }

    registry.fill(HIST("hEventCounter"), counter++);
    registry.fill(HIST("hRho"), centrality, collision.rho());
    registry.fill(HIST("hCentralityAnalyzed"), centrality);
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      registry.fill(HIST("hJetPt"), jet.pt());
      registry.fill(HIST("hJetPt_bkgsub"), jet.pt() - collision.rho() * jet.area());
      registry.fill(HIST("hJetEta"), jet.eta());
      registry.fill(HIST("hJetPhi"), jet.phi());
      registry.fill(HIST("hJetArea"), jet.area());

      float dPhi{RecoDecay::constrainAngle(jet.phi() - psi.psi2, -o2::constants::math::PI)};
      registry.fill(HIST("hCentJetPtdPhiq2"), centrality, jet.pt() - (collision.rho() * jet.area()), dPhi, qPerc[0]);

      if (cfgrhoPhi) {
        auto rhoLocal = evalRho(rhoFit.get(), jetR, jet.phi(), collision.rho());
        registry.fill(HIST("hRhoPhi"), centrality, rhoLocal);
        registry.fill(HIST("hCentJetPtdPhiq2RhoPhi"), centrality, jet.pt() - (rhoLocal * jet.area()), dPhi, qPerc[0]);
        registry.fill(HIST("hCentPhi"), centrality, rhoFit->Eval(jet.phi()));
        registry.fill(HIST("hdPhiRhoPhi"), dPhi, rhoLocal);
      }
    }
    registry.fill(HIST("hEventCounter"), counter++);

    if (centrality < 30 || centrality > 50)
      return;
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process ese collisions", true);

  void processESEEPData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors, aod::QPercentileFT0Cs>::iterator const& collision,
                        soa::Filtered<aod::ChargedJets> const&,
                        aod::JetTracks const&)
  {

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;

    if (cfgEvSelOccupancy && !isOccupancyAccepted(collision))
      return;

    [[maybe_unused]] const auto psi{procEP<PsiFillerEP>(collision)};
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEEPData, "process ese collisions for filling EP and EPR", false);

  void processESEBackground(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors, aod::QPercentileFT0Cs>>::iterator const& collision,
                            soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                            aod::JetTracks const& tracks)
  {
    bkgFluctuationsRandomCone(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEBackground, "process random cones with event plane and ESE", false);

  void processESEOccupancy(soa::Join<aod::JetCollisions, aod::Qvectors, aod::QPercentileFT0Cs>::iterator const& collision,
                           soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks)
  {
    float count{0.5};
    registry.fill(HIST("hEventCounterOcc"), count++);
    const auto psi{procEP<PsiFillerOcc>(collision)};
    const auto qPerc{collision.qPERCFT0C()};

    auto occupancy{collision.trackOccupancyInTimeRange()};
    registry.fill(HIST("hPsiOccupancy"), collision.centFT0M(), psi.psi2, occupancy);
    registry.fill(HIST("hOccupancy"), occupancy);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;
    registry.fill(HIST("hEventCounterOcc"), count++);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;

      registry.fill(HIST("hTrackPt"), collision.centFT0M(), track.pt(), qPerc[0], occupancy);
      if (track.pt() < cfgOccupancyPtCut->at(0) || track.pt() > cfgOccupancyPtCut->at(1))
        continue;
      registry.fill(HIST("hTrackEta"), collision.centFT0M(), track.eta(), occupancy);
      registry.fill(HIST("hTrackPhi"), collision.centFT0M(), track.phi(), occupancy);
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEOccupancy, "process occupancy", false);

  void processMCParticleLevel(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator const& mcCollision,
                              soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                              soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetEventWeights, aod::ChargedMCParticleLevelJetConstituents>> const& jets,
                              aod::JetParticles const&)
  {
    float counter{0.5f};
    registry.fill(HIST("mcp/hEventCounter"), counter++);
    if (mcCollision.size() < 1) {
      return;
    }
    registry.fill(HIST("mcp/hEventCounter"), counter++);
    if (collisions.size() != 1) {
      return;
    }

    registry.fill(HIST("mcp/hEventCounter"), counter++);
    auto centrality{-1};
    bool fOccupancy = true;
    bool eventSel = true;
    for (const auto& col : collisions) {
      if (cfgisPbPb)
        centrality = col.centFT0M();
      if (cfgEvSelOccupancy && !isOccupancyAccepted(col))
        fOccupancy = false;
      if (!jetderiveddatautilities::selectCollision(col, eventSelectionBits))
        eventSel = false;
    }
    if (cfgEvSelOccupancy && !fOccupancy)
      return;
    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    if (!eventSel)
      return;

    registry.fill(HIST("mcp/hEventCounter"), counter++);

    registry.fill(HIST("mcp/hCentralitySel"), centrality);
    jetLoopMCP<aod::JetParticles>(jets, centrality, mcCollision.rho());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleLevel, "jets on particle level MC", false);

  void processMCDetectorLevel(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
                              ChargedMCDJets const& mcdjets,
                              aod::JetTracks const&,
                              aod::JetParticles const&)
  {
    float counter{0.5f};
    registry.fill(HIST("mcd/hEventCounter"), counter++);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;
    registry.fill(HIST("mcd/hEventCounter"), counter++);
    if (cfgEvSelOccupancy && !isOccupancyAccepted(collision))
      return;
    registry.fill(HIST("mcd/hEventCounter"), counter++);

    if (!(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }

    auto centrality = cfgisPbPb ? collision.centFT0M() : -1;

    registry.fill(HIST("mcd/hCentralitySel"), centrality);
    jetLoopMCD<aod::JetTracks>(mcdjets, centrality, collision.rho());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCDetectorLevel, "jets on detector level", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator const& mcCol,
                               soa::SmallGroups<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>> const& collisions,
                               ChargedMCDJets const& mcdjets,
                               JetMCPTable const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    float counter{0.5f};
    registry.fill(HIST("mcm/hMCEventCounter"), counter++);
    if (mcCol.size() < 1) {
      return;
    }
    if (collisions.size() != 1) {
      return;
    }
    registry.fill(HIST("mcm/hMCEventCounter"), counter++);
    if (!(std::abs(mcCol.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("mcm/hMCEventCounter"), counter++);

    for (const auto& collision : collisions) {
      float secCount{0.5f};
      registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);
      if (!(std::abs(collision.posZ()) < vertexZCut)) {
        return;
      }
      registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);

      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
        return;
      registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);

      if (cfgEvSelOccupancy && !isOccupancyAccepted(collision))
        return;
      registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);

      auto centrality = cfgisPbPb ? collision.centFT0M() : -1;
      if (cfgisPbPb)
        if (centrality < centRange->at(0) || centrality > centRange->at(1))
          registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);

      registry.fill(HIST("mcm/hCentralityAnalyzed"), centrality);
      matchedJetLoop<JetMCPTable, aod::JetTracks>(mcdjets.sliceBy(mcdjetsPerJCollision, collision.globalIndex()), centrality, collision.rho(), mcCol.rho());

      registry.fill(HIST("mcm/hMCDMatchedEventCounter"), secCount++);
    }
    registry.fill(HIST("mcm/hMCEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCChargedMatched, "jet MC process: geometrically matched MCP and MCD for response matrix and efficiency", false);

  // template <bool FillAllPsi, bool FillHist, typename EPCol>
  template <EventPlaneFiller P, typename EPCol>
  EventPlane procEP(EPCol const& vec)
  {
    constexpr std::array<float, 2> AmpCut{1e-8, 0.0};
    auto computeEP = [&AmpCut](const std::vector<float>& vec, auto det, float n) { return vec[2] > AmpCut[det] ? (1.0 / n) * std::atan2(vec[1], vec[0]) : 999.; };
    std::map<std::string, float> epMap;
    std::map<std::string, float> ep3Map;
    auto vec1{qVecNoESE<DetID::FT0A, P.hist>(vec)};
    epMap["FT0A"] = computeEP(vec1, 0, 2.0);
    ep3Map["FT0A"] = computeEP(vec1, 0, 3.0);
    if constexpr (P.psi) {
      auto vec2{qVecNoESE<DetID::FT0C, false>(vec)};
      epMap["FT0C"] = computeEP(vec2, 0, 2.0);
      ep3Map["FT0C"] = computeEP(vec2, 0, 3.0);
      epMap["FV0A"] = computeEP(qVecNoESE<DetID::FV0A, false>(vec), 0, 2.0);
      epMap["TPCpos"] = computeEP(qVecNoESE<DetID::TPCpos, false>(vec), 1, 2.0);
      epMap["TPCneg"] = computeEP(qVecNoESE<DetID::TPCneg, false>(vec), 1, 2.0);
      if constexpr (P.hist)
        fillEP(/*std::make_index_sequence<5>{},*/ vec, epMap);
      auto cosPsi = [](float psiX, float psiY) { return (static_cast<double>(psiX) == 999. || static_cast<double>(psiY) == 999.) ? 999. : std::cos(2.0 * (psiX - psiY)); };
      std::array<float, 3> epCorrContainer{};
      epCorrContainer[0] = cosPsi(epMap.at(cfgEPRefA), epMap.at(cfgEPRefC));
      epCorrContainer[1] = cosPsi(epMap.at(cfgEPRefA), epMap.at(cfgEPRefB));
      epCorrContainer[2] = cosPsi(epMap.at(cfgEPRefB), epMap.at(cfgEPRefC));
      if constexpr (P.hist)
        fillEPCos(/*std::make_index_sequence<3>{},*/ vec, epCorrContainer);
    }
    EventPlane localPlane;
    localPlane.psi2 = epMap.at(cfgEPRefA);
    localPlane.psi3 = ep3Map.at(cfgEPRefA);
    return localPlane;
    // return epMap.at(cfgEPRefA);
  }
  template </*std::size_t... Idx,*/ typename collision>
  void fillEPCos(/*const std::index_sequence<Idx...>&,*/ const collision& col, const std::array<float, 3>& Corr)
  {
    // static constexpr std::string CosList[] = {"hCosPsi2AmC", "hCosPsi2AmB", "hCosPsi2BmC"};
    // (registry.fill(HIST(CosList[Idx]), col.centrality(), Corr[Idx], col.qPERCFT0C()[0]), ...);
    registry.fill(HIST("hCosPsi2AmC"), col.centFT0M(), Corr[0], col.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2AmB"), col.centFT0M(), Corr[1], col.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2BmC"), col.centFT0M(), Corr[2], col.qPERCFT0C()[0]);
  }

  template </*std::size_t... Idx,*/ typename collision>
  void fillEP(/*const std::index_sequence<Idx...>&,*/ const collision& col, const std::map<std::string, float>& epMap)
  {
    // static constexpr std::string_view EpList[] = {"hPsi2FT0A", "hPsi2FV0A", "hPsi2FT0C", "hPsi2TPCpos", "hPsi2TPCneg"};
    // (registry.fill(HIST(EpList[Idx]), col.centrality(), epMap.at(std::string(RemovePrefix(EpList[Idx])))), ...);
    registry.fill(HIST("hPsi2FT0A"), col.centFT0M(), epMap.at("FT0A"));
    registry.fill(HIST("hPsi2FV0A"), col.centFT0M(), epMap.at("FV0A"));
    registry.fill(HIST("hPsi2FT0C"), col.centFT0M(), epMap.at("FT0C"));
    registry.fill(HIST("hPsi2TPCpos"), col.centFT0M(), epMap.at("TPCpos"));
    registry.fill(HIST("hPsi2TPCneg"), col.centFT0M(), epMap.at("TPCneg"));
  }
  constexpr std::string_view RemovePrefix(std::string_view str)
  {
    constexpr std::string_view Prefix{"hPsi2"};
    return str.substr(Prefix.size());
  }

  constexpr int detIDN(const DetID id)
  {
    switch (id) {
      case DetID::FT0C:
        return 0;
      case DetID::FT0A:
        return 1;
      case DetID::FT0M:
        return 2;
      case DetID::FV0A:
        return 3;
      case DetID::TPCpos:
        return 4;
      case DetID::TPCneg:
        return 5;
      case DetID::TPCall:
        return 6;
    }
    return -1;
  }

  template <DetID id, bool fill, typename Col>
  std::vector<float> qVecNoESE(Col collision)
  {
    const int nmode{2};
    int detId{detIDN(id)};
    int detInd{detId * 4 + cfgnTotalSystem * 4 * (nmode - 2)};
    if constexpr (fill) {
      if (collision.qvecAmp()[detInd] > 1e-8) {
        registry.fill(HIST("hQvecUncorV2"), collision.centFT0M(), collision.qvecRe()[detInd], collision.qvecIm()[detInd]);
        registry.fill(HIST("hQvecRectrV2"), collision.centFT0M(), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1]);
        registry.fill(HIST("hQvecTwistV2"), collision.centFT0M(), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2]);
        registry.fill(HIST("hQvecFinalV2"), collision.centFT0M(), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3]);
        registry.fill(HIST("hEPUncorV2"), collision.centFT0M(), 0.5 * std::atan2(collision.qvecIm()[detInd], collision.qvecRe()[detInd]));
        registry.fill(HIST("hEPRectrV2"), collision.centFT0M(), 0.5 * std::atan2(collision.qvecIm()[detInd + 1], collision.qvecRe()[detInd + 1]));
        registry.fill(HIST("hEPTwistV2"), collision.centFT0M(), 0.5 * std::atan2(collision.qvecIm()[detInd + 2], collision.qvecRe()[detInd + 2]));
      }
    }
    std::vector<float> qVec{};
    qVec.push_back(collision.qvecRe()[detInd + cfgnCorrLevel]);
    qVec.push_back(collision.qvecIm()[detInd + cfgnCorrLevel]);
    qVec.push_back(collision.qvecAmp()[detId]);
    return qVec;
  }

  template <typename col>
  bool isOccupancyAccepted(const col& collision)
  {
    auto occupancy{collision.trackOccupancyInTimeRange()};
    if (occupancy < cfgCutOccupancy->at(0) || occupancy > cfgCutOccupancy->at(1))
      return false;
    else
      return true;
  }

  template <typename Cent>
  bool isCentralitySelected(const Cent& centrality)
  {
    if (centrality < centRange->at(0) || centrality > centRange->at(1))
      return false;
    else
      return true;
  }

  template <bool fillHist, typename Col, typename TTracks, typename Jets>
  std::unique_ptr<TF1> fitRho(const Col& col, const EventPlane& ep, TTracks const& tracks, Jets const& jets)
  {
    float leadingJetPt = 0.0;
    float leadingJetEta = 0.0;
    // float leadingJetPhi = 0.0;
    for (const auto& jet : jets) {
      if (jet.pt() > leadingJetPt) {
        leadingJetPt = jet.pt();
        leadingJetEta = jet.eta();
        // leadingJetPhi = jet.phi();
      }
    }

    int nTrk{0};
    if (jets.size() > 0) {
      for (const auto& track : tracks) {
        if constexpr (fillHist)
          registry.fill(HIST("hTrackCounter"), 0.5);
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetR) && track.pt() >= 0.2 && track.pt() <= 5) {
          nTrk++;
          if constexpr (fillHist)
            registry.fill(HIST("hTrackCounter"), 1.5);
        }
      }
    }
    if (nTrk < 1)
      return nullptr;

    auto hPhiPt = std::unique_ptr<TH1F>(new TH1F("h_ptsum_sumpt_fit", "h_ptsum_sumpt fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI));
    for (const auto& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetR) && track.pt() >= 0.2 && track.pt() <= 5) {
        hPhiPt->Fill(track.phi(), track.pt());
        if constexpr (fillHist) {
          registry.fill(HIST("hTrackCounter"), 2.5);
          registry.fill(HIST("hPhiPtsum"), track.phi(), track.pt());
        }
      }
    }
    auto modulationFit = std::unique_ptr<TF1>(new TF1("fit_rholoc", "[0] * (1. + 2. * ([1] * std::cos(2. * (x - [2])) + [3] * std::cos(3. * (x - [4]))))", 0, o2::constants::math::TwoPI));

    modulationFit->SetParameter(0, 1.0);
    modulationFit->SetParameter(1, 0.01);
    modulationFit->SetParameter(3, 0.01);

    modulationFit->FixParameter(2, (ep.psi2 < 0) ? RecoDecay::constrainAngle(ep.psi2) : ep.psi2);
    modulationFit->FixParameter(4, (ep.psi3 < 0) ? RecoDecay::constrainAngle(ep.psi3) : ep.psi3);

    hPhiPt->Fit(modulationFit.get(), "Q", "", 0, o2::constants::math::TwoPI);

    if constexpr (fillHist) {
      registry.fill(HIST("hfitPar0"), col.centFT0M(), modulationFit->GetParameter(0));
      registry.fill(HIST("hfitPar1"), col.centFT0M(), modulationFit->GetParameter(1));
      registry.fill(HIST("hfitPar2"), col.centFT0M(), modulationFit->GetParameter(2));
      registry.fill(HIST("hfitPar3"), col.centFT0M(), modulationFit->GetParameter(3));
      registry.fill(HIST("hfitPar4"), col.centFT0M(), modulationFit->GetParameter(4));
    }

    if (modulationFit->GetParameter(0) <= 0)
      return nullptr;

    double chi2{0.};
    for (int i{0}; i < hPhiPt->GetXaxis()->GetNbins(); i++) {
      if (hPhiPt->GetBinContent(i + 1) <= 0.)
        continue;
      chi2 += std::pow((hPhiPt->GetBinContent(i + 1) - modulationFit->Eval(hPhiPt->GetXaxis()->GetBinCenter(1 + i))), 2) / hPhiPt->GetBinContent(i + 1);
    }

    int nDF{1};
    int numParams{2};
    nDF = static_cast<int>(modulationFit->GetXaxis()->GetNbins()) - numParams;
    if (nDF <= 0)
      return nullptr;

    auto cDF = 1. - TMath::Gamma(nDF, chi2);

    if constexpr (fillHist) {
      registry.fill(HIST("hPValueCentCDF"), col.centFT0M(), cDF);
      registry.fill(HIST("hCentChi2Ndf"), col.centFT0M(), chi2 / (static_cast<float>(nDF)));
    }

    return modulationFit;
  }

  template <typename phiT, typename rhoT>
  double evalRho(TF1* fit, const float& radius, phiT const& jetPhi, rhoT const& colRho)
  {
    double integralValue{fit->Integral(jetPhi - radius, jetPhi + radius)};
    double rhoLocal{colRho / (2 * radius * fit->GetParameter(0)) * integralValue};
    return rhoLocal;
  }

  template <typename TCollisions, typename TJets, typename TTracks>
  void bkgFluctuationsRandomCone(TCollisions const& collision, TJets const& jets, TTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (cfgEvSelOccupancy && !isOccupancyAccepted(collision))
      return;

    const auto psi{procEP<PsiFillerEse>(collision)};
    auto qPerc{collision.qPERCFT0C()};
    if (qPerc[0] < 0)
      return;

    TRandom3 randomNumber(0);
    float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
    float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, -o2::constants::math::PI);
        float dEta = track.eta() - randomConeEta;
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("hCentRhoRandomCone"), collision.centFT0M(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho());

    randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, -o2::constants::math::PI);
        float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("hCentRhoRandomConeRandomTrackDir"), collision.centFT0M(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho());

    if (jets.size() > 0) {
      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, -o2::constants::math::PI);
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while ((randomConeLeadJetDeltaR <= 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
        jetWasInCone = true;
        randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
        randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, -o2::constants::math::PI);
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePt = 0.0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, -o2::constants::math::PI);
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
      }
    }
    auto rho = collision.rho();
    std::unique_ptr<TF1> rhoFit{nullptr};
    if (cfgrhoPhi) {
      rhoFit = fitRho<false>(collision, psi, tracks, jets);
      if (!rhoFit)
        return;
      rho = evalRho(rhoFit.get(), randomConeR, randomConePhi, rho);
    }
    float dPhi{RecoDecay::constrainAngle(randomConePhi - psi.psi2, -o2::constants::math::PI)};
    registry.fill(HIST("hCentRhoRandomConewoLeadingJet"), collision.centFT0M(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rho, dPhi, qPerc[0]);

    double randomConePtWithoutOneLeadJet = 0;
    double randomConePtWithoutTwoLeadJet = 0;
    if (jets.size() > 1) { // if there are no jets, or just one, in the acceptance (from the jetfinder cuts) then one cannot find 2 leading jets
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, -o2::constants::math::PI);
          float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;
          if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            if (!isTrackInJet(track, jets.iteratorAt(0))) {
              randomConePtWithoutOneLeadJet += track.pt();
              if (!isTrackInJet(track, jets.iteratorAt(1))) {
                randomConePtWithoutTwoLeadJet += track.pt();
              }
            }
          }
        }
      }
    }
    registry.fill(HIST("hCentRhoRandomConeRndTrackDirwoOneLeadingJet"), collision.centFT0M(), randomConePtWithoutOneLeadJet - o2::constants::math::PI * randomConeR * randomConeR * rho, dPhi, qPerc[0]);
    registry.fill(HIST("hCentRhoRandomConeRndTrackDirwoTwoLeadingJet"), collision.centFT0M(), randomConePtWithoutTwoLeadJet - o2::constants::math::PI * randomConeR * randomConeR * rho, dPhi, qPerc[0]);
  }
  template <typename TTracks, typename TJets>
  bool isTrackInJet(TTracks const& track, TJets const& jet)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  // static constexpr std::string_view LevelJets[] = {"mcd/", "mcp/"};
  // enum JetType { MCP = 0,
  //                MCD = 1
  //               };
  // template <JetType jetLvl, typename Jets>
  template <typename JTracks, typename Jets>
  void jetLoopMCD(const Jets& jets, const float& centrality, const float& rho)
  {
    float weight = 1.0;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JTracks>(jet)) {
        continue;
      }
      auto pt = jet.pt();
      if (cfgbkgSubMC) {
        pt = jet.pt() - (rho * jet.area());
      }
      if (cfgUseMCEventWeights) {
        weight = jet.eventWeight();
      }
      registry.fill(/*HIST(LevelJets[jetLvl]) +*/ HIST("mcd/hJetSparse"), centrality, pt, jet.eta(), jet.phi(), weight); /* detector level mcm*/
    }
  }

  template <typename JTracks, typename Jets>
  void jetLoopMCP(const Jets& jets, const float& centrality, const float& rho)
  {
    bool mcLevelIsParticleLevel = true;
    float weight = 1.0;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JTracks>(jet, mcLevelIsParticleLevel)) {
        continue;
      }
      auto pt = jet.pt();
      if (cfgbkgSubMC) {
        pt = jet.pt() - (rho * jet.area());
      }
      if (cfgUseMCEventWeights) {
        weight = jet.eventWeight();
      }
      registry.fill(/*HIST(LevelJets[jetLvl]) +*/ HIST("mcp/hJetSparse"), centrality, pt, jet.eta(), jet.phi(), weight); /* detector level mcm*/
    }
  }

  template <typename MCPTab, typename JTracks, typename Jets>
  void matchedJetLoop(const Jets& jets, const float& centrality, const float& rho, const float& rho2)
  {
    float weight = 1.0;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JTracks>(jet)) {
        continue;
      }
      float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
      if (jet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }

      auto pt = jet.pt();
      if (cfgbkgSubMC) {
        pt = jet.pt() - (rho * jet.area());
      }
      if (cfgUseMCEventWeights) {
        weight = jet.eventWeight();
      }
      registry.fill(HIST("mcm/hJetSparse"), centrality, pt, jet.eta(), jet.phi(), weight); /* detector level mcm*/

      if (jet.has_matchedJetGeo()) {
        registry.fill(HIST("mcm/hDetSparseMatch"), centrality, pt, jet.eta(), jet.phi(), weight);
        for (const auto& matchedJet : jet.template matchedJetGeo_as<MCPTab>()) {
          if (matchedJet.pt() > pTHatMaxMCD * pTHat)
            continue;
          auto matchedpt = matchedJet.pt();
          if (cfgbkgSubMC) {
            matchedpt = matchedJet.pt() - (rho2 * matchedJet.area());
          }
          registry.fill(HIST("mcm/hPartSparseMatch"), centrality, matchedpt, matchedJet.eta(), matchedJet.phi(), weight);
          registry.fill(HIST("mcm/hMatchedJetsPtDelta"), matchedJet.pt(), jet.pt() - matchedJet.pt(), weight);
          registry.fill(HIST("mcm/hMatchedJetsPhiDelta"), matchedJet.phi(), jet.phi() - matchedJet.phi(), weight);
          registry.fill(HIST("mcm/hMatchedJetsEtaDelta"), matchedJet.eta(), jet.eta() - matchedJet.eta(), weight);
          registry.fill(HIST("mcm/hGenMatchedJetsPtDeltadPt"), centrality, matchedpt, (pt - matchedpt) / matchedpt, weight);
          registry.fill(HIST("mcm/hRecoMatchedJetsPtDeltadPt"), centrality, pt, (pt - matchedpt) / pt, weight);
          registry.fill(HIST("mcm/hRespMcDMcPMatch"), centrality, pt, matchedpt, weight);
        }
      }
    }
  }

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet, bool mcLevelIsParticleLevel = false)
  {
    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > -98.0);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double pt = constituent.pt();

        if ((!checkLeadConstituentMinPtForMcpJets && mcLevelIsParticleLevel) || (checkConstituentMinPt && pt >= leadingConstituentPtMin)) { // if the jet is mcp level and checkLeadConstituentMinPtForMcpJets is true, then the pt of the leading track of that jet does not need to be below the defined leadingConstituentPtMin cut
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }
    return true;
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc)}; }
