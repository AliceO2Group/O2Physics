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
/// \author Yubiao Wang <yubiao.wang@cern.ch>
/// \file jetChargedV2.cxx
/// \brief This file contains the implementation for the Charged Jet v2 analysis in the ALICE experiment

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <THn.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector2.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetChargedV2 {
  using JetBkgRhoMcCollisions = soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>;
  using ChargedMCDMatchedJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
  using ChargedMCPMatchedJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
  using ChargedMCDMatchedJetsWeighted = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>;
  using ChargedMCPMatchedJetsWeighted = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>;

  HistogramRegistry registry;
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<int> centralityMode{"centralityMode", 0, "0 = FT0C (default), 1 = FT0M, 2 = FT0A"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", 99, "additional cut on dcaZ to PV for tracks; uniformTracks in particular don't cut on this at all"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 1000., "maximum pT acceptance for tracks"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> jetPtMin{"jetPtMin", 0.15, "minimum pT acceptance for jets"};
  Configurable<float> jetPtMax{"jetPtMax", 200.0, "maximum pT acceptance for jets"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.9, "minimum eta acceptance for jets"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.9, "maximum eta acceptance for jets"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "number of bins for eta axes"};
  Configurable<float> jetRadius{"jetRadius", 0.2, "jet resolution parameters"};
  Configurable<float> randomConeLeadJetDeltaR{"randomConeLeadJetDeltaR", -99.0, "min distance between leading jet axis and random cone (RC) axis; if negative, min distance is set to automatic value of R_leadJet+R_RC "};

  Configurable<float> localRhoFitPtMin{"localRhoFitPtMin", 0.2, "Minimum track pT used for local rho fluctuation fit"};
  Configurable<float> localRhoFitPtMax{"localRhoFitPtMax", 5, "Maximum track pT used for local rho fluctuation fit"};

  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};

  //=====================< evt pln >=====================//
  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "event selection"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0M", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  ConfigurableAxis cfgAxisQvecF{"cfgAxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgAxisQvec{"cfgAxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {90, 0, 90}, ""};

  ConfigurableAxis cfgAxisVnCent{"cfgAxisVnCent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 50, 70, 100}, " % "};

  ConfigurableAxis cfgAxisEvtfit{"cfgAxisEvtfit", {10000, 0, 10000}, ""};
  EventPlaneHelper helperEP;
  int detId;
  int refAId;
  int refBId;

  //=====================< jetSpectraConfig to this analysis >=====================//
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<bool> checkLeadConstituentPtForMcpJets{"checkLeadConstituentPtForMcpJets", false, "flag to choose whether particle level jets should have their lead track pt above leadingConstituentPtMin to be accepted; off by default, as leadingConstituentPtMin cut is only applied on MCD jets for the Pb-Pb analysis using pp MC anchored to Pb-Pb for the response matrix"};
  Configurable<bool> cfgChkFitQuality{"cfgChkFitQuality", false, "check fit quality"};
  Configurable<bool> subtractMCPBackground{"subtractMCPBackground", true, "subtract MCP Background with General Purpose anchored MC"};
  Configurable<bool> useMedianRho{"useMedianRho", false, "use median rho for subtract MCP Background"};
  Configurable<bool> useLocalRho{"useLocalRho", false, "use local rho for subtract MCP Background"};

  float configSwitchLow = -98.0;
  float configSwitchHigh = 9998.0;
  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  template <typename T>
  int getDetId(const T& name)
  {
    if (name.value == "BPos" || name.value == "BNeg" || name.value == "BTot") {
      LOGF(warning, "Using deprecated label: %s. Please use TPCpos, TPCneg, TPCall instead.", name.value);
    }
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos" || name.value == "BPos") {
      return 4;
    } else if (name.value == "TPCneg" || name.value == "BNeg") {
      return 5;
    } else if (name.value == "TPCall" || name.value == "BTot") {
      return 6;
    } else {
      return 0;
    }
  }
  //=====================< evt p615ln | end >=====================//

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  double evtnum = 0;
  double accptTrack = 0;
  double fitTrack = 0;
  float collQvecAmpDetId = 1e-8;
  TH1F* hPtsumSumptFit = nullptr;
  TH1F* hPtsumSumptFitMCP = nullptr;
  TF1* fFitModulationV2v3 = 0x0;
  TF1* fFitModulationV2v3P = 0x0;
  TH1F* hPtsumSumptFitRM = nullptr;
  TF1* fFitModulationRM = 0x0;

  void init(o2::framework::InitContext&)
  {
    detId = getDetId(cfgDetName);
    refAId = getDetId(cfgRefAName);
    refBId = getDetId(cfgRefBName);
    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }
    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (std::abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }

    auto jetPtTemp = 0.0;
    jetPtBins.push_back(jetPtTemp);
    jetPtBinsRhoAreaSub.push_back(jetPtTemp);
    while (jetPtTemp < jetPtMax) {
      double jetPtTempA = 100.0;
      double jetPtTempB = 100.0;
      if (jetPtTemp < jetPtTempA) {
        jetPtTemp += 1.0;
        jetPtBins.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(-jetPtTemp);
      } else if (jetPtTemp < jetPtTempB) {
        jetPtTemp += 5.0;
        jetPtBins.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(-jetPtTemp);

      } else {
        jetPtTemp += 10.0;
        jetPtBins.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(-jetPtTemp);
      }
    }
    std::sort(jetPtBinsRhoAreaSub.begin(), jetPtBinsRhoAreaSub.end());

    //< MCAxis >//
    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec trackPtAxis = {200, -0.5, 199.5, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec trackEtaAxis = {nBinsEta, -1.0, 1.0, "#eta"};
    AxisSpec phiAxis = {160, -1.0, 7.0, "#varphi"};
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {nBinsEta, -1.0, 1.0, "#eta"};

    AxisSpec axisPt = {40, 0.0, 4.0};
    AxisSpec axisEta = {32, -0.8, 0.8};
    AxisSpec axixCent = {20, 0, 100};
    AxisSpec axisChID = {220, 0, 220};

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (doprocessCollisions || doprocessCollisionsWeighted) {
      registry.add("h_collisions", "number of events;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.get<TH1>(HIST("h_collisions"))->GetXaxis()->SetBinLabel(1, "allColl");
      registry.get<TH1>(HIST("h_collisions"))->GetXaxis()->SetBinLabel(2, "qualitySel");
      registry.get<TH1>(HIST("h_collisions"))->GetXaxis()->SetBinLabel(3, "centralitycut");
      registry.get<TH1>(HIST("h_collisions"))->GetXaxis()->SetBinLabel(4, "occupancycut");
      registry.add("h_collisions_weighted", "number of events;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.get<TH1>(HIST("h_collisions_weighted"))->GetXaxis()->SetBinLabel(1, "allColl");
      registry.get<TH1>(HIST("h_collisions_weighted"))->GetXaxis()->SetBinLabel(2, "qualitySel");
      registry.get<TH1>(HIST("h_collisions_weighted"))->GetXaxis()->SetBinLabel(3, "centralitycut");
      registry.get<TH1>(HIST("h_collisions_weighted"))->GetXaxis()->SetBinLabel(4, "occupancycut");
      registry.add("h_coll_phat", "collision #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_coll_phat_weighted", "collision #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collisions_Zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_fakecollisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}); // is not filled if running on data
      registry.add("h2_centrality_collisions", "event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 30000}}});
    }

    if (doprocessMCCollisions || doprocessMCCollisionsWeighted) {
      registry.add("h_mccollisions", "number of mc events; event status; entries", {HistType::kTH1F, {{10, 0.0, 10}}});
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(2, "noRecoColl");
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(3, "splitColl");
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(4, "recoEvtSel");
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(5, "centralitycut");
      registry.get<TH1>(HIST("h_mccollisions"))->GetXaxis()->SetBinLabel(6, "occupancycut");
      registry.add("h_mccollisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(2, "noRecoColl");
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(3, "splitColl");
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(4, "recoEvtSel");
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(5, "centralitycut");
      registry.get<TH1>(HIST("h_mccollisions_weighted"))->GetXaxis()->SetBinLabel(6, "occupancycut");
      registry.add("h_mccoll_phat", "mc collision #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_mccoll_phat_weighted", "mc collision #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_mccollisions_zvertex", "position of mc collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h2_centrality_mccollisions", "mc event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}});
    }

    // test plots in JJMC
    if (doprocessSigmaPtMCD) {
      registry.add("h2_centrality_jet_pt", "centrality vs. jet pT;centrality; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH2F, {centralityAxis, jetPtAxis}});
      registry.add("h2_centrality_jet_eta", "centrality vs. jet eta;centrality; #eta; counts", {HistType::kTH2F, {centralityAxis, jetEtaAxis}});
      registry.add("h2_centrality_jet_phi", "centrality vs. jet phi;centrality; #varphi; counts", {HistType::kTH2F, {centralityAxis, phiAxis}});
      registry.add("h2_jet_pt_jet_area", "jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_jet_pt_jet_ntracks", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_track_pt", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}});
      registry.add("h3_jet_pt_jet_eta_jet_phi", "jet pt vs. eta vs. phi", {HistType::kTH3F, {jetPtAxis, jetEtaAxis, phiAxis}});
    }
    // test end plots in JJMC

    if (doprocessInOutJetV2 || doprocessInOutJetV2MCD || doprocessSigmaPtData || doprocessSigmaPtMCD || doprocessSigmaPtAreaSubMCD) {
      // test
      registry.add("h_mcd_pt_before_matching_mcdprocess", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      //=====================< evt pln plot >=====================//
      AxisSpec axisCent{cfgAxisCent, "centrality"};
      AxisSpec axisQvec{cfgAxisQvec, "Q"};
      AxisSpec axisQvecF{cfgAxisQvecF, "Q"};
      AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};
      for (uint i = 0; i < cfgnMods->size(); i++) {
        histosQA.add(Form("histQvecUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
        histosQA.add(Form("histQvecRectrV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
        histosQA.add(Form("histQvecTwistV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
        histosQA.add(Form("histQvecFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

        histosQA.add(Form("histEvtPlUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRectrV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlTwistV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

        histosQA.add(Form("histEvtPlRes_SigRefAV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_SigRefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_RefARefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      }
      histosQA.add("histCent", "Centrality TrkProcess", HistType::kTH1F, {axisCent});

      //< fit quality >//
      registry.add("h_PvalueCDF_CombinFit", "cDF #chi^{2}; entries", {HistType::kTH1F, {{50, 0, 1}}});
      registry.add("h2_PvalueCDFCent_CombinFit", "p-value cDF vs centrality; centrality; p-value", {HistType::kTH2F, {{100, 0, 100}, {40, 0, 1}}});
      registry.add("h2_Chi2Cent_CombinFit", "Chi2 vs centrality; centrality; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 5}}});
      registry.add("h2_PChi2_CombinFit", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});

      registry.add("h2_PChi2_CombinFitA", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});
      registry.add("h2_PChi2_CombinFitB", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});

      registry.add("h_evtnum_NTrk", "eventNumber vs Number of Track ; #eventNumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

      registry.add("h_v2obs_centrality", "fitparameter v2obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});
      registry.add("h_v3obs_centrality", "fitparameter v3obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});
      registry.add("h_fitparaRho_evtnum", "fitparameter #rho_{0} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_fitparaPsi2_evtnum", "fitparameter #Psi_{2} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_fitparaPsi3_evtnum", "fitparameter #Psi_{3} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_evtnum_centrlity", "eventNumber vs centrality ; #eventNumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

      registry.add("h2_phi_rholocal", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho(#varphi) ", {HistType::kTH2F, {{40, 0., o2::constants::math::TwoPI}, {210, -10.0, 200.0}}});
      registry.add("h2_rholocal_cent", "#varphi vs #rho(#varphi); #cent;  #rho(#varphi) ", {HistType::kTH2F, {{100, 0., 100}, {210, -10.0, 200.0}}});
      //< \sigma p_T at local rho test plot | end >

      registry.add("h_jet_pt_rhoareasubtracted", "jet pT rhoareasubtracted;#it{p}_{T,jet} (GeV/#it{c}); entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      registry.add("leadJetPt", "leadJet Pt ", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("leadJetPhi", "leadJet constituent #phi ", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("leadJetEta", "leadJet constituent #eta ", {HistType::kTH1F, {{100, -1.0, 1.0}}});

      //< RC test plots >//
      registry.add("h3_centrality_deltapT_RandomCornPhi_rhorandomconewithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
      registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphi", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});

      registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}(#varphi); #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
      //< bkg sub plot | end >//
      //< median rho >//
      registry.add("h_jet_pt_in_plane_v2", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_out_of_plane_v2", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_in_plane_v3", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_out_of_plane_v3", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      registry.add("h2_centrality_jet_pt_in_plane_v2", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_pt_out_of_plane_v2", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      //< rho(phi) >//
      registry.add("h_jet_pt_inclusive_v2_rho", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_in_plane_v2_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_out_of_plane_v2_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_in_plane_v3_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_pt_out_of_plane_v3_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      registry.add("h2_centrality_jet_pt_in_plane_v2_rho", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_pt_out_of_plane_v2_rho", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});

      registry.add("h2_centrality_jet_pt_rhoareasubtracted", "centrality vs. jet pT;centrality; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH2F, {centralityAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_eta_rhoareasubtracted", "centrality vs. jet eta;centrality; #eta; counts", {HistType::kTH2F, {centralityAxis, jetEtaAxis}});
      registry.add("h2_centrality_jet_phi_rhoareasubtracted", "centrality vs. jet phi;centrality; #varphi; counts", {HistType::kTH2F, {centralityAxis, phiAxis}});
      registry.add("h2_jet_pt_track_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, trackPtAxis}});
    }

    if (doprocessSigmaPtMCP || doprocessSigmaPtAreaSubMCP) {
      registry.add("h_jet_pt_part", "partvjet pT;#it{p}_{T,jet}^{part} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_eta_part", "part jet #eta;#eta^{part}; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi_part", "part jet #varphi;#phi^{part}; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h2_jet_pt_part_jet_area_part", "part jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet}^{part} (GeV/#it{c}); Area_{jet}^{part}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_jet_pt_part_jet_ntracks_part", "part jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet, tracks}^{part}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_jet_eta_jet_phi_part", "part jet pt vs. eta vs. phi", {HistType::kTH3F, {jetPtAxis, jetEtaAxis, phiAxis}});

      registry.add("h_jet_pt_part_rhoareasubtracted", "part jet corr pT;#it{p}_{T,jet}^{part} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_eta_part_rhoareasubtracted", "part jet #eta;#eta^{part}; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi_part_rhoareasubtracted", "part jet #varphi;#varphi^{part}; counts", {HistType::kTH1F, {phiAxis}});

      registry.add("leadJetPtMCP", "MCP leadJet Pt ", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("leadJetPhiMCP", "MCP leadJet constituent #phi ", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("leadJetEtaMCP", "MCP leadJet constituent #eta ", {HistType::kTH1F, {{100, -1.0, 1.0}}});

      registry.add("h_mcp_v2obs_centrality", "fitparameter v2obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});
      registry.add("h_mcp_v3obs_centrality", "fitparameter v3obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});
      registry.add("h_mcp_fitparaRho_evtnum", "fitparameter #rho_{0} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_mcp_fitparaPsi2_evtnum", "fitparameter #Psi_{2} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_mcp_fitparaPsi3_evtnum", "fitparameter #Psi_{3} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

      registry.add("h2_mcp_phi_rholocal", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho(#varphi) ", {HistType::kTH2F, {{40, 0., o2::constants::math::TwoPI}, {210, -10.0, 200.0}}});
      registry.add("h2_mcp_centrality_rholocal", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho(#varphi) ", {HistType::kTH2F, {{120, -10.0, 110.0}, {210, -10.0, 200.0}}});
      registry.add("h_mcp_jet_pt_rhoareasubtracted", "jet pT rholocal;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mcp_jet_pt_rholocal", "jet pT rholocal;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      //< MCP fit test >//
      registry.add("h_mcp_evtnum_centrlity", "eventNumber vs centrality ; #eventNumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_ep2_evtnum", "fitparameter #rho_{0} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
      registry.add("h_ep3_evtnum", "fitparameter #Psi_{2} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

      registry.add("h_mcp_jet_pt_in_plane_v2_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mcp_jet_pt_out_of_plane_v2_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mcp_jet_pt_in_plane_v3_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mcp_jet_pt_out_of_plane_v3_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      registry.add("h2_mcp_centrality_jet_pt_in_plane_v2_rho", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      registry.add("h2_mcp_centrality_jet_pt_out_of_plane_v2_rho", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      registry.add("h2_mcp_centrality_jet_pt_in_plane_v3_rho", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
      registry.add("h2_mcp_centrality_jet_pt_out_of_plane_v3_rho", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});

      registry.add("h3_mcp_centrality_deltapT_RandomCornPhi_localrhovsphi", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});

      registry.add("h_mcColl_counts", " number of mc events; event status; entries", {HistType::kTH1F, {{10, 0, 10}}});
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(4, "splitColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(5, "recoEvtSel");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(6, "centralitycut");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(7, "occupancycut");
      registry.add("h_mcColl_rho", "mc collision rho;#rho (GeV/#it{c}); counts", {HistType::kTH1F, {{500, 0.0, 500.0}}});
      registry.add("h_mc_zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
    }
    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
    }
    if (doprocessJetsMatchedSubtracted) {
      registry.add("h_mc_collisions_matched", "mc collisions status;event status;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
      registry.add("h_mcd_events_matched", "mcd event status;event status;entries", {HistType::kTH1F, {{6, 0.0, 6.0}}});
      registry.add("h_mc_rho_matched", "mc collision rho;#rho (GeV/#it{c}); counts", {HistType::kTH1F, {{500, -100.0, 500.0}}});
      registry.add("h_accept_Track_Match", "all and accept track;Track;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("h_mcd_pt_before_matching", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mcd_pt_after_matching", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

      if (useMedianRho) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});

        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_mcdetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_mcpetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_in", "jet mcp corr pT vs. corr delta pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_in", "jet mcd corr pT vs. corr delta pT / jet mcd corr pt in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_in", "jet mcp corr pT vs. jet mcd corr pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});

        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_mcdetaconstraint", "corr pT mcd vs. corr cpT mcp out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_mcpetaconstraint", "corr pT mcd vs. corr cpT mcp out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_out", "jet mcp corr pT vs. corr delta pT / jet mcp corr pt out-of-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_out", "jet mcd corr pT vs. corr delta pT / jet mcd corr pt out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_out", "jet mcp corr pT vs. jet mcd corr pT / jet mcp corr pt out-of-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
      }
      if (useLocalRho) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_incl_rhoareasubtracted_mcdetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_incl_rhoareasubtracted_mcpetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_incl_rhoareasubtracted", "jet mcp corr pT vs. corr delta pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_incl_rhoareasubtracted", "jet mcd corr pT vs. corr delta pT / jet mcd corr pt in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_incl_rhoareasubtracted", "jet mcp corr pT vs. jet mcd corr pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});

        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_rhoareasubtracted_mcdetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_rhoareasubtracted_mcpetaconstraint", "corr pT mcd vs. corr cpT mcp in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_in_rhoareasubtracted", "jet mcp corr pT vs. corr delta pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_in_rhoareasubtracted", "jet mcd corr pT vs. corr delta pT / jet mcd corr pt in-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_in_rhoareasubtracted", "jet mcp corr pT vs. jet mcd corr pT / jet mcp corr pt in-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});

        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_rhoareasubtracted_mcdetaconstraint", "corr pT mcd vs. corr cpT mcp out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_rhoareasubtracted_mcpetaconstraint", "corr pT mcd vs. corr cpT mcp out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_out_rhoareasubtracted", "jet mcp corr pT vs. corr delta pT / jet mcp corr pt out-of-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_out_rhoareasubtracted", "jet mcd corr pT vs. corr delta pT / jet mcd corr pt out-of-plane;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_out_rhoareasubtracted", "jet mcp corr pT vs. jet mcd corr pT / jet mcp corr pt out-of-plane;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
      }
    }

    registry.add("h_accept_Track", "all and accept track;Track;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(1, "acceptTrk");
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(2, "acceptTrkInFit");
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(3, "beforeSumptFit");
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(4, "afterSumptFit");
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(5, "getNtrk");
    registry.get<TH1>(HIST("h_accept_Track"))->GetXaxis()->SetBinLabel(6, "getNtrkMCP");

    registry.add("h_track_pt", "track #it{p}_{T} ; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH1F, {trackPtAxis}});
    registry.add("h2_track_eta_track_phi", "track eta vs. track phi; #eta; #phi; counts", {HistType::kTH2F, {trackEtaAxis, phiAxis}});
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  Filter mcEventCuts = (nabs(aod::jmccollision::posZ) < vertexZCut);

  Preslice<aod::JetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;
  Preslice<ChargedMCDMatchedJets> mcdjetsPerJCollision = o2::aod::jet::collisionId;
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet, bool mcLevelIsParticleLevel = false)
  {
    double jetAreaFractionMinAcc = -98.0;
    double leadingConstituentPtMinAcc = -98.0;
    double leadingConstituentPtMaxAcc = 9998.0;
    if (jetAreaFractionMin > jetAreaFractionMinAcc) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > leadingConstituentPtMinAcc);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < leadingConstituentPtMaxAcc);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }
    if (mcLevelIsParticleLevel && !checkLeadConstituentPtForMcpJets) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
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

  //< collisions accept >//
  template <typename TColl>
  bool applyCollisionCuts(TColl const& collision, bool fillHistograms = false, bool isWeighted = false, float eventWeight = 1.0)
  {
    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }

    if (fillHistograms) {
      registry.fill(HIST("h_collisions"), 0.5);
      registry.fill(HIST("h2_centrality_collisions"), centrality, 0.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_collisions"), 1.5);
      registry.fill(HIST("h2_centrality_collisions"), centrality, 1.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    }

    if (centrality < centralityMin || centralityMax < centrality) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_collisions"), 2.5);
      registry.fill(HIST("h2_centrality_collisions"), centrality, 2.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);
    }

    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_collisions"), 3.5);
      registry.fill(HIST("h2_centrality_collisions"), centrality, 3.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_collisions_weighted"), 3.5, eventWeight);
    }

    return true;
  }

  template <typename TMCColl, typename TCollisions>
  bool applyMCCollisionCuts(TMCColl const& mccollision, TCollisions const& collisions, bool fillHistograms = false, bool isWeighted = false, float eventWeight = 1.0)
  {
    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = mccollision.centFT0M();
        break;
      case 2:
        centrality = mccollision.centFT0A();
        break;
      default:
        centrality = mccollision.centFT0C();
        break;
    }

    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 0.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 0.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 0.5, eventWeight);
    }

    if (collisions.size() < 1) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 1.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 1.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 1.5, eventWeight);
    }

    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 2.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 2.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 2.5, eventWeight);
    }

    bool hasSel8Coll = false;
    bool centralityIsGood = false;
    bool occupancyIsGood = false;
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) {
        hasSel8Coll = true;
      }
      if ((trackOccupancyInTimeRangeMin < collisions.begin().trackOccupancyInTimeRange()) && (collisions.begin().trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
        occupancyIsGood = true;
      }

      if ((centralityMin < centrality) && (centrality < centralityMax)) {
        centralityIsGood = true;
      }
    } else {
      for (auto const& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
          hasSel8Coll = true;
        }
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
        }

        float centrality = -1.0;
        switch (centralityMode) {
          case 1:
            centrality = collision.centFT0M();
            break;
          case 2:
            centrality = collision.centFT0A();
            break;
          default:
            centrality = collision.centFT0C();
            break;
        }
        if ((centralityMin < centrality) && (centrality < centralityMax)) {
          centralityIsGood = true;
        }
      }
    }

    if (!hasSel8Coll) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 3.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 3.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 3.5, eventWeight);
    }

    if (!centralityIsGood) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 4.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 4.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 4.5, eventWeight);
    }

    if (!occupancyIsGood) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("h_mccollisions"), 5.5);
      registry.fill(HIST("h2_centrality_mccollisions"), centrality, 5.5, eventWeight);
      if (isWeighted)
        registry.fill(HIST("h_mccollisions_weighted"), 5.5, eventWeight);
    }

    return true;
  }

  template <typename TJets>
  void fillJetSpectraHistograms(TJets const& jet, float centrality, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h2_centrality_jet_pt"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_centrality_jet_eta"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_jet_area"), jet.pt(), jet.area(), weight);
      registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracksIds().size(), weight);
      registry.fill(HIST("h3_jet_pt_jet_eta_jet_phi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }

    for (const auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
    }
  }

  double chiSquareCDF(int nDF, double x)
  {
    return TMath::Gamma(nDF / 2., x / 2.);
  }

  // leading jet fill
  template <typename T>
  void fillLeadingJetQA(T const& jets, double& leadingJetPt, double& leadingJetPhi, double& leadingJetEta)
  {
    for (const auto& jet : jets) {
      if (jet.pt() > leadingJetPt) {
        leadingJetPt = jet.pt();
        leadingJetEta = jet.eta();
        leadingJetPhi = jet.phi();
      }
    }
    registry.fill(HIST("leadJetPt"), leadingJetPt);
    registry.fill(HIST("leadJetPhi"), leadingJetPhi);
    registry.fill(HIST("leadJetEta"), leadingJetEta);
  }

  // create h_ptsum_sumpt_fit, with number of Track
  template <typename T, typename U>
  void getNtrk(T const& tracks, U const& jets, int& nTrk, double& evtnum, double& leadingJetEta)
  {
    if (jets.size() > 0) {
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
          registry.fill(HIST("h_accept_Track"), 4.5);
          nTrk += 1;
        }
      }
      registry.fill(HIST("h_evtnum_NTrk"), evtnum, nTrk);
    }
  }

  // fill nTrk plot for fit rho(varphi)
  template <typename U, typename J>
  void fillNtrkCheck(U const& tracks, J const& jets, TH1F* hPtsumSumptFit, double& leadingJetEta)
  {
    if (jets.size() > 0) {
      for (auto const& trackfit : tracks) {
        registry.fill(HIST("h_accept_Track"), 0.5);
        if (jetderiveddatautilities::selectTrack(trackfit, trackSelection) && (std::fabs(trackfit.eta() - leadingJetEta) > jetRadius) && trackfit.pt() >= localRhoFitPtMin && trackfit.pt() <= localRhoFitPtMax) {
          registry.fill(HIST("h_accept_Track"), 1.5);
        }
      }

      for (auto const& track : tracks) {
        registry.fill(HIST("h_accept_Track"), 2.5);
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
          registry.fill(HIST("h_accept_Track"), 3.5);
          hPtsumSumptFit->Fill(track.phi(), track.pt());
        }
      }
    }
  }

  // MCP leading jet fill
  template <typename T>
  void fillLeadingJetQAMCP(T const& jets, double& leadingJetPt, double& leadingJetPhi, double& leadingJetEta)
  {
    for (const auto& jet : jets) {
      if (jet.pt() > leadingJetPt) {
        leadingJetPt = jet.pt();
        leadingJetEta = jet.eta();
        leadingJetPhi = jet.phi();
      }
    }
    registry.fill(HIST("leadJetPtMCP"), leadingJetPt);
    registry.fill(HIST("leadJetPhiMCP"), leadingJetPhi);
    registry.fill(HIST("leadJetEtaMCP"), leadingJetEta);
  }

  // Run General_Purpose MC MCP
  template <typename U, typename J>
  void fitFncAreaSubMCP(U const& collision, J const& jets, TH1F* hPtsumSumptFitMCP, bool mcLevelIsParticleLevel, float weight = 1.0)
  {
    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }

    double ep2 = 0.;
    double ep3 = 0.;
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      } else if (nmode == cfgNmodB) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      }
    }

    const char* fitFunctionV2v3P = "[0] * (1. + 2. * ([1] * std::cos(2. * (x - [2])) + [3] * std::cos(3. * (x - [4]))))";
    fFitModulationV2v3P = new TF1("fit_kV3", fitFunctionV2v3P, 0, o2::constants::math::TwoPI);
    //=========================< set parameter >=========================//
    fFitModulationV2v3P->SetParameter(0, 1.);
    fFitModulationV2v3P->SetParameter(1, 0.01);
    fFitModulationV2v3P->SetParameter(3, 0.01);

    double ep2fix = 0.;
    double ep3fix = 0.;

    if (ep2 < 0) {
      ep2fix = RecoDecay::constrainAngle(ep2);
      fFitModulationV2v3P->FixParameter(2, ep2fix);
    } else {
      fFitModulationV2v3P->FixParameter(2, ep2);
    }
    if (ep3 < 0) {
      ep3fix = RecoDecay::constrainAngle(ep3);
      fFitModulationV2v3P->FixParameter(4, ep3fix);
    } else {
      fFitModulationV2v3P->FixParameter(4, ep3);
    }

    hPtsumSumptFitMCP->Fit(fFitModulationV2v3P, "Q", "ep", 0, o2::constants::math::TwoPI);

    double temppara[5];
    temppara[0] = fFitModulationV2v3P->GetParameter(0);
    temppara[1] = fFitModulationV2v3P->GetParameter(1);
    temppara[2] = fFitModulationV2v3P->GetParameter(2);
    temppara[3] = fFitModulationV2v3P->GetParameter(3);
    temppara[4] = fFitModulationV2v3P->GetParameter(4);
    if (temppara[0] == 0) {
      return;
    }
    registry.fill(HIST("h_mcp_fitparaRho_evtnum"), evtnum, temppara[0]);
    registry.fill(HIST("h_mcp_fitparaPsi2_evtnum"), evtnum, temppara[2]);
    registry.fill(HIST("h_mcp_fitparaPsi3_evtnum"), evtnum, temppara[4]);
    registry.fill(HIST("h_mcp_v2obs_centrality"), centrality, temppara[1]);
    registry.fill(HIST("h_mcp_v3obs_centrality"), centrality, temppara[3]);
    registry.fill(HIST("h_mcp_evtnum_centrlity"), evtnum, centrality);

    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
          continue;
        }
        if (jet.r() != round(selectedJetsRadius * 100.0f)) {
          continue;
        }

        double integralValue = fFitModulationV2v3P->Integral(jet.phi() - jetRadius, jet.phi() + jetRadius);
        double rholocal = collision.rho() / (2 * jetRadius * temppara[0]) * integralValue;
        registry.fill(HIST("h2_mcp_phi_rholocal"), jet.phi() - ep2, rholocal, weight);
        registry.fill(HIST("h2_mcp_centrality_rholocal"), centrality, rholocal, weight);
        if (nmode == cfgNmodA) {
          registry.fill(HIST("h_mcp_jet_pt_rholocal"), jet.pt() - (rholocal * jet.area()), weight);

          double phiMinusPsi2;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          phiMinusPsi2 = jet.phi() - ep2;
          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_mcp_jet_pt_in_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_in_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), weight);
          } else {
            registry.fill(HIST("h_mcp_jet_pt_out_of_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_out_of_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), weight);
          }
        } else if (nmode == cfgNmodB) {
          double phiMinusPsi3;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_mcp_jet_pt_in_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_in_plane_v3_rho"), centrality, jet.pt() - (rholocal * jet.area()), weight);
          } else {
            registry.fill(HIST("h_mcp_jet_pt_out_of_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_out_of_plane_v3_rho"), centrality, jet.pt() - (rholocal * jet.area()), weight);
          }
        }
      }
    }
  }

  // Run jet-jet MC MCP leading jet fill
  template <typename U, typename J>
  void fitFncMCP(U const& collision, J const& jets, bool mcLevelIsParticleLevel, float weight = 1.0)
  {
    double ep2 = 0.;
    double ep3 = 0.;
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      } else if (nmode == cfgNmodB) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      }
    }
    registry.fill(HIST("h_mcp_evtnum_centrlity"), evtnum, collision.centFT0M());
    registry.fill(HIST("h_ep2_evtnum"), evtnum, ep2);
    registry.fill(HIST("h_ep3_evtnum"), evtnum, ep3);

    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
          continue;
        }
        if (jet.r() != round(selectedJetsRadius * 100.0f)) {
          continue;
        }

        if (nmode == cfgNmodA) {
          registry.fill(HIST("h_mcp_jet_pt_rhoareasubtracted"), jet.pt(), weight);

          double phiMinusPsi2;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          phiMinusPsi2 = jet.phi() - ep2;
          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_mcp_jet_pt_in_plane_v2_rho"), jet.pt(), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_in_plane_v2_rho"), collision.centFT0M(), jet.pt(), weight);
          } else {
            registry.fill(HIST("h_mcp_jet_pt_out_of_plane_v2_rho"), jet.pt(), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_out_of_plane_v2_rho"), collision.centFT0M(), jet.pt(), weight);
          }
        } else if (nmode == cfgNmodB) {
          double phiMinusPsi3;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_mcp_jet_pt_in_plane_v3_rho"), jet.pt(), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_in_plane_v3_rho"), collision.centFT0M(), jet.pt(), weight);
          } else {
            registry.fill(HIST("h_mcp_jet_pt_out_of_plane_v3_rho"), jet.pt(), weight);
            registry.fill(HIST("h2_mcp_centrality_jet_pt_out_of_plane_v3_rho"), collision.centFT0M(), jet.pt(), weight);
          }
        }
      }
    }
  }

  // Run General_Purpose MC MCP
  template <typename TJets>
  void fillMCPAreaSubHistograms(TJets const& jet, float rho = 0.0, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      // fill mcp jet histograms
      double jetcorrpt = jet.pt() - (rho * jet.area());
      registry.fill(HIST("h_jet_pt_part_rhoareasubtracted"), jetcorrpt, weight);
      if (jetcorrpt > 0) {
        registry.fill(HIST("h_jet_eta_part_rhoareasubtracted"), jet.eta(), weight);
        registry.fill(HIST("h_jet_phi_part_rhoareasubtracted"), jet.phi(), weight);
      }
    }
  }

  // Run jet-jet MC MCP
  template <typename TJets>
  void fillMCPHistograms(TJets const& jet, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      // fill mcp jet histograms
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h3_jet_pt_jet_eta_jet_phi_part"), jet.pt(), jet.eta(), jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_area_part"), jet.pt(), jet.area(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_ntracks_part"), jet.pt(), jet.tracksIds().size(), weight);
    }
  }

  template <typename TTracks>
  void fillTrackHistograms(TTracks const& track, float weight = 1.0)
  {
    registry.fill(HIST("h_track_pt"), track.pt(), weight);
    registry.fill(HIST("h2_track_eta_track_phi"), track.eta(), track.phi(), weight);
  }

  template <typename TBase, typename TTag>
  void fillGeoMatchedHistograms(TBase const& jetMCD, double ep2, float rho, float mcrho = 0.0, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetMCD.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jetMCD.has_matchedJetGeo()) {
      for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TTag>>()) {
        if (jetMCP.pt() > pTHatMaxMCD * pTHat) {
          continue;
        }
        if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
          int evtPlnAngleA = 7;
          int evtPlnAngleB = 3;
          int evtPlnAngleC = 5;
          double phiMinusPsi2 = jetMCD.phi() - ep2;
          double corrTagjetpt = jetMCP.pt() - (mcrho * jetMCP.area());
          double corrBasejetpt = jetMCD.pt() - (rho * jetMCD.area());
          double dcorrpt = corrTagjetpt - corrBasejetpt;
          if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_in"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            } else {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_out"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            }
          }
          if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_in"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_in"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            } else {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_out"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_out"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            }
          }
        }
      }
    }
  }

  template <typename TBase, typename TTag>
  void fillGeoMatchedCorrHistograms(TBase const& jetMCD, TF1* fFitModulationRM, float tempparaA, double ep2, float rho, bool subtractMCPBackground, float mcrho, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetMCD.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jetMCD.r() == round(selectedJetsRadius * 100.0f) && jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      double corrBasejetpt = jetMCD.pt() - rho * jetMCD.area();
      registry.fill(HIST("h_mcd_pt_before_matching"), corrBasejetpt, weight);
    }

    if (jetMCD.has_matchedJetGeo()) {
      for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TTag>>()) {
        if (jetMCP.pt() > pTHatMaxMCD * pTHat) {
          continue;
        }
        if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
          if (jetMCD.r() == round(selectedJetsRadius * 100.0f) && jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            double corrBasejetpt = jetMCD.pt() - rho * jetMCD.area();
            registry.fill(HIST("h_mcd_pt_after_matching"), corrBasejetpt, weight);
          }
        }
        if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
          int evtPlnAngleA = 7;
          int evtPlnAngleB = 3;
          int evtPlnAngleC = 5;
          double integralValue = fFitModulationRM->Integral(jetMCD.phi() - jetRadius, jetMCD.phi() + jetRadius);
          double rholocal = rho / (2 * jetRadius * tempparaA) * integralValue;
          double corrBasejetpt = jetMCD.pt() - (rholocal * jetMCD.area());
          double corrTagjetpt;
          if (subtractMCPBackground) {
            double integralValueMCP = fFitModulationRM->Integral(jetMCP.phi() - jetRadius, jetMCP.phi() + jetRadius);
            double rholocalMCP = mcrho / (2 * jetRadius * tempparaA) * integralValueMCP;
            corrTagjetpt = jetMCP.pt() - (rholocalMCP * jetMCP.area());
          } else {
            corrTagjetpt = jetMCP.pt() - (mcrho * jetMCP.area());
          }
          double dcorrpt = corrTagjetpt - corrBasejetpt;
          double phiMinusPsi2 = jetMCD.phi() - ep2;
          if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_incl_rhoareasubtracted_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_incl_rhoareasubtracted"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_rhoareasubtracted_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_in_rhoareasubtracted"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            } else {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_rhoareasubtracted_mcdetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_out_rhoareasubtracted"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
            }
          }
          if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_incl_rhoareasubtracted_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_incl_rhoareasubtracted"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_incl_rhoareasubtracted"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_in_rhoareasubtracted_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_in_rhoareasubtracted"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_in_rhoareasubtracted"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            } else {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_out_rhoareasubtracted_mcpetaconstraint"), corrBasejetpt, corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_out_rhoareasubtracted"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_out_rhoareasubtracted"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
            }
          }
        }
      }
    }
  }

  template <typename TBase, typename TTag>
  void fillMatchedHistograms(TBase const& jetMCD, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetMCD.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    // fill geometry matched histograms
    if (jetMCD.has_matchedJetGeo()) {
      for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TTag>>()) {
        if (jetMCP.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
          continue;
        }
        if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
          double dpt = jetMCD.pt() - jetMCP.pt();
          if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo"), jetMCD.pt(), dpt / jetMCD.pt(), weight);
          }
          if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo"), jetMCP.pt(), dpt / jetMCP.pt(), weight);
            registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo"), jetMCP.pt(), jetMCD.pt() / jetMCP.pt(), weight);
          }
        }
      }
    }
  }

  void processCollisions(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    bool fillHistograms = true;
    bool isWeighted = false;
    if (!applyCollisionCuts(collision, fillHistograms, isWeighted)) {
      return;
    }

    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ());
  }
  PROCESS_SWITCH(JetChargedV2, processCollisions, "collisions Data and MCD", true);

  void processCollisionsWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                                 aod::JetMcCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      registry.fill(HIST("h_fakecollisions"), 0.5);
    }
    bool fillHistograms = true;
    bool isWeighted = true;
    float eventWeight = collision.weight();
    if (!applyCollisionCuts(collision, fillHistograms, isWeighted, eventWeight)) {
      return;
    }

    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ(), eventWeight);

    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    registry.fill(HIST("h_coll_phat"), pTHat);
    registry.fill(HIST("h_coll_phat_weighted"), pTHat, eventWeight);
  }
  PROCESS_SWITCH(JetChargedV2, processCollisionsWeighted, "weighted collisions for MCD", false);

  void processMCCollisions(soa::Filtered<aod::JetMcCollisions>::iterator const& mccollision, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  {
    bool fillHistograms = true;
    bool isWeighted = false;
    if (!applyMCCollisionCuts(mccollision, collisions, fillHistograms, isWeighted)) {
      return;
    }

    registry.fill(HIST("h_mccollisions_zvertex"), mccollision.posZ());
  }
  PROCESS_SWITCH(JetChargedV2, processMCCollisions, "collisions MCP", false);

  void processMCCollisionsWeighted(soa::Filtered<aod::JetMcCollisions>::iterator const& mccollision, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  {
    bool fillHistograms = true;
    bool isWeighted = true;
    float eventWeight = mccollision.weight();
    if (!applyMCCollisionCuts(mccollision, collisions, fillHistograms, isWeighted, eventWeight)) {
      return;
    }

    registry.fill(HIST("h_mccollisions_zvertex"), mccollision.posZ(), eventWeight);

    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    registry.fill(HIST("h_mccoll_phat"), pTHat);
    registry.fill(HIST("h_mccoll_phat_weighted"), pTHat, eventWeight);
  }
  PROCESS_SWITCH(JetChargedV2, processMCCollisionsWeighted, "weighted collisions for MCP", false);

  void processInOutJetV2(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                         aod::JetTracks const&)
  {
    bool fillHistograms = true;
    bool isWeighted = false;
    if (!applyCollisionCuts(collision, fillHistograms, isWeighted)) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
    }

    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    histosQA.fill(HIST("histCent"), collision.cent());

    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      int refAInd = refAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      int refBInd = refBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId || collision.qvecAmp()[refAId] < collQvecAmpDetId || collision.qvecAmp()[refBId] < collQvecAmpDetId) {
          histosQA.fill(HIST("histQvecUncorV2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
          histosQA.fill(HIST("histQvecRectrV2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
          histosQA.fill(HIST("histQvecTwistV2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
          histosQA.fill(HIST("histQvecFinalV2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

          histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());

          histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
        }
      } else if (nmode == cfgNmodB) {
        histosQA.fill(HIST("histQvecUncorV3"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
        histosQA.fill(HIST("histQvecRectrV3"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
        histosQA.fill(HIST("histQvecTwistV3"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
        histosQA.fill(HIST("histQvecFinalV3"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

        histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRectrV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlTwistV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());

        histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
      }
    }

    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      if (nmode == cfgNmodA) {
        double phiMinusPsi2;
        if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
          continue;
        }
        float ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
        for (auto const& jet : jets) {
          if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (!isAcceptedJet<aod::JetTracks>(jet)) {
            continue;
          }
          if (jet.r() != round(selectedJetsRadius * 100.0f)) {
            continue;
          }
          registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jet.pt() - (collision.rho() * jet.area()), 1.0);

          phiMinusPsi2 = jet.phi() - ep2;
          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2"), centrality, jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2"), centrality, jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      } else if (nmode == cfgNmodB) {
        double phiMinusPsi3;
        float ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
        for (auto const& jet : jets) {
          if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (!isAcceptedJet<aod::JetTracks>(jet)) {
            continue;
          }
          if (jet.r() != round(selectedJetsRadius * 100.0f)) {
            continue;
          }
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetChargedV2, processInOutJetV2, "Jet V2 in and out of plane", false);

  void processInOutJetV2MCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                            soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                            aod::JetTracks const&)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
    }
    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    histosQA.fill(HIST("histCent"), collision.cent());
    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      int refAInd = refAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      int refBInd = refBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId || collision.qvecAmp()[refAId] < collQvecAmpDetId || collision.qvecAmp()[refBId] < collQvecAmpDetId) {
          histosQA.fill(HIST("histQvecUncorV2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
          histosQA.fill(HIST("histQvecRectrV2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
          histosQA.fill(HIST("histQvecTwistV2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
          histosQA.fill(HIST("histQvecFinalV2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

          histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());

          histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
        }
      } else if (nmode == cfgNmodB) {
        histosQA.fill(HIST("histQvecUncorV3"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
        histosQA.fill(HIST("histQvecRectrV3"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
        histosQA.fill(HIST("histQvecTwistV3"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
        histosQA.fill(HIST("histQvecFinalV3"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

        histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRectrV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlTwistV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());

        histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.cent());
      }
    }

    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      if (nmode == cfgNmodA) {
        double phiMinusPsi2;
        if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
          continue;
        }
        float ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
        for (auto const& jet : jets) {
          if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (!isAcceptedJet<aod::JetTracks>(jet)) {
            continue;
          }
          if (jet.r() != round(selectedJetsRadius * 100.0f)) {
            continue;
          }
          registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jet.pt() - (collision.rho() * jet.area()), 1.0);

          phiMinusPsi2 = jet.phi() - ep2;
          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2"), centrality, jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2"), centrality, jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      } else if (nmode == cfgNmodB) {
        double phiMinusPsi3;
        float ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
        for (auto const& jet : jets) {
          if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (!isAcceptedJet<aod::JetTracks>(jet)) {
            continue;
          }
          if (jet.r() != round(selectedJetsRadius * 100.0f)) {
            continue;
          }
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetChargedV2, processInOutJetV2MCD, "Jet V2 in and out of plane MCD", false);

  void processSigmaPtData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                          soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                          aod::JetTracks const& tracks)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }

    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
    }

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    fillLeadingJetQA(jets, leadingJetPt, leadingJetPhi, leadingJetEta);

    int nTrk = 0;
    getNtrk(tracks, jets, nTrk, evtnum, leadingJetEta);
    if (nTrk <= 0) {
      return;
    }
    hPtsumSumptFit = new TH1F("h_ptsum_sumpt_fit", "h_ptsum_sumpt fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI);

    fillNtrkCheck(tracks, jets, hPtsumSumptFit, leadingJetEta);

    double ep2 = 0.;
    double ep3 = 0.;
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      } else if (nmode == cfgNmodB) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      }
    }

    const char* fitFunctionV2v3 = "[0] * (1. + 2. * ([1] * std::cos(2. * (x - [2])) + [3] * std::cos(3. * (x - [4]))))";
    fFitModulationV2v3 = new TF1("fit_kV3", fitFunctionV2v3, 0, o2::constants::math::TwoPI);
    //=========================< set parameter >=========================//
    fFitModulationV2v3->SetParameter(0, 1.);
    fFitModulationV2v3->SetParameter(1, 0.01);
    fFitModulationV2v3->SetParameter(3, 0.01);

    double ep2fix = 0.;
    double ep3fix = 0.;

    if (ep2 < 0) {
      ep2fix = RecoDecay::constrainAngle(ep2);
      fFitModulationV2v3->FixParameter(2, ep2fix);
    } else {
      fFitModulationV2v3->FixParameter(2, ep2);
    }
    if (ep3 < 0) {
      ep3fix = RecoDecay::constrainAngle(ep3);
      fFitModulationV2v3->FixParameter(4, ep3fix);
    } else {
      fFitModulationV2v3->FixParameter(4, ep3);
    }

    hPtsumSumptFit->Fit(fFitModulationV2v3, "Q", "ep", 0, o2::constants::math::TwoPI);

    double temppara[5];
    temppara[0] = fFitModulationV2v3->GetParameter(0);
    temppara[1] = fFitModulationV2v3->GetParameter(1);
    temppara[2] = fFitModulationV2v3->GetParameter(2);
    temppara[3] = fFitModulationV2v3->GetParameter(3);
    temppara[4] = fFitModulationV2v3->GetParameter(4);

    registry.fill(HIST("h_fitparaRho_evtnum"), evtnum, temppara[0]);
    registry.fill(HIST("h_fitparaPsi2_evtnum"), evtnum, temppara[2]);
    registry.fill(HIST("h_fitparaPsi3_evtnum"), evtnum, temppara[4]);
    registry.fill(HIST("h_v2obs_centrality"), centrality, temppara[1]);
    registry.fill(HIST("h_v3obs_centrality"), centrality, temppara[3]);
    registry.fill(HIST("h_evtnum_centrlity"), evtnum, centrality);

    if (temppara[0] == 0) {
      return;
    }

    int nDF = 1;
    int numOfFreePara = 2;
    nDF = static_cast<int>(fFitModulationV2v3->GetXaxis()->GetNbins()) - numOfFreePara;
    if (nDF == 0 || static_cast<float>(nDF) <= 0.)
      return;
    double chi2 = 0.;
    for (int i = 0; i < hPtsumSumptFit->GetXaxis()->GetNbins(); i++) {
      if (hPtsumSumptFit->GetBinContent(i + 1) <= 0.)
        continue;
      chi2 += std::pow((hPtsumSumptFit->GetBinContent(i + 1) - fFitModulationV2v3->Eval(hPtsumSumptFit->GetXaxis()->GetBinCenter(1 + i))), 2) / hPtsumSumptFit->GetBinContent(i + 1);
    }

    double chiSqr = 999.;
    double cDF = 1.;

    chiSqr = chi2;
    cDF = 1. - chiSquareCDF(nDF, chiSqr);

    int evtCentAreaMin = 0;
    int evtCentAreaMax = 5;
    int evtMidAreaMin = 30;
    int evtMidAreaMax = 50;

    if (cfgChkFitQuality) {
      registry.fill(HIST("h_PvalueCDF_CombinFit"), cDF);
      registry.fill(HIST("h2_PvalueCDFCent_CombinFit"), centrality, cDF);
      registry.fill(HIST("h2_Chi2Cent_CombinFit"), centrality, chiSqr / (static_cast<float>(nDF)));
      registry.fill(HIST("h2_PChi2_CombinFit"), cDF, chiSqr / (static_cast<float>(nDF)));
      if (centrality >= evtCentAreaMin && centrality <= evtCentAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitA"), cDF, chiSqr / (static_cast<float>(nDF)));

      } else if (centrality >= evtMidAreaMin && centrality <= evtMidAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitB"), cDF, chiSqr / (static_cast<float>(nDF)));
      }
      if (centrality >= evtCentAreaMin && centrality <= evtCentAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitA"), cDF, chiSqr / (static_cast<float>(nDF)));

      } else if (centrality >= evtMidAreaMin && centrality <= evtMidAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitB"), cDF, chiSqr / (static_cast<float>(nDF)));
      }
    }

    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<aod::JetTracks>(jet)) {
          continue;
        }
        if (jet.r() != round(selectedJetsRadius * 100.0f)) {
          continue;
        }

        double integralValue = fFitModulationV2v3->Integral(jet.phi() - jetRadius, jet.phi() + jetRadius);
        double rholocal = collision.rho() / (2 * jetRadius * temppara[0]) * integralValue;
        registry.fill(HIST("h2_rholocal_cent"), centrality, rholocal, 1.0);

        if (nmode == cfgNmodA) {
          double phiMinusPsi2;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          phiMinusPsi2 = jet.phi() - ep2;

          registry.fill(HIST("h2_phi_rholocal"), jet.phi() - ep2, rholocal, 1.0);
          registry.fill(HIST("h_jet_pt_inclusive_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);

          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), 1.0);
          }
        } else if (nmode == cfgNmodB) {
          double phiMinusPsi3;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
          }
        }
      }
    }
    // RCpT
    for (uint i = 0; i < cfgnMods->size(); i++) {
      TRandom3 randomNumber(0);
      float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
      float randomConePt = 0;
      double integralValueRC = fFitModulationV2v3->Integral(randomConePhi - randomConeR, randomConePhi + randomConeR);
      double rholocalRC = collision.rho() / (2 * randomConeR * temppara[0]) * integralValueRC;

      int nmode = cfgnMods->at(i);
      if (nmode == cfgNmodA) {
        double rcPhiPsi2;
        rcPhiPsi2 = randomConePhi - ep2;

        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphi"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);

        // removing the leading jet from the random cone
        if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
          float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
          float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

          bool jetWasInCone = false;
          while ((randomConeLeadJetDeltaR <= 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
            jetWasInCone = true;
            randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
            randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
            dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
          }
          if (jetWasInCone) {
            randomConePt = 0.0;
            for (auto const& track : tracks) {
              if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > randomConeR)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
                float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
                float dEta = track.eta() - randomConeEta;
                if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
                  randomConePt += track.pt();
                }
              }
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutleadingjet"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_rhorandomconewithoutleadingjet"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), rcPhiPsi2, 1.0);
      } else if (nmode == cfgNmodB) {
        continue;
      }
    }
    delete hPtsumSumptFit;
    delete fFitModulationV2v3;
    evtnum += 1;
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPtData, "Sigma pT and bkg as fcn of phi", false);

  void processSigmaPtMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         aod::JetTracks const&)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }

    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetSpectraHistograms(jet, centrality);
    }
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPtMCD, "jet spectra for MCD", false);

  void processSigmaPtAreaSubMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                                soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                                aod::JetTracks const& tracks)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }

    float centrality = -1.0;
    switch (centralityMode) {
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      default:
        centrality = collision.centFT0C();
        break;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
    }

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    fillLeadingJetQA(jets, leadingJetPt, leadingJetPhi, leadingJetEta);

    int nTrk = 0;
    getNtrk(tracks, jets, nTrk, evtnum, leadingJetEta);
    if (nTrk <= 0) {
      return;
    }
    hPtsumSumptFit = new TH1F("h_ptsum_sumpt_fit", "h_ptsum_sumpt fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI);

    fillNtrkCheck(tracks, jets, hPtsumSumptFit, leadingJetEta);

    double ep2 = 0.;
    double ep3 = 0.;
    int cfgNmodA = 2;
    int cfgNmodB = 3;
    int evtPlnAngleA = 7;
    int evtPlnAngleB = 3;
    int evtPlnAngleC = 5;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == cfgNmodA) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      } else if (nmode == cfgNmodB) {
        if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      }
    }

    const char* fitFunctionV2v3 = "[0] * (1. + 2. * ([1] * std::cos(2. * (x - [2])) + [3] * std::cos(3. * (x - [4]))))";
    fFitModulationV2v3 = new TF1("fit_kV3", fitFunctionV2v3, 0, o2::constants::math::TwoPI);
    //=========================< set parameter >=========================//
    fFitModulationV2v3->SetParameter(0, 1.);
    fFitModulationV2v3->SetParameter(1, 0.01);
    fFitModulationV2v3->SetParameter(3, 0.01);

    double ep2fix = 0.;
    double ep3fix = 0.;

    if (ep2 < 0) {
      ep2fix = RecoDecay::constrainAngle(ep2);
      fFitModulationV2v3->FixParameter(2, ep2fix);
    } else {
      fFitModulationV2v3->FixParameter(2, ep2);
    }
    if (ep3 < 0) {
      ep3fix = RecoDecay::constrainAngle(ep3);
      fFitModulationV2v3->FixParameter(4, ep3fix);
    } else {
      fFitModulationV2v3->FixParameter(4, ep3);
    }

    hPtsumSumptFit->Fit(fFitModulationV2v3, "Q", "ep", 0, o2::constants::math::TwoPI);

    double temppara[5];
    temppara[0] = fFitModulationV2v3->GetParameter(0);
    temppara[1] = fFitModulationV2v3->GetParameter(1);
    temppara[2] = fFitModulationV2v3->GetParameter(2);
    temppara[3] = fFitModulationV2v3->GetParameter(3);
    temppara[4] = fFitModulationV2v3->GetParameter(4);

    registry.fill(HIST("h_fitparaRho_evtnum"), evtnum, temppara[0]);
    registry.fill(HIST("h_fitparaPsi2_evtnum"), evtnum, temppara[2]);
    registry.fill(HIST("h_fitparaPsi3_evtnum"), evtnum, temppara[4]);
    registry.fill(HIST("h_v2obs_centrality"), centrality, temppara[1]);
    registry.fill(HIST("h_v3obs_centrality"), centrality, temppara[3]);
    registry.fill(HIST("h_evtnum_centrlity"), evtnum, centrality);

    if (temppara[0] == 0) {
      return;
    }

    int nDF = 1;
    int numOfFreePara = 2;
    nDF = static_cast<int>(fFitModulationV2v3->GetXaxis()->GetNbins()) - numOfFreePara;
    if (nDF == 0 || static_cast<float>(nDF) <= 0.)
      return;
    double chi2 = 0.;
    for (int i = 0; i < hPtsumSumptFit->GetXaxis()->GetNbins(); i++) {
      if (hPtsumSumptFit->GetBinContent(i + 1) <= 0.)
        continue;
      chi2 += std::pow((hPtsumSumptFit->GetBinContent(i + 1) - fFitModulationV2v3->Eval(hPtsumSumptFit->GetXaxis()->GetBinCenter(1 + i))), 2) / hPtsumSumptFit->GetBinContent(i + 1);
    }

    double chiSqr = 999.;
    double cDF = 1.;

    chiSqr = chi2;
    cDF = 1. - chiSquareCDF(nDF, chiSqr);

    int evtCentAreaMin = 0;
    int evtCentAreaMax = 5;
    int evtMidAreaMin = 30;
    int evtMidAreaMax = 50;

    if (cfgChkFitQuality) {
      registry.fill(HIST("h_PvalueCDF_CombinFit"), cDF);
      registry.fill(HIST("h2_PvalueCDFCent_CombinFit"), centrality, cDF);
      registry.fill(HIST("h2_Chi2Cent_CombinFit"), centrality, chiSqr / (static_cast<float>(nDF)));
      registry.fill(HIST("h2_PChi2_CombinFit"), cDF, chiSqr / (static_cast<float>(nDF)));
      if (centrality >= evtCentAreaMin && centrality <= evtCentAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitA"), cDF, chiSqr / (static_cast<float>(nDF)));

      } else if (centrality >= evtMidAreaMin && centrality <= evtMidAreaMax) {
        registry.fill(HIST("h2_PChi2_CombinFitB"), cDF, chiSqr / (static_cast<float>(nDF)));
      }
    }

    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<aod::JetTracks>(jet)) {
          continue;
        }

        if (jet.r() == round(selectedJetsRadius * 100.0f) && jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          double corrBasejetpt = jet.pt() - collision.rho() * jet.area();
          registry.fill(HIST("h_mcd_pt_before_matching_mcdprocess"), corrBasejetpt);
        }

        double integralValue = fFitModulationV2v3->Integral(jet.phi() - jetRadius, jet.phi() + jetRadius);
        double rholocal = collision.rho() / (2 * jetRadius * temppara[0]) * integralValue;
        registry.fill(HIST("h2_rholocal_cent"), centrality, rholocal, 1.0);

        if (nmode == cfgNmodA) {
          double phiMinusPsi2;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          phiMinusPsi2 = jet.phi() - ep2;
          registry.fill(HIST("h2_phi_rholocal"), jet.phi() - ep2, rholocal, 1.0);
          registry.fill(HIST("h_jet_pt_inclusive_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);

          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi2 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2_rho"), centrality, jet.pt() - (rholocal * jet.area()), 1.0);
          }
        } else if (nmode == cfgNmodB) {
          double phiMinusPsi3;
          if (collision.qvecAmp()[detId] < collQvecAmpDetId) {
            continue;
          }
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleA * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= evtPlnAngleB * o2::constants::math::PIQuarter && phiMinusPsi3 < evtPlnAngleC * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
          }
        }
      }
    }
    // RCpT
    for (uint i = 0; i < cfgnMods->size(); i++) {
      TRandom3 randomNumber(0);
      float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
      float randomConePt = 0;
      double integralValueRC = fFitModulationV2v3->Integral(randomConePhi - randomConeR, randomConePhi + randomConeR);
      double rholocalRC = collision.rho() / (2 * randomConeR * temppara[0]) * integralValueRC;

      int nmode = cfgnMods->at(i);
      if (nmode == cfgNmodA) {
        double rcPhiPsi2;
        rcPhiPsi2 = randomConePhi - ep2;

        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphi"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);

        // removing the leading jet from the random cone
        if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
          float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
          float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

          bool jetWasInCone = false;
          while ((randomConeLeadJetDeltaR <= 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
            jetWasInCone = true;
            randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
            randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
            dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
          }
          if (jetWasInCone) {
            randomConePt = 0.0;
            for (auto const& track : tracks) {
              if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > randomConeR)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
                float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
                float dEta = track.eta() - randomConeEta;
                if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
                  randomConePt += track.pt();
                }
              }
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutleadingjet"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_rhorandomconewithoutleadingjet"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), rcPhiPsi2, 1.0);
      } else if (nmode == cfgNmodB) {
        continue;
      }
    }
    delete hPtsumSumptFit;
    delete fFitModulationV2v3;
    evtnum += 1;
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPtAreaSubMCD, "jet spectra with rho-area subtraction for MCD", false);

  void processSigmaPtMCP(soa::Filtered<aod::JetMcCollisions>::iterator const& mccollision,
                         soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                         aod::JetParticles const&)
  {
    bool mcLevelIsParticleLevel = true;

    if (!applyMCCollisionCuts(mccollision, collisions)) {
      return;
    }
    registry.fill(HIST("h_mccollisions_zvertex"), mccollision.posZ());

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
        continue;
      }
      fillMCPHistograms(jet);
    }
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPtMCP, "jet charged v2 for MC particle level", false);

  void processSigmaPtAreaSubMCP(JetBkgRhoMcCollisions::iterator const& mccollision,
                                soa::SmallGroups<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::Qvectors>> const& collisions,
                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                aod::JetTracks const& tracks,
                                aod::JetParticles const&)
  {
    bool mcLevelIsParticleLevel = true;

    if (!applyMCCollisionCuts(mccollision, collisions)) {
      return;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
        continue;
      }
      fillMCPAreaSubHistograms(jet, mccollision.rho());
    }

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    int nTrk = 0;
    for (auto const& collision : collisions) {
      auto collTracks = tracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      fillLeadingJetQAMCP(jets, leadingJetPt, leadingJetPhi, leadingJetEta);

      if (jets.size() > 0) {
        for (auto const& track : collTracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
            nTrk += 1;
          }
        }
      }
      if (nTrk <= 0) {
        return;
      }
      hPtsumSumptFitMCP = new TH1F("h_ptsum_sumpt_fit", "h_ptsum_sumpt fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI);

      if (jets.size() > 0) {
        for (auto const& track : collTracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
            hPtsumSumptFitMCP->Fill(track.phi(), track.pt());
          }
        }
      }
      fitFncAreaSubMCP(collision, jets, hPtsumSumptFitMCP, mcLevelIsParticleLevel);
    }

    delete hPtsumSumptFitMCP;
    delete fFitModulationV2v3P;
    evtnum += 1;
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPtAreaSubMCP, "jet spectra with area-based subtraction for MC particle level", false);

  void processJetsMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          ChargedMCDMatchedJets const& mcdjets,
                          ChargedMCPMatchedJets const&,
                          aod::JetTracks const&, aod::JetParticles const&)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }

    for (const auto& mcdjet : mcdjets) {
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetChargedV2, processJetsMatched, "matched mcp and mcd jets", false);

  void processJetsMatchedWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                  ChargedMCDMatchedJetsWeighted const& mcdjets,
                                  ChargedMCPMatchedJetsWeighted const&,
                                  aod::JetTracks const&, aod::JetParticles const&)
  {
    bool fillHistograms = false;
    bool isWeighted = true;
    float eventWeight = collision.weight();
    if (!applyCollisionCuts(collision, fillHistograms, isWeighted, eventWeight)) {
      return;
    }

    for (const auto& mcdjet : mcdjets) {
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetChargedV2, processJetsMatchedWeighted, "matched mcp and mcd jets with weighted events", false);

  void processJetsMatchedSubtracted(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                                    JetBkgRhoMcCollisions const&,
                                    ChargedMCDMatchedJets const& mcdjets,
                                    ChargedMCPMatchedJets const&,
                                    aod::JetTracks const& tracks, aod::JetParticles const&)
  {
    if (!applyCollisionCuts(collision)) {
      return;
    }

    double mcrho = collision.has_mcCollision() ? collision.mcCollision_as<JetBkgRhoMcCollisions>().rho() : -1;
    registry.fill(HIST("h_mc_rho_matched"), mcrho);

    for (const auto& mcdjet : mcdjets) {
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }

      double leadingJetPt = -1;
      double leadingJetEta = -1;

      if (mcdjet.pt() > leadingJetPt) {
        leadingJetPt = mcdjet.pt();
        leadingJetEta = mcdjet.eta();
      }
      int nTrk = 0;
      if (mcdjet.size() > 0) {
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
            registry.fill(HIST("h_accept_Track"), 4.5);
            nTrk += 1;
          }
        }
      }
      if (nTrk <= 0) {
        return;
      }
      hPtsumSumptFitRM = new TH1F("h_ptsum_sumpt_fit_RM", "h_ptsum_sumpt_RM fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI);
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= localRhoFitPtMin && track.pt() <= localRhoFitPtMax) {
          registry.fill(HIST("h_accept_Track_Match"), 0.5);
          hPtsumSumptFitRM->Fill(track.phi(), track.pt());
        }
      }

      double ep2 = 0.;
      double ep3 = 0.;
      int cfgNmodA = 2;
      int cfgNmodB = 3;

      for (uint i = 0; i < cfgnMods->size(); i++) {
        int nmode = cfgnMods->at(i);
        int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
        if (nmode == cfgNmodA) {
          if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
            ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
          }
        } else if (nmode == cfgNmodB) {
          if (collision.qvecAmp()[detId] > collQvecAmpDetId) {
            ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
          }
        }
      }

      const char* fitFunctionRM = "[0] * (1. + 2. * ([1] * std::cos(2. * (x - [2])) + [3] * std::cos(3. * (x - [4]))))";
      fFitModulationRM = new TF1("fit_kV3", fitFunctionRM, 0, o2::constants::math::TwoPI);
      //=========================< set parameter >=========================//
      fFitModulationRM->SetParameter(0, 1.);
      fFitModulationRM->SetParameter(1, 0.01);
      fFitModulationRM->SetParameter(3, 0.01);

      double ep2fix = 0.;
      double ep3fix = 0.;

      if (ep2 < 0) {
        ep2fix = RecoDecay::constrainAngle(ep2);
        fFitModulationRM->FixParameter(2, ep2fix);
      } else {
        fFitModulationRM->FixParameter(2, ep2);
      }
      if (ep3 < 0) {
        ep3fix = RecoDecay::constrainAngle(ep3);
        fFitModulationRM->FixParameter(4, ep3fix);
      } else {
        fFitModulationRM->FixParameter(4, ep3);
      }

      hPtsumSumptFitRM->Fit(fFitModulationRM, "Q", "ep", 0, o2::constants::math::TwoPI);

      float tempparaA;
      tempparaA = fFitModulationRM->GetParameter(0);

      if (tempparaA == 0) {
        return;
      }
      if (useMedianRho) {
        fillGeoMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, ep2, collision.rho(), mcrho);
      }
      if (useLocalRho) {
        fillGeoMatchedCorrHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, fFitModulationRM, tempparaA, ep2, collision.rho(), subtractMCPBackground, mcrho);
      }
      delete hPtsumSumptFitRM;
      delete fFitModulationRM;
      evtnum += 1;
    }
  }
  PROCESS_SWITCH(JetChargedV2, processJetsMatchedSubtracted, "matched mcp and mcd jets", false);

  void processTracksQA(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                       soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      fillTrackHistograms(track);
    }
  }
  PROCESS_SWITCH(JetChargedV2, processTracksQA, "QA for charged tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetChargedV2>(cfgc)};
}
