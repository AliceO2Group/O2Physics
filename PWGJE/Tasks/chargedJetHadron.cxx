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

/// \file chargedJetHadron.cxx
/// \brief Charged-particle jet - hadron correlation task
/// \author Yongzhen Hou <yongzhen.hou@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/AnalysisDataModel.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THn.h>

#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::math;
const double assoHadronPtCut = 2.0;
static const float kNaN = std::numeric_limits<float>::quiet_NaN();

struct ChargedJetHadron {

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", 0.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 100.0, "maximum centrality"};
  Configurable<float> leadingjetptMin{"leadingjetptMin", 20.0, "minimum leadingjetpt"};
  Configurable<float> subleadingjetptMin{"subleadingjetptMin", 10.0, "minimum subleadingjetpt"};
  Configurable<float> dijetDphiCut{"dijetDphiCut", 0.5, "minimum dijetDphiCut"};
  Configurable<float> assoHadronPtMaxCut{"assoHadronPtMaxCut", 10.0, "maximum associate hadron pt cut"};
  Configurable<float> etaGapdw{"etaGapdw", 0.5, "dijet eta gap low threshold"};
  Configurable<float> etaGapup{"etaGapup", 1.0, "dijet eta gap high threshold"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.7, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<bool> checkLeadConstituentPtForMcpJets{"checkLeadConstituentPtForMcpJets", false, "flag to choose whether particle level jets should have their lead track pt above leadingConstituentPtMin to be accepted; off by default, as leadingConstituentPtMin cut is only applied on MCD jets for the Pb-Pb analysis using pp MC anchored to Pb-Pb for the response matrix"};
  Configurable<bool> doDijetRaa{"doDijetRaa", false, "0: all axis fill of thnsparse, 1: partial filling of thnsparse"};
  Configurable<bool> doEventWeighted{"doEventWeighted", false, "0: weight is 1 for MB Sample, 1: weight from Jet-Jet Sample"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 0, "0:FT0C; 1:FT0A; 2:FT0M"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "number of events mixed in ME process"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0, 100., 300., 600., 1000., 2000., 5000., 8000.}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  ConfigurableAxis binsCentrality{"binsCentrality", {VARIABLE_WIDTH, 0.0, 10., 30., 50, 70., 100.}, "Mixing bins - centrality"};

  // Filter ..................
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut &&
                            aod::jcollision::trackOccupancyInTimeRange >= trackOccupancyInTimeRangeMin && aod::jcollision::trackOccupancyInTimeRange <= trackOccupancyInTimeRangeMax &&
                            ((skipMBGapEvents.node() == false) || (aod::jcollision::getSubGeneratorId != static_cast<int>(jetderiveddatautilities::JCollisionSubGeneratorId::mbGap))));
  Filter mcCollisionFilter = (nabs(aod::jmccollision::posZ) < vertexZCut &&
                              ((skipMBGapEvents.node() == false) || (aod::jmccollision::getSubGeneratorId != static_cast<int>(jetderiveddatautilities::JCollisionSubGeneratorId::mbGap))));
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta >= trackEtaMin && aod::jtrack::eta <= trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta >= trackEtaMin && aod::jmcparticle::eta <= trackEtaMax);
  Filter eventCuts;

  SliceCache cache;
  using FilterCollisions = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>;
  using FilterCollision = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
  using FilterMcpCollision = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator;
  using FilterMcpCollisions = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>;
  using FilterJetTracks = soa::Filtered<aod::JetTracks>;
  using CorrChargedJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
  using CorrChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
  using CorrChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

  /*
  using BinningTypePP = ColumnBinningPolicy<aod::jcollision::PosZ, aod::jcollision::MultFT0M>;
  using BinningTypeMCPP = ColumnBinningPolicy<aod::jmccollision::PosZ, aod::jmccollision::MultFT0M>;
  using BinningType = ColumnBinningPolicy<aod::jcollision::PosZ, aod::jcollision::MultFT0C>;
  using BinningTypeMC = ColumnBinningPolicy<aod::jmccollision::PosZ, aod::jmccollision::MultFT0C>;
  BinningTypePP corrBinning{{binsZVtx, binsMultiplicity}, true};
  BinningTypeMCPP corrBinningMC{{binsZVtx, binsMultiplicityMc}, true};
  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};
  BinningTypeMC corrBinningMC{{binsZVtx, binsMultiplicityMc}, true};
  */

  using BinningTypePP = ColumnBinningPolicy<aod::jcollision::PosZ, aod::jcollision::CentFT0M>;
  using BinningType = ColumnBinningPolicy<aod::jcollision::PosZ, aod::jcollision::CentFT0C>;
  using BinningTypeMC = ColumnBinningPolicy<aod::jmccollision::PosZ, aod::jmccollision::CentFT0M>;
  BinningType corrBinningPP{{binsZVtx, binsCentrality}, true};
  BinningType corrBinning{{binsZVtx, binsCentrality}, true};
  BinningTypeMC corrBinningMC{{binsZVtx, binsCentrality}, true};

  HistogramRegistry registry; // histogram registry
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (cfgCentEstimator == 0) {
      eventCuts = (aod::jcollision::centFT0C >= centralityMin && aod::jcollision::centFT0C < centralityMax);
    } else if (cfgCentEstimator == 1) {
      eventCuts = (aod::jcollision::centFT0A >= centralityMin && aod::jcollision::centFT0A < centralityMax);
    } else {
      eventCuts = (aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);
    }

    AxisSpec centralityAxis = {110, -5., 105., "Centrality"};
    AxisSpec trackPtAxis = {200, 0.0, 200.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {40, -1.0, 1.0, "#eta"};
    AxisSpec phiAxis = {65, -0.2, 6.3, "#varphi"};
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {280, -80., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetmultetaAxis = {4, -0.5, 0.5, "#Delta#eta"};
    AxisSpec detaAxis = {64, -1.6, 1.6, "#Delta#eta"};
    AxisSpec dphiAxis = {70, -1.7, 5.3, "#Delta#varphi"};
    AxisSpec drAxis = {60, 0.0, 1.5, "#Delta#it{R}"};
    AxisSpec axisBdtScore = {100, 0., 1., "Bdt score"};

    if (doprocessCollisionsQCData || doprocessCollisionsQCMCD) {
      registry.add("h_collisions", "event status;event status; entries", {HistType::kTH1F, {{7, 0.0, 7.0}}});
      registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{7, 0.0, 7.0}}});
      registry.add("h_fakecollisions", "event status;event status; entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 30000}}});
      registry.add("h_collisions_Zvertex", "position of collision; #it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_collisions_multFT0", " multiplicity using multFT0; entries", {HistType::kTH1F, {{300, 0, 60000}}});
      registry.add("h2_track_eta_track_phi", "track #eta vs. track #phi; #eta; #phi; counts", {HistType::kTH2F, {etaAxis, phiAxis}});
      registry.add("h2_track_eta_pt", "track #eta vs. track #it{p}_{T}; #eta; #it{p}_{T,track} (GeV/#it{c}; counts", {HistType::kTH2F, {etaAxis, trackPtAxis}});
      registry.add("h2_track_phi_pt", "track #phi vs. track #it{p}_{T}; #phi; #it{p}_{T,track} (GeV/#it{c}; counts", {HistType::kTH2F, {phiAxis, trackPtAxis}});
    }

    if (doprocessSpectraAreaSubData || doprocessSpectraAreaSubMCD) {
      if (doprocessSpectraAreaSubMCD && doEventWeighted) {
        registry.add("h_jet_phat", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
        registry.add("h_jet_phat_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      }
      registry.add("h_jet_pt", "jet pT; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_eta", "jet eta; #eta_{jet}; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_jet_phi", "jet phi; #phi_{jet}; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_jet_area", "jet Area_{jet}; Area_{jet}; counts", {HistType::kTH1F, {{150, 0., 1.5}}});
      registry.add("h_jet_ntracks", "jet N_{jet tracks}; N_{jet, tracks}; counts", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h_jet_pt_rhoareasubtracted", "jet pt; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_eta_rhoareasubtracted", "jet eta; #eta_{jet}; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_jet_phi_rhoareasubtracted", "jet phi; #phi_{jet}; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_jet_area_rhoareasubtracted", "jet Area_{jet}; Area_{jet}; counts", {HistType::kTH1F, {{150, 0., 1.5}}});
      registry.add("h_jet_ntracks_rhoareasubtracted", "jet N_{jet tracks}; N_{jet,tracks}; counts", {HistType::kTH1F, {{200, 0., 200.}}});
    }

    //========jet-hadron correlations======================
    if (doprocessJetHadron || doprocessMixJetHadron || doprocessJetHadronMCD || doprocessMixJetHadronMCD) {
      registry.add("h_trigjet_corrpt", "trigger jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("thn_jeth_correlations", "jet-h correlations; jetpT; trackpT; jeth#Delta#eta; jeth#Delta#varphi; jeth#Delta#it{R}", HistType::kTHnSparseF, {jetPtAxis, trackPtAxis, detaAxis, dphiAxis, drAxis});
      registry.add("h_jeth_event_stats", "Same event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_jeth_event_stats"))->GetXaxis()->SetBinLabel(1, "Total jets");
      registry.get<TH1>(HIST("h_jeth_event_stats"))->GetXaxis()->SetBinLabel(2, "Total jets with pTHat cut");
      registry.get<TH1>(HIST("h_jeth_event_stats"))->GetXaxis()->SetBinLabel(3, "Total jets with cuts");
      registry.get<TH1>(HIST("h_jeth_event_stats"))->GetXaxis()->SetBinLabel(4, "Total j-h pairs");
      registry.get<TH1>(HIST("h_jeth_event_stats"))->GetXaxis()->SetBinLabel(5, "Total j-h pairs with accepted");
      registry.add("h_mixtrigjet_corrpt", "trigger jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("thn_mixjeth_correlations", "ME: jet-h correlations; jetpT; trackpT; jeth#Delta#eta; jeth#Delta#varphi; jeth#Delta#it{R}", HistType::kTHnSparseF, {jetPtAxis, trackPtAxis, detaAxis, dphiAxis, drAxis});
      registry.add("h_mixjeth_event_stats", "Mixed event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_mixjeth_event_stats"))->GetXaxis()->SetBinLabel(1, "Total mixed events");
      registry.get<TH1>(HIST("h_mixjeth_event_stats"))->GetXaxis()->SetBinLabel(2, "Total jets");
      registry.get<TH1>(HIST("h_mixjeth_event_stats"))->GetXaxis()->SetBinLabel(3, "Total jets with cuts");
      registry.get<TH1>(HIST("h_mixjeth_event_stats"))->GetXaxis()->SetBinLabel(4, "Total j-h pairs");
      registry.get<TH1>(HIST("h_mixjeth_event_stats"))->GetXaxis()->SetBinLabel(5, "Total j-h pairs with accepted");
    }

    if (doprocessHFJetCorrelation) {
      registry.add("h_d0jet_pt", "D0 jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_d0jet_corrpt", "D0 jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_d0jet_eta", "D0 jet eta;#eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_d0jet_phi", "D0 jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_d0_pt", ";p_{T,D^{0}};dN/dp_{T,D^{0}}", {HistType::kTH1F, {{200, 0., 10.}}});
      registry.add("h_d0_mass", ";m_{D^{0}} (GeV/c^{2});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("h_d0_eta", ";#eta_{D^{0}} (GeV/c^{2});dN/d#eta_{D^{0}}", {HistType::kTH1F, {{200, -5., 5.}}});
      registry.add("h_d0_phi", ";#varphi_{D^{0}} (GeV/c^{2});dN/d#varphi_{D^{0}}", {HistType::kTH1F, {{200, -10., 10.}}});
      registry.add("h_d0_bdtScorePrompt", "D0 BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
      registry.add("h_d0_bdtScoreBkg", "D0 BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
      registry.add("h_d0bar_bdtScorePrompt", "D0bar BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
      registry.add("h_d0bar_bdtScoreBkg", "D0bar BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
      registry.add("h2_d0jet_detadphi", "D^{0}-jets deta vs dphi; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
    }
    //========leading jet-hadron correlations======================
    if (doprocessLeadingJetHadron || doprocessLeadinJetHadronMCD) {
      registry.add("h_centrality", "centrality distributions; centrality; counts", {HistType::kTH1F, {centralityAxis}});
      registry.add("h_inclusivejet_corrpt", "inclusive jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_leadjet_pt", "leading jet pT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_leadjet_corrpt", "leading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_subleadjet_pt", "subleading jet pT;#it{p}_{T,subleadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_subleadjet_corrpt", "subleading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_leadjet_eta", "leading jet eta;#eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_leadjet_phi", "leading jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_subleadjet_eta", "subleading jet eta;#eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_subleadjet_phi", "subleading jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h2_dijet_detanoflip_dphi", "dijet #Delta#eta no flip vs #Delta#varphi; #Delta#eta_{noflip}; #Delta#varphi; counts", {HistType::kTH2F, {detaAxis, {63, 0, 6.3}}});
      registry.add("h2_dijet_Asymmetry", "dijet Asymmetry; #it{p}_{T,subleadingjet} (GeV/#it{c}); #it{X}_{J}; counts", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {40, 0, 1.0}}});
      registry.add("h3_dijet_deta_pt", "dijet #Delta#eta flip vs #it{p}_{T,jet1-jet2}; #Delta#eta_{flip}; #it{p}_{T,jet1} (GeV/#it{c}); #it{p}_{T,jet2} (GeV/#it{c})", {HistType::kTH3F, {{16, 0, 1.6}, jetPtAxis, jetPtAxis}});
      registry.add("h2_dijet_TimeEtaThan0_pt", "dijet #eta_{jet1}#eta_{jet1} > 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_dijet_TimeEtaLess0_pt", "dijet #eta_{jet1}#eta_{jet1} < 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});

      registry.add("h_jeth_detatot", "jet-hadron tot #Delta#eta;#Delta#eta; counts", {HistType::kTH1F, {detaAxis}});
      registry.add("h_jeth_deta", "jet-hadron #Delta#eta;#Delta#eta; counts", {HistType::kTH1F, {detaAxis}});
      registry.add("h_jeth_dphi", "jet-hadron #Delta#varphi;#Delta#varphi; counts", {HistType::kTH1F, {dphiAxis}});
      registry.add("h2_jeth_deta_dphi", "jeth deta vs dphi; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsup_deta_dphi", "jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| > 1.0; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsmd_deta_dphi", "jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| #in (0.5, 1.0); #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsdw_deta_dphi", "jeth deta vs dphi with physical cuts  |#Delta#eta_{jet1,2}| < 0.5; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("thn_ljeth_correlations", "leading jet-h correlations; leadingjetpT; subleadingjetpT; trackpT; timedijeteta; #Delta#eta_{jet1,2}; track #eta; jeth#Delta#eta; jeth#Delta#varphi", HistType::kTHnSparseF, {jetPtAxis, jetPtAxis, {10, 0., 10.}, jetmultetaAxis, {16, 0, 1.6}, etaAxis, detaAxis, dphiAxis});
    }

    if (doprocessMixLeadingJetHadron || doprocessMixLeadinJetHadronMCD) {
      registry.add("h_mixdijet_pair_counts_cut", "ME: number of pairs with leadingjet & subleadingjet cut pair; jet pairs; counts", {HistType::kTH1F, {{10, 0, 10}}});
      registry.add("h_mixleadjet_corrpt", "ME: leading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mixsubleadjet_corrpt", "ME: subleading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mixleadjet_eta", "ME: leading jet eta; #eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_mixsubleadjet_eta", "ME: subleading jet eta; #eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h2_mixdijet_detanoflip_dphi", "ME: dijet #Delta#eta no flip vs #Delta#varphi; #Delta#eta_{noflip}; #Delta#varphi; counts", {HistType::kTH2F, {detaAxis, {63, 0, 6.3}}});
      registry.add("h2_mixdijet_Asymmetry", "ME: dijet Asymmetry; #it{p}_{T,subleadingjet} (GeV/#it{c}); #it{X}_{J}; counts", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {40, 0, 1.0}}});
      registry.add("h3_mixdijet_deta_pt", "ME: dijet #Delta#eta flip vs #it{p}_{T,jet1-jet2}; #Delta#eta_{flip}; #Delta#varphi; counts", {HistType::kTH3F, {{16, 0, 1.6}, jetPtAxis, jetPtAxis}});
      registry.add("h2_mixdijet_TimeEtaThan0_pt", "dijet #eta_{jet1}#eta_{jet1} > 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_mixdijet_TimeEtaLess0_pt", "dijet #eta_{jet1}#eta_{jet1} < 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h_mixjeth_detatot", "ME: jet-hadron correlations; no flip #Delta#eta", {HistType::kTH1F, {detaAxis}});
      registry.add("h_mixjeth_deta", "ME: jet-hadron correlations; #Delta#eta", {HistType::kTH1F, {detaAxis}});
      registry.add("h_mixjeth_dphi", "ME: jet-hadron correlations; #Delta#phi", {HistType::kTH1F, {dphiAxis}});
      registry.add("h2_mixjeth_deta_dphi", "ME: jet-hadron correlations; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsup_deta_dphi", "ME: jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| > 1.0; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsmd_deta_dphi", "ME: jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| #in (0.5, 1.0); #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsdw_deta_dphi", "ME: jeth deta vs dphi with physical cuts  |#Delta#eta_{jet1,2}| < 0.5; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("thn_mixljeth_correlations", "ME: leading jet-h correlations; leadingjetpT; subleadingjetpT; trackpT; timedijeteta; #Delta#eta_{jet1,2}; track #eta; jeth#Delta#eta; jeth#Delta#varphi; poolBin", HistType::kTHnSparseF, {jetPtAxis, jetPtAxis, {10, 0., 10.}, jetmultetaAxis, {16, 0, 1.6}, etaAxis, detaAxis, dphiAxis, {15, 0, 15}});

      registry.add("h_mix_event_stats", "Mixed event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_mix_event_stats"))->GetXaxis()->SetBinLabel(1, "Total mixed events");
      registry.get<TH1>(HIST("h_mix_event_stats"))->GetXaxis()->SetBinLabel(2, "Total dijets");
      registry.get<TH1>(HIST("h_mix_event_stats"))->GetXaxis()->SetBinLabel(3, "Total dijets with cuts");
      registry.get<TH1>(HIST("h_mix_event_stats"))->GetXaxis()->SetBinLabel(4, "Total Lj-h pairs");
      registry.get<TH1>(HIST("h_mix_event_stats"))->GetXaxis()->SetBinLabel(5, "Total Lj-h pairs with cut");
    }

    if (doprocessCollisionsQCMCP) {
      registry.add("h_mcColl_counts", " number of mc events; event status; entries", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(4, "recoEvtSel");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(5, "centralitycut");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(6, "occupancycut");
      registry.add("h_mcdColl_mult", " mcd multiplicity global tracks; entries", {HistType::kTH1F, {{300, 0, 60000}}});

      registry.add("h_mcpColl_Zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_mcpColl_centrality", "mcp collision centrality; centrality; counts", {HistType::kTH1F, {centralityAxis}});
      registry.add("h_mcpColl_mult", " mcp multiplicity global tracks; entries", {HistType::kTH1F, {{300, 0, 60000}}});
      registry.add("h2_particle_eta_phi", "particle #eta vs. particle #phi; #eta; #phi; counts", {HistType::kTH2F, {etaAxis, phiAxis}});
      registry.add("h2_particle_eta_pt", "particle #eta vs. particle #it{p}_{T}; #eta; #it{p}_{T,particle} (GeV/#it{c}; counts", {HistType::kTH2F, {etaAxis, trackPtAxis}});
      registry.add("h2_particle_phi_pt", "particle #phi vs. particle #it{p}_{T}; #phi; #it{p}_{T,particle} (GeV/#it{c}; counts", {HistType::kTH2F, {phiAxis, trackPtAxis}});
    }
    if (doprocessSpectraAreaSubMCP) {
      if (doEventWeighted) {
        registry.add("h_mcColl_counts_weight", " number of weighted mc events; event status; entries", {HistType::kTH1F, {{7, 0., 7.}}});
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(1, "McColl");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(2, "xsectGen");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(3, "event weight");
      }
      registry.add("h_mcColl_rho", "mc collision rho;#rho (GeV/#it{c}); counts", {HistType::kTH1F, {{500, 0.0, 500.0}}});
      registry.add("h_jet_pt_part", "partvjet pT;#it{p}_{T,jet}^{part} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_eta_part", "part jet #eta;#eta^{part}; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_jet_phi_part", "part jet #varphi;#phi^{part}; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_jet_area_part", "part jet Area_{jet}; Area_{jet}^{part}; counts", {HistType::kTH1F, {{150, 0., 1.5}}});
      registry.add("h_jet_ntracks_part", "part jet N_{jet tracks}; N_{jet, tracks}^{part}; counts", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h_jet_pt_part_rhoareasubtracted", "part jet corr pT;#it{p}_{T,jet}^{part} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_eta_part_rhoareasubtracted", "part jet #eta;#eta^{part}; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_jet_phi_part_rhoareasubtracted", "part jet #varphi;#varphi^{part}; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_jet_area_part_rhoareasubtracted", "part jet Area_{jet}; Area_{jet}^{part}; counts", {HistType::kTH1F, {{150, 0., 1.5}}});
      registry.add("h_jet_ntracks_part_rhoareasubtracted", "part jet N_{jet tracks}; N_{jet, tracks}^{part}; counts", {HistType::kTH1F, {{200, -0.5, 199.5}}});
    }

    if (doprocessJetHadronMCP || doprocessMixJetHadronMCP) {
      //.........MCP: jet-hadron correlations.................
      registry.add("h_trigjet_corrpt_part", "trigger jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("thn_jeth_correlations_part", "MCP: jet-h correlations; jetpT; trackpT; jeth#Delta#eta; jeth#Delta#varphi; jeth#Delta#it{R}", HistType::kTHnSparseF, {jetPtAxis, trackPtAxis, detaAxis, dphiAxis, drAxis});

      registry.add("h_jeth_event_stats_part", "MCP: same event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_jeth_event_stats_part"))->GetXaxis()->SetBinLabel(1, "Total jets");
      registry.get<TH1>(HIST("h_jeth_event_stats_part"))->GetXaxis()->SetBinLabel(2, "Total jets with pTHat cut");
      registry.get<TH1>(HIST("h_jeth_event_stats_part"))->GetXaxis()->SetBinLabel(3, "Total jets with cuts");
      registry.get<TH1>(HIST("h_jeth_event_stats_part"))->GetXaxis()->SetBinLabel(4, "Total j-h pairs");
      registry.get<TH1>(HIST("h_jeth_event_stats_part"))->GetXaxis()->SetBinLabel(5, "Total j-h pairs with accepted");

      registry.add("h_mixtrigjet_corrpt_part", "trigger jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("thn_mixjeth_correlations_part", "mcpME: jet-h correlations; jetpT; trackpT; jeth#Delta#eta; jeth#Delta#varphi; jeth#Delta#it{R}", HistType::kTHnSparseF, {jetPtAxis, trackPtAxis, detaAxis, dphiAxis, drAxis});
      registry.add("h_mixjeth_event_stats_part", "MCP: mixed event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_mixjeth_event_stats_part"))->GetXaxis()->SetBinLabel(1, "Total mixed events");
      registry.get<TH1>(HIST("h_mixjeth_event_stats_part"))->GetXaxis()->SetBinLabel(2, "Total jets");
      registry.get<TH1>(HIST("h_mixjeth_event_stats_part"))->GetXaxis()->SetBinLabel(3, "Total jets with cuts");
      registry.get<TH1>(HIST("h_mixjeth_event_stats_part"))->GetXaxis()->SetBinLabel(4, "Total j-h pairs");
      registry.get<TH1>(HIST("h_mixjeth_event_stats_part"))->GetXaxis()->SetBinLabel(5, "Total j-h pairs with accepted");
    }

    if (doprocessLeadingJetHadronMCP) {
      //.........SE leading jet correlations...............
      registry.add("h_leadjet_pt_part", "MCP: leading jet pT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_leadjet_corrpt_part", "MCP: leading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_leadjet_eta_part", "MCP: leading jet eta;#eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_leadjet_phi_part", "MCP: leading jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h_subleadjet_pt_part", "MCP: subleading jet pT;#it{p}_{T,subleadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_subleadjet_corrpt_part", "MCP: subleading jet corrpT; #it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_subleadjet_eta_part", "MCP: subleading jet eta;#eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_subleadjet_phi_part", "MCP: subleading jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h2_dijet_detanoflip_dphi_part", "MCP: dijet #Delta#eta no flip vs #Delta#varphi; #Delta#eta_{noflip}; #Delta#varphi; counts", {HistType::kTH2F, {detaAxis, {63, 0, 6.3}}});
      registry.add("h2_dijet_Asymmetry_part", "MCP: dijet Asymmetry; #it{p}_{T,subleadingjet} (GeV/#it{c}); #it{X}_{J}; counts", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {40, 0, 1.0}}});
      registry.add("h3_dijet_deta_pt_part", "MCP: dijet #Delta#eta flip vs #it{p}_{T,jet1-jet2}; #Delta#eta_{flip}; #Delta#varphi; counts", {HistType::kTH3F, {{16, 0, 1.6}, jetPtAxis, jetPtAxis}});
      registry.add("h2_dijet_TimeEtaThan0_pt_part", "dijet #eta_{jet1}#eta_{jet1} > 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_dijet_TimeEtaLess0_pt_part", "dijet #eta_{jet1}#eta_{jet1} < 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});

      registry.add("h_jeth_detatot_part", "MCP: jet-hadron tot #Delta#eta;#Delta#eta; counts", {HistType::kTH1F, {detaAxis}});
      registry.add("h_jeth_deta_part", "MCP: jet-hadron #Delta#eta;#Delta#eta; counts", {HistType::kTH1F, {detaAxis}});
      registry.add("h_jeth_dphi_part", "MCP: jet-hadron #Delta#varphi;#Delta#varphi; counts", {HistType::kTH1F, {dphiAxis}});
      registry.add("h2_jeth_deta_dphi_part", "MCP: jeth deta vs dphi; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsup_deta_dphi_part", "MCP: jeth deta vs dphi with physical cuts |#Delta#eta_{jet}| > 1.0; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsmd_deta_dphi_part", "MCP: jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| #in (0.5, 1.0); #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_jeth_physicalcutsdw_deta_dphi_part", "MCP: jeth deta vs dphi with physical cuts  |#Delta#eta_{jet1,2}| < 0.5; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("thn_ljeth_correlations_part", "MCP: leading jet-h correlations; leadingjetpT; subleadingjetpT; trackpT; timedijeteta; #Delta#eta_{jet1,2}; track #eta; jeth#Delta#eta; jeth#Delta#varphi", HistType::kTHnSparseF, {jetPtAxis, jetPtAxis, {10, 0., 10.}, jetmultetaAxis, {16, 0, 1.6}, etaAxis, detaAxis, dphiAxis});
    }

    if (doprocessMixLeadingJetHadronMCP) {
      //...........mcp mixed events: leading jet correlations.................
      registry.add("h_mixdijet_pair_counts_cut_part", "mcpME: number of pairs with leadingjet & subleadingjet cut pair; jet pairs; counts", {HistType::kTH1F, {{10, 0, 10}}});
      registry.add("h_mixleadjet_corrpt_part", "mcpME: leading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mixsubleadjet_corrpt_part", "mcpME: subleading jet corrpT;#it{p}_{T,leadingjet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_mixleadjet_eta_part", "mcpME: leading jet eta; #eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h_mixsubleadjet_eta_part", "mcpME: subleading jet eta; #eta; counts", {HistType::kTH1F, {etaAxis}});
      registry.add("h2_mixdijet_detanoflip_dphi_part", "mcpME: dijet #Delta#eta no flip vs #Delta#varphi; #Delta#eta_{noflip}; #Delta#varphi; counts", {HistType::kTH2F, {detaAxis, {63, 0, 6.3}}});
      registry.add("h2_mixdijet_Asymmetry_part", "mcpME: dijet Asymmetry; #it{p}_{T,subleadingjet} (GeV/#it{c}); #it{X}_{J}; counts", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {40, 0, 1.0}}});
      registry.add("h3_mixdijet_deta_pt_part", "mcpME: dijet #Delta#eta flip vs #it{p}_{T,jet1-jet2}; #Delta#eta_{flip}; #Delta#varphi; counts", {HistType::kTH3F, {{16, 0, 1.6}, jetPtAxis, jetPtAxis}});
      registry.add("h2_mixdijet_TimeEtaThan0_pt_part", "dijet #eta_{jet1}#eta_{jet1} > 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_mixdijet_TimeEtaLess0_pt_part", "dijet #eta_{jet1}#eta_{jet1} < 0", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});

      registry.add("h_mixjeth_detatot_part", "mcpME: jet-hadron correlations; no flip #Delta#eta", {HistType::kTH1F, {detaAxis}});
      registry.add("h_mixjeth_deta_part", "mcpME: jet-hadron correlations; #Delta#eta", {HistType::kTH1F, {detaAxis}});
      registry.add("h_mixjeth_dphi_part", "mcpME: jet-hadron correlations; #Delta#phi", {HistType::kTH1F, {dphiAxis}});
      registry.add("h2_mixjeth_deta_dphi_part", "mcpME: jet-hadron correlations; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsup_deta_dphi_part", "mcpME: jeth deta vs dphi with physical cuts |#Delta#eta_{jet}| > 1.0; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsmd_deta_dphi_part", "mcpME: jeth deta vs dphi with physical cuts |#Delta#eta_{jet1,2}| #in (0.5, 1.0); #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("h2_mixjeth_physicalcutsdw_deta_dphi_part", "mcpME: jeth deta vs dphi with physical cuts  |#Delta#eta_{jet1,2}| < 0.5; #Delta#eta; #Delta#phi", {HistType::kTH2F, {detaAxis, dphiAxis}});
      registry.add("thn_mixljeth_correlations_part", "mcpME: leading jet-h correlations; leadingjetpT; subleadingjetpT; trackpT; timedijeteta; #Delta#eta_{jet1,2}; track #eta; jeth#Delta#eta; jeth#Delta#varphi", HistType::kTHnSparseF, {jetPtAxis, jetPtAxis, {10, 0., 10.}, jetmultetaAxis, {16, 0, 1.6}, etaAxis, detaAxis, dphiAxis});

      registry.add("h_mixevent_stats_part", "MCP: mixed event statistics; Event pair type; counts", {HistType::kTH1F, {{7, 0., 7.}}});
      registry.get<TH1>(HIST("h_mixevent_stats_part"))->GetXaxis()->SetBinLabel(1, "Total mixed events");
      registry.get<TH1>(HIST("h_mixevent_stats_part"))->GetXaxis()->SetBinLabel(2, "Total dijets");
      registry.get<TH1>(HIST("h_mixevent_stats_part"))->GetXaxis()->SetBinLabel(3, "Total dijets with cuts");
      registry.get<TH1>(HIST("h_mixevent_stats_part"))->GetXaxis()->SetBinLabel(4, "Total Lj-h pairs");
      registry.get<TH1>(HIST("h_mixevent_stats_part"))->GetXaxis()->SetBinLabel(5, "Total Lj-h pairs with cut");
    }
    if (!(acceptSplitCollisions == NonSplitOnly || acceptSplitCollisions == SplitOkCheckAnyAssocColl || acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly)) {
      LOGF(fatal, "Configurable acceptSplitCollisions has wrong input value; stopping workflow");
    }
  }
  // ==========================================================
  // getcentrality getmultiplicity
  // ==========================================================
  template <typename TCollision>
  float getCentrality(const TCollision& coll) const
  {
    if (cfgCentEstimator.value == 0)
      return coll.centFT0C();
    if (cfgCentEstimator.value == 1)
      return coll.centFT0A();
    return coll.centFT0M();
  }

  template <typename TCollision>
  float getMultiplicity(const TCollision& coll) const
  {
    if (cfgCentEstimator.value == 0)
      return coll.multFT0C();
    if (cfgCentEstimator.value == 1)
      return coll.multFT0A();
    return coll.multFT0C() + coll.multFT0A();
  }
  // ==========================================================
  // event selection, vertexZ, occupancy, centrality
  // ==========================================================
  template <typename TCollision>
  bool isGoodCollision(const TCollision& coll) const
  {
    if (!jetderiveddatautilities::selectCollision(coll, eventSelectionBits, skipMBGapEvents.value))
      return false;
    const auto occ = coll.trackOccupancyInTimeRange();
    if (occ < trackOccupancyInTimeRangeMin.value || occ > trackOccupancyInTimeRangeMax.value)
      return false;
    float cent = getCentrality(coll);
    if (cent < centralityMin.value || cent > centralityMax.value)
      return false;
    if (std::abs(coll.posZ()) > vertexZCut.value)
      return false;

    return true;
  }

  template <typename TMcCollision, typename TCollisions>
  bool applyMCCollisionCuts(const TMcCollision& mccollision, const TCollisions& collisions) const
  {
    // MC z-vertex cut, must have associated collisions, split-collision rule
    if (std::abs(mccollision.posZ()) > vertexZCut.value)
      return false;

    if (collisions.size() < 1)
      return false;
    if (acceptSplitCollisions.value == NonSplitOnly && collisions.size() > 1)
      return false;

    // At least one associated collision must pass all cuts
    if (acceptSplitCollisions.value == SplitOkCheckFirstAssocCollOnly)
      return isGoodCollision(*collisions.begin());

    for (auto const& collision : collisions) {
      if (isGoodCollision(collision)) {
        return true;
      }
    }
    return false;
  }

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(const TJets& jet, bool mcLevelIsParticleLevel = false)
  {
    double jetAreaLimit = -98.0, ptMinDefault = -98.0, ptMaxDefault = 9998.0;
    // Jet area cut
    if (jetAreaFractionMin > jetAreaLimit && jet.area() < jetAreaFractionMin * PI * (jet.r() * 0.01) * (jet.r() * 0.01))
      return false;
    // Leading constituent pT selection
    const float ptMin = leadingConstituentPtMin;
    const float ptMax = leadingConstituentPtMax;
    if ((ptMin <= ptMinDefault && ptMax >= ptMaxDefault) || (mcLevelIsParticleLevel && !checkLeadConstituentPtForMcpJets))
      return true;
    // Loop jet constituents
    for (const auto& constituent : jet.template tracks_as<TTracks>()) {
      const double pt = constituent.pt();
      if (pt >= ptMin && pt <= ptMax)
        return true;
    }
    return false;
  }
  // ==========================================================
  template <typename TTracks>
  void fillTrackHistograms(TTracks const& track, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (track.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    if (doEventWeighted) {
      registry.fill(HIST("h_jet_phat"), pTHat);
      registry.fill(HIST("h_jet_phat_weighted"), pTHat, weight);
    }
    registry.fill(HIST("h2_track_eta_track_phi"), track.eta(), track.phi(), weight);
    registry.fill(HIST("h2_track_eta_pt"), track.eta(), track.pt(), weight);
    registry.fill(HIST("h2_track_phi_pt"), track.phi(), track.pt(), weight);
  }

  template <typename TParticles>
  void fillParticleHistograms(const TParticles& particle, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (particle.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    registry.fill(HIST("h2_particle_eta_phi"), particle.eta(), particle.phi(), weight);
    registry.fill(HIST("h2_particle_eta_pt"), particle.eta(), particle.pt(), weight);
    registry.fill(HIST("h2_particle_phi_pt"), particle.phi(), particle.pt(), weight);
  }

  template <typename TJets>
  void fillJetHistograms(TJets const& jet, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_area"), jet.area(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracksIds().size(), weight);
    }
  }

  template <typename TJets>
  void fillJetAreaSubHistograms(TJets const& jet, float rho, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    double jetcorrpt = jet.pt() - (rho * jet.area());
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jetcorrpt, weight);
      if (jetcorrpt > 0) {
        registry.fill(HIST("h_jet_eta_rhoareasubtracted"), jet.eta(), weight);
        registry.fill(HIST("h_jet_phi_rhoareasubtracted"), jet.phi(), weight);
        registry.fill(HIST("h_jet_area_rhoareasubtracted"), jet.area(), weight);
        registry.fill(HIST("h_jet_ntracks_rhoareasubtracted"), jet.tracksIds().size(), weight);
      }
    }
  }

  template <typename TJets>
  void fillMCPHistograms(TJets const& jet, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_area_part"), jet.area(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size(), weight);
    }
  }

  template <typename TJets>
  void fillMCPAreaSubHistograms(TJets const& jet, float rho = 0.0, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
      return;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      double jetcorrpt = jet.pt() - (rho * jet.area());
      registry.fill(HIST("h_jet_pt_part_rhoareasubtracted"), jetcorrpt, weight);
      if (jetcorrpt > 0) {
        registry.fill(HIST("h_jet_eta_part_rhoareasubtracted"), jet.eta(), weight);
        registry.fill(HIST("h_jet_phi_part_rhoareasubtracted"), jet.phi(), weight);
        registry.fill(HIST("h_jet_area_part_rhoareasubtracted"), jet.area(), weight);
        registry.fill(HIST("h_jet_ntracks_part_rhoareasubtracted"), jet.tracksIds().size(), weight);
      }
    }
  }

  // ==========================================================
  //..........jet - hadron correlations........................
  // ==========================================================
  template <typename TCollision, typename TJets, typename TTracks>
  void fillJetHadronHistograms(const TCollision& collision, const TJets& jets, const TTracks& tracks, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    using TracksTable = std::decay_t<decltype(tracks)>;
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<TracksTable>(jet))
        continue;
      registry.fill(HIST("h_jeth_event_stats"), 1);
      if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
        continue;
      registry.fill(HIST("h_jeth_event_stats"), 2);
      double ptCorr = jet.pt() - jet.area() * collision.rho();
      if (ptCorr < subleadingjetptMin)
        continue;
      registry.fill(HIST("h_trigjet_corrpt"), ptCorr);
      registry.fill(HIST("h_jeth_event_stats"), 3);
      for (auto const& track : tracks) {
        registry.fill(HIST("h_jeth_event_stats"), 4);
        if (!jetderiveddatautilities::selectTrack(track, trackSelection))
          continue;
        registry.fill(HIST("h_jeth_event_stats"), 5);
        double deta = track.eta() - jet.eta();
        double dphi = track.phi() - jet.phi();
        dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
        double dr = std::sqrt(deta * deta + dphi * dphi);
        registry.fill(HIST("thn_jeth_correlations"), ptCorr, track.pt(), deta, dphi, dr, weight);
      }
    }
  }

  //.......mixed events.................................
  template <typename TCollisions, typename TJets, typename TTracks>
  void fillMixJetHadronHistograms(const TCollisions& collisions, const TJets& jets, const TTracks& tracks, float weight = 1.0)
  {
    using TracksTable = std::decay_t<decltype(tracks)>;
    auto tracksTuple = std::make_tuple(jets, tracks);
    Pair<TCollisions, TJets, TTracks, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache}; //// -1 is the number of the bin to skip

    for (const auto& [c1, jets1, c2, tracks2] : pairData) {
      weight = doEventWeighted ? c1.weight() : 1.f;
      const float pTHat = 10.f / std::pow(weight, 1.f / pTHatExponent);
      registry.fill(HIST("h_mixjeth_event_stats"), 1);
      if (!isGoodCollision(c1) || !isGoodCollision(c2))
        continue;
      for (auto const& jet : jets1) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<TracksTable>(jet))
          continue; // for pp Reason: Trying to dereference index with a wrong type in tracks_as<T> for base target "JTracks"
        if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
          continue;
        registry.fill(HIST("h_mixjeth_event_stats"), 2);
        double ptCorr = jet.pt() - jet.area() * c1.rho();
        if (ptCorr < subleadingjetptMin)
          continue;
        registry.fill(HIST("h_mixtrigjet_corrpt"), ptCorr);
        registry.fill(HIST("h_mixjeth_event_stats"), 3);
        for (auto const& track : tracks2) {
          registry.fill(HIST("h_mixjeth_event_stats"), 4);
          if (!jetderiveddatautilities::selectTrack(track, trackSelection))
            continue;
          registry.fill(HIST("h_mixjeth_event_stats"), 5);
          double deta = track.eta() - jet.eta();
          double dphi = track.phi() - jet.phi();
          dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
          double dr = std::sqrt(deta * deta + dphi * dphi);
          registry.fill(HIST("thn_mixjeth_correlations"), ptCorr, track.pt(), deta, dphi, dr, weight);
        }
      }
    }
  }

  //........MCP..jet - hadron correlations..........................................
  template <typename TmcCollision, typename TJets, typename TParticles>
  void fillMCPJetHadronHistograms(const TmcCollision& mccollision, const TJets& jets, const TParticles& particles, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    using ParticlesTable = std::decay_t<decltype(particles)>;
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<ParticlesTable>(jet, true)) {
        continue;
      }
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
        return;
      registry.fill(HIST("h_jeth_event_stats_part"), 2);
      double ptCorr = jet.pt() - jet.area() * mccollision.rho();
      if (ptCorr < subleadingjetptMin)
        continue;
      registry.fill(HIST("h_trigjet_corrpt_part"), ptCorr);
      registry.fill(HIST("h_jeth_event_stats_part"), 3);
      for (auto const& particle : particles) {
        registry.fill(HIST("h_jeth_event_stats_part"), 4);
        double deta = particle.eta() - jet.eta();
        double dphi = particle.phi() - jet.phi();
        dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
        double dr = std::sqrt(deta * deta + dphi * dphi);
        registry.fill(HIST("thn_jeth_correlations_part"), ptCorr, particle.pt(), deta, dphi, dr, weight);
      }
    }
  }
  //......mixed events......................
  template <typename TmcCollisions, typename TCollisions, typename TJets, typename TParticles>
  void fillMCPMixJetHadronHistograms(const TmcCollisions& mccollisions, const TCollisions& collisions, const TJets& jets, const TParticles& particles, float weight = 1.0)
  {
    using ParticlesTable = std::decay_t<decltype(particles)>;
    auto particlesTuple = std::make_tuple(jets, particles);
    Pair<TmcCollisions, TJets, TParticles, BinningTypeMC> pairMCData{corrBinningMC, numberEventsMixed, -1, mccollisions, particlesTuple, &cache};

    for (const auto& [c1, jets1, c2, particles2] : pairMCData) {
      weight = doEventWeighted ? c1.weight() : 1.f;
      const float pTHat = 10.f / std::pow(weight, 1.f / pTHatExponent);
      registry.fill(HIST("h_mixjeth_event_stats_part"), 1);
      if (!applyMCCollisionCuts(c1, collisions) || !applyMCCollisionCuts(c2, collisions))
        continue;

      for (auto const& jet : jets1) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
          continue;
        if (!isAcceptedJet<ParticlesTable>(jet, true))
          continue;
        if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
          return;
        registry.fill(HIST("h_mixjeth_event_stats_part"), 2);

        double ptCorr = jet.pt() - jet.area() * c1.rho();
        if (ptCorr < subleadingjetptMin)
          continue;
        registry.fill(HIST("h_mixtrigjet_corrpt_part"), ptCorr);
        registry.fill(HIST("h_mixjeth_event_stats_part"), 3);

        for (auto const& particle : particles2) {
          registry.fill(HIST("h_mixjeth_event_stats_part"), 4);
          double deta = particle.eta() - jet.eta();
          double dphi = particle.phi() - jet.phi();
          dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
          double dr = std::sqrt(deta * deta + dphi * dphi);
          registry.fill(HIST("thn_mixjeth_correlations_part"), ptCorr, particle.pt(), deta, dphi, dr, weight);
        }
      }
    }
  }

  // ==========================================================
  //..........leading jet - hadron correlations................
  // ==========================================================
  template <typename TCollision, typename TJets, typename TTracks>
  void fillLeadingJetHadronHistograms(const TCollision& collision, const TJets& jets, const TTracks& tracks, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    using TracksTable = std::decay_t<decltype(tracks)>;
    registry.fill(HIST("h_centrality"), getCentrality(collision));
    typename TJets::iterator leadingJet;
    typename TJets::iterator subleadingJet;
    bool hasLeading = false;
    bool hasSubleading = false;
    double ptLeadingCorr = -1.0;
    double ptSubleadingCorr = -1.0;

    for (auto it = jets.begin(); it != jets.end(); ++it) {
      const auto& jet = *it;
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (!isAcceptedJet<TracksTable>(jet))
        continue;
      if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
        continue;
      double ptCorr = jet.pt() - jet.area() * collision.rho();
      registry.fill(HIST("h_inclusivejet_corrpt"), ptCorr, weight);
      if (ptCorr > ptLeadingCorr) {
        if (hasLeading) {
          subleadingJet = leadingJet;
          ptSubleadingCorr = ptLeadingCorr;
          hasSubleading = true;
        }
        leadingJet = it;
        ptLeadingCorr = ptCorr;
        hasLeading = true;
      } else if (ptCorr > ptSubleadingCorr) {
        subleadingJet = it;
        ptSubleadingCorr = ptCorr;
        hasSubleading = true;
      }
    }
    if (!hasLeading || !hasSubleading)
      return; // fewer than 2 jets

    double phiLead = leadingJet.phi();
    double phiSub = subleadingJet.phi();
    double deltaPhiJets = phiLead - phiSub;
    deltaPhiJets = RecoDecay::constrainAngle(deltaPhiJets, -PIHalf);
    if (std::abs(deltaPhiJets) < dijetDphiCut * PI)
      return;
    // === Step2: eta ordering etajet1 > etajet2) ===
    double etaJet1Raw = leadingJet.eta();
    double etaJet2Raw = subleadingJet.eta();
    double multEta1Eta2 = etaJet1Raw * etaJet2Raw;
    double deltaEtaJetsNoflip = etaJet1Raw - etaJet2Raw;
    double inverse = (etaJet1Raw > 0) ? 1.0 : -1.0; // Dr.Yang suggestion
    double flip = (etaJet1Raw > etaJet2Raw) ? 1.0 : -1.0;
    double etajet1 = flip * etaJet1Raw;      // leading jet eta after flip
    double etajet2 = flip * etaJet2Raw;      // subleading jet eta after flip
    double deltaEtaJets = etajet1 - etajet2; // >= 0
    registry.fill(HIST("h_leadjet_pt"), leadingJet.pt(), weight);
    registry.fill(HIST("h_subleadjet_pt"), subleadingJet.pt(), weight);
    registry.fill(HIST("h_leadjet_corrpt"), ptLeadingCorr, weight);
    registry.fill(HIST("h_subleadjet_corrpt"), ptSubleadingCorr, weight);

    if (ptLeadingCorr < leadingjetptMin || ptSubleadingCorr < subleadingjetptMin)
      return;
    registry.fill(HIST("h_leadjet_eta"), etaJet1Raw, weight);
    registry.fill(HIST("h_subleadjet_eta"), etaJet2Raw, weight);
    registry.fill(HIST("h_leadjet_phi"), phiLead, weight);
    registry.fill(HIST("h_subleadjet_phi"), phiSub, weight);
    registry.fill(HIST("h2_dijet_detanoflip_dphi"), deltaEtaJetsNoflip, deltaPhiJets, weight);
    registry.fill(HIST("h2_dijet_Asymmetry"), ptSubleadingCorr, ptSubleadingCorr / ptLeadingCorr, weight);
    registry.fill(HIST("h3_dijet_deta_pt"), deltaEtaJets, ptLeadingCorr, ptSubleadingCorr, weight);
    if (multEta1Eta2 > 0)
      registry.fill(HIST("h2_dijet_TimeEtaThan0_pt"), ptLeadingCorr, ptSubleadingCorr, weight);
    else if (multEta1Eta2 < 0)
      registry.fill(HIST("h2_dijet_TimeEtaLess0_pt"), ptLeadingCorr, ptSubleadingCorr, weight);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;
      double hpt = track.pt();
      double heta = inverse * (track.eta()); // Dr.Yang
      double detatot = track.eta() - etaJet1Raw;
      double deta = flip * (track.eta() - etaJet1Raw); // always relative to leadingJet (after flip)
      double dphi = track.phi() - leadingJet.phi();
      dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
      if (hpt > assoHadronPtMaxCut)
        continue;
      registry.fill(HIST("h_jeth_detatot"), detatot, weight);
      registry.fill(HIST("h_jeth_deta"), deta, weight);
      registry.fill(HIST("h_jeth_dphi"), dphi, weight);
      if (doDijetRaa)
        registry.fill(HIST("thn_ljeth_correlations"), ptLeadingCorr, ptSubleadingCorr, kNaN, multEta1Eta2, deltaEtaJetsNoflip, kNaN, kNaN, kNaN, weight);
      else
        registry.fill(HIST("thn_ljeth_correlations"), ptLeadingCorr, ptSubleadingCorr, hpt, multEta1Eta2, deltaEtaJets, heta, deta, dphi, weight);
      if (hpt < assoHadronPtCut) {
        registry.fill(HIST("h2_jeth_deta_dphi"), deta, dphi, weight);
        if (std::abs(deltaEtaJets) >= etaGapup)
          registry.fill(HIST("h2_jeth_physicalcutsup_deta_dphi"), deta, dphi, weight);
        if (std::abs(deltaEtaJets) >= etaGapdw && std::abs(deltaEtaJets) < etaGapup)
          registry.fill(HIST("h2_jeth_physicalcutsmd_deta_dphi"), deta, dphi, weight);
        if (std::abs(deltaEtaJets) < etaGapdw)
          registry.fill(HIST("h2_jeth_physicalcutsdw_deta_dphi"), deta, dphi, weight);
      }
    }
  }

  //.......mixed events leadingjet-hadrons.........................
  template <typename TCollisions, typename TJets, typename TTracks>
  void fillMixLeadingJetHadronHistograms(const TCollisions& collisions, const TJets& jets, const TTracks& tracks, float weight = 1.0)
  {
    using TracksTable = std::decay_t<decltype(tracks)>;
    auto tracksTuple = std::make_tuple(jets, tracks);
    Pair<TCollisions, TJets, TTracks, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, jets1, c2, tracks2] : pairData) {
      weight = doEventWeighted ? c1.weight() : 1.f;
      const float pTHat = 10.f / std::pow(weight, 1.f / pTHatExponent);
      registry.fill(HIST("h_mix_event_stats"), 1);
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), getCentrality(c2)));
      if (!isGoodCollision(c1) || !isGoodCollision(c2))
        continue;

      typename TJets::iterator leadingJet;
      typename TJets::iterator subleadingJet;
      bool hasLeading = false;
      bool hasSubleading = false;
      double ptLeadingCorr = -1.0;
      double ptSubleadingCorr = -1.0;

      for (auto it = jets1.begin(); it != jets1.end(); ++it) {
        const auto& jet = *it;
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
          continue;
        if (!isAcceptedJet<TracksTable>(jet))
          continue;
        if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin)
          continue;

        double ptCorr = jet.pt() - jet.area() * c1.rho();
        if (ptCorr > ptLeadingCorr) {
          if (hasLeading) {
            subleadingJet = leadingJet;
            ptSubleadingCorr = ptLeadingCorr;
            hasSubleading = true;
          }
          leadingJet = it;
          ptLeadingCorr = ptCorr;
          hasLeading = true;
        } else if (ptCorr > ptSubleadingCorr) {
          subleadingJet = it;
          ptSubleadingCorr = ptCorr;
          hasSubleading = true;
        }
      }
      if (!hasLeading || !hasSubleading)
        continue;

      double phiLead = leadingJet.phi();
      double phiSub = subleadingJet.phi();
      double deltaPhiJets = phiLead - phiSub;
      deltaPhiJets = RecoDecay::constrainAngle(deltaPhiJets, -PIHalf);
      if (std::abs(deltaPhiJets) < dijetDphiCut * PI)
        continue;

      registry.fill(HIST("h_mix_event_stats"), 2);
      double etaJet1Raw = leadingJet.eta();
      double etaJet2Raw = subleadingJet.eta();
      double multEta1Eta2 = etaJet1Raw * etaJet2Raw;
      double flip = (etaJet1Raw > etaJet2Raw) ? 1.0 : -1.0;
      double inverse = (etaJet1Raw > 0) ? 1.0 : -1.0; // Dr.Yang suggestion
      double etajet1 = flip * etaJet1Raw;
      double etajet2 = flip * etaJet2Raw;
      double deltaEtaJetsNoflip = etaJet1Raw - etaJet2Raw;
      double deltaEtaJets = etajet1 - etajet2;
      registry.fill(HIST("h_mixleadjet_corrpt"), ptLeadingCorr, weight);
      registry.fill(HIST("h_mixsubleadjet_corrpt"), ptSubleadingCorr, weight);

      if (ptLeadingCorr < leadingjetptMin || ptSubleadingCorr < subleadingjetptMin)
        continue;
      registry.fill(HIST("h_mix_event_stats"), 3);
      registry.fill(HIST("h_mixdijet_pair_counts_cut"), 2);
      registry.fill(HIST("h_mixleadjet_eta"), etaJet1Raw, weight);
      registry.fill(HIST("h_mixsubleadjet_eta"), etaJet2Raw, weight);
      registry.fill(HIST("h2_mixdijet_detanoflip_dphi"), deltaEtaJetsNoflip, deltaPhiJets, weight);
      registry.fill(HIST("h2_mixdijet_Asymmetry"), ptSubleadingCorr, ptSubleadingCorr / ptLeadingCorr, weight);
      registry.fill(HIST("h3_mixdijet_deta_pt"), deltaEtaJets, ptLeadingCorr, ptSubleadingCorr, weight);
      if (multEta1Eta2 > 0)
        registry.fill(HIST("h2_mixdijet_TimeEtaThan0_pt"), ptLeadingCorr, ptSubleadingCorr, weight);
      else if (multEta1Eta2 < 0)
        registry.fill(HIST("h2_mixdijet_TimeEtaLess0_pt"), ptLeadingCorr, ptSubleadingCorr, weight);

      for (auto const& track : tracks2) {
        registry.fill(HIST("h_mix_event_stats"), 4);
        if (!jetderiveddatautilities::selectTrack(track, trackSelection))
          continue;
        registry.fill(HIST("h_mix_event_stats"), 5);
        double hpt = track.pt();
        double heta = inverse * (track.eta()); // Dr.Yang
        double detatot = track.eta() - etaJet1Raw;
        double deta = flip * (track.eta() - etajet1);
        double dphi = track.phi() - phiLead;
        dphi = RecoDecay::constrainAngle(dphi, -PIHalf);

        if (hpt > assoHadronPtMaxCut)
          continue;
        registry.fill(HIST("h_mixjeth_detatot"), detatot, weight);
        registry.fill(HIST("h_mixjeth_deta"), deta, weight);
        registry.fill(HIST("h_mixjeth_dphi"), dphi, weight);
        if (doDijetRaa)
          registry.fill(HIST("thn_mixljeth_correlations"), ptLeadingCorr, ptSubleadingCorr, kNaN, multEta1Eta2, deltaEtaJetsNoflip, kNaN, kNaN, kNaN, weight);
        else
          registry.fill(HIST("thn_mixljeth_correlations"), ptLeadingCorr, ptSubleadingCorr, hpt, multEta1Eta2, deltaEtaJets, heta, deta, dphi, poolBin, weight);
        if (hpt < assoHadronPtCut) {
          registry.fill(HIST("h2_mixjeth_deta_dphi"), deta, dphi, weight);
          if (std::abs(deltaEtaJets) >= etaGapup)
            registry.fill(HIST("h2_mixjeth_physicalcutsup_deta_dphi"), deta, dphi, weight);
          if (std::abs(deltaEtaJets) >= etaGapdw && std::abs(deltaEtaJets) < etaGapup)
            registry.fill(HIST("h2_mixjeth_physicalcutsmd_deta_dphi"), deta, dphi, weight);
          if (std::abs(deltaEtaJets) < etaGapdw)
            registry.fill(HIST("h2_mixjeth_physicalcutsdw_deta_dphi"), deta, dphi, weight);
        }
      }
    }
  }

  //........MCP..leading jet - hadron correlations.....................
  template <typename TmcCollision, typename TJets, typename TParticles>
  void fillMCPLeadingJetHadronHistograms(const TmcCollision& mccollision, const TJets& jets, const TParticles& particles, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    using ParticlesTable = std::decay_t<decltype(particles)>;
    typename TJets::iterator leadingJet;
    typename TJets::iterator subleadingJet;
    bool hasLeading = false;
    bool hasSubleading = false;
    double ptLeadingCorr = -1.0;
    double ptSubleadingCorr = -1.0;

    for (auto it = jets.begin(); it != jets.end(); ++it) {
      const auto& jet = *it;
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (!isAcceptedJet<ParticlesTable>(jet, true))
        continue;
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
        return;

      double ptCorr = jet.pt() - jet.area() * mccollision.rho();
      if (ptCorr > ptLeadingCorr) {
        if (hasLeading) {
          subleadingJet = leadingJet;
          ptSubleadingCorr = ptLeadingCorr;
          hasSubleading = true;
        }
        leadingJet = it;
        ptLeadingCorr = ptCorr;
        hasLeading = true;
      } else if (ptCorr > ptSubleadingCorr) {
        subleadingJet = it;
        ptSubleadingCorr = ptCorr;
        hasSubleading = true;
      }
    }
    if (!hasLeading || !hasSubleading)
      return;
    double deltaPhiJets = RecoDecay::constrainAngle(leadingJet.phi() - subleadingJet.phi(), -PIHalf);
    if (std::abs(deltaPhiJets) < dijetDphiCut * PI)
      return;

    double etaJet1Raw = leadingJet.eta();
    double etaJet2Raw = subleadingJet.eta();
    double multEta1Eta2 = etaJet1Raw * etaJet2Raw;
    double deltaEtaJetsNoflip = etaJet1Raw - etaJet2Raw;
    double inverse = (etaJet1Raw > 0) ? 1.0 : -1.0; // Dr.Yang suggestion
    double flip = (etaJet1Raw > etaJet2Raw) ? 1.0 : -1.0;
    double etajet1 = flip * etaJet1Raw;      // leading jet eta after flip
    double etajet2 = flip * etaJet2Raw;      // subleading jet eta after flip
    double deltaEtaJets = etajet1 - etajet2; // >= 0
    registry.fill(HIST("h_leadjet_pt_part"), leadingJet.pt(), weight);
    registry.fill(HIST("h_subleadjet_pt_part"), subleadingJet.pt(), weight);
    registry.fill(HIST("h_leadjet_corrpt_part"), ptLeadingCorr, weight);
    registry.fill(HIST("h_subleadjet_corrpt_part"), ptSubleadingCorr, weight);

    if (ptLeadingCorr < leadingjetptMin || ptSubleadingCorr < subleadingjetptMin)
      return;
    registry.fill(HIST("h_leadjet_eta_part"), etaJet1Raw, weight);
    registry.fill(HIST("h_subleadjet_eta_part"), etaJet2Raw, weight);
    registry.fill(HIST("h_leadjet_phi_part"), leadingJet.phi(), weight);
    registry.fill(HIST("h_subleadjet_phi_part"), subleadingJet.phi(), weight);
    registry.fill(HIST("h2_dijet_detanoflip_dphi_part"), deltaEtaJetsNoflip, deltaPhiJets, weight);
    registry.fill(HIST("h2_dijet_Asymmetry_part"), ptSubleadingCorr, ptSubleadingCorr / ptLeadingCorr, weight);
    registry.fill(HIST("h3_dijet_deta_pt_part"), deltaEtaJets, ptLeadingCorr, ptSubleadingCorr, weight);
    if (multEta1Eta2 > 0)
      registry.fill(HIST("h2_dijet_TimeEtaThan0_pt_part"), ptLeadingCorr, ptSubleadingCorr, weight);
    else if (multEta1Eta2 < 0)
      registry.fill(HIST("h2_dijet_TimeEtaLess0_pt_part"), ptLeadingCorr, ptSubleadingCorr, weight);

    for (auto const& particle : particles) {
      double hpt = particle.pt();
      double heta = inverse * particle.eta(); // Dr.Yang
      double detatot = particle.eta() - etaJet1Raw;
      double deta = flip * (particle.eta() - etaJet1Raw); // always relative to leadingJet (after flip)
      double dphi = particle.phi() - leadingJet.phi();
      dphi = RecoDecay::constrainAngle(dphi, -PIHalf);
      if (hpt > assoHadronPtMaxCut)
        continue;
      registry.fill(HIST("h_jeth_detatot_part"), detatot, weight);
      registry.fill(HIST("h_jeth_deta_part"), deta, weight);
      registry.fill(HIST("h_jeth_dphi_part"), dphi, weight);
      if (doDijetRaa)
        registry.fill(HIST("thn_ljeth_correlations_part"), ptLeadingCorr, ptSubleadingCorr, kNaN, multEta1Eta2, deltaEtaJetsNoflip, kNaN, kNaN, kNaN, weight);
      else
        registry.fill(HIST("thn_ljeth_correlations_part"), ptLeadingCorr, ptSubleadingCorr, hpt, multEta1Eta2, deltaEtaJets, heta, deta, dphi, weight);
      if (hpt < assoHadronPtCut) {
        registry.fill(HIST("h2_jeth_deta_dphi_part"), deta, dphi);
        if (std::abs(deltaEtaJets) >= etaGapup)
          registry.fill(HIST("h2_jeth_physicalcutsup_deta_dphi_part"), deta, dphi, weight);
        if (std::abs(deltaEtaJets) >= etaGapdw && std::abs(deltaEtaJets) < etaGapup)
          registry.fill(HIST("h2_jeth_physicalcutsmd_deta_dphi_part"), deta, dphi, weight);
        if (std::abs(deltaEtaJets) < etaGapdw)
          registry.fill(HIST("h2_jeth_physicalcutsdw_deta_dphi_part"), deta, dphi, weight);
      }
    }
  }

  //..........MCP..mixed events.........................................
  template <typename TmcCollisions, typename TCollisions, typename TJets, typename TParticles>
  void fillMCPMixLeadingJetHadronHistograms(const TmcCollisions& mccollisions, const TCollisions& collisions, const TJets& jets, const TParticles& particles, float weight = 1.0)
  {
    using ParticlesTable = std::decay_t<decltype(particles)>;
    auto particlesTuple = std::make_tuple(jets, particles);
    Pair<TmcCollisions, TJets, TParticles, BinningTypeMC> pairMCData{corrBinningMC, numberEventsMixed, -1, mccollisions, particlesTuple, &cache};

    for (const auto& [c1, jets1, c2, particles2] : pairMCData) {
      weight = doEventWeighted ? c1.weight() : 1.f;
      const float pTHat = 10.f / std::pow(weight, 1.f / pTHatExponent);
      registry.fill(HIST("h_mixevent_stats_part"), 1);
      // int poolBin = corrBinningMC.getBin(std::make_tuple(c2.posZ(), getMultiplicity(c2)));
      if (!applyMCCollisionCuts(c1, collisions) || !applyMCCollisionCuts(c2, collisions))
        continue;

      typename TJets::iterator leadingJet;
      typename TJets::iterator subleadingJet;
      bool hasLeading = false;
      bool hasSubleading = false;
      double ptLeadingCorr = -1.0;
      double ptSubleadingCorr = -1.0;

      for (auto it = jets1.begin(); it != jets1.end(); ++it) {
        const auto& jet = *it;
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
          continue;
        if (!isAcceptedJet<ParticlesTable>(jet, true))
          continue;
        if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin)
          return;

        double ptCorr = jet.pt() - jet.area() * c1.rho();
        if (ptCorr > ptLeadingCorr) {
          if (hasLeading) {
            subleadingJet = leadingJet;
            ptSubleadingCorr = ptLeadingCorr;
            hasSubleading = true;
          }
          leadingJet = it;
          ptLeadingCorr = ptCorr;
          hasLeading = true;
        } else if (ptCorr > ptSubleadingCorr) {
          subleadingJet = it;
          ptSubleadingCorr = ptCorr;
          hasSubleading = true;
        }
      }
      if (!hasLeading || !hasSubleading)
        continue;

      double deltaPhiJets = leadingJet.phi() - subleadingJet.phi();
      deltaPhiJets = RecoDecay::constrainAngle(deltaPhiJets, -PIHalf);
      if (std::abs(deltaPhiJets) < dijetDphiCut * PI)
        continue;
      registry.fill(HIST("h_mixevent_stats_part"), 2);

      double etaJet1Raw = leadingJet.eta();
      double etaJet2Raw = subleadingJet.eta();
      double multEta1Eta2 = etaJet1Raw * etaJet2Raw;
      double deltaEtaJetsNoflip = etaJet1Raw - etaJet2Raw;
      double inverse = (etaJet1Raw > 0) ? 1.0 : -1.0; // Dr.Yang suggestion
      double flip = (etaJet1Raw > etaJet2Raw) ? 1.0 : -1.0;
      double etajet1 = flip * etaJet1Raw;
      double etajet2 = flip * etaJet2Raw;
      double deltaEtaJets = etajet1 - etajet2;
      registry.fill(HIST("h_mixleadjet_corrpt_part"), ptLeadingCorr, weight);
      registry.fill(HIST("h_mixsubleadjet_corrpt_part"), ptSubleadingCorr, weight);

      if (ptLeadingCorr < leadingjetptMin || ptSubleadingCorr < subleadingjetptMin)
        return;
      registry.fill(HIST("h_mixevent_stats_part"), 3);
      registry.fill(HIST("h_mixdijet_pair_counts_cut_part"), 1);
      registry.fill(HIST("h_mixleadjet_eta_part"), etaJet1Raw, weight);
      registry.fill(HIST("h_mixsubleadjet_eta_part"), etaJet2Raw, weight);
      registry.fill(HIST("h2_mixdijet_detanoflip_dphi_part"), deltaEtaJetsNoflip, deltaPhiJets, weight);
      registry.fill(HIST("h2_mixdijet_Asymmetry_part"), ptSubleadingCorr, ptSubleadingCorr / ptLeadingCorr, weight);
      registry.fill(HIST("h3_mixdijet_deta_pt_part"), deltaEtaJets, ptLeadingCorr, ptSubleadingCorr, weight);
      if (multEta1Eta2 > 0)
        registry.fill(HIST("h2_mixdijet_TimeEtaThan0_pt_part"), ptLeadingCorr, ptSubleadingCorr, weight);
      else if (multEta1Eta2 < 0)
        registry.fill(HIST("h2_mixdijet_TimeEtaLess0_pt_part"), ptLeadingCorr, ptSubleadingCorr, weight);

      for (auto const& particle : particles2) {
        registry.fill(HIST("h_mixevent_stats_part"), 4);
        double hpt = particle.pt();
        double heta = inverse * particle.eta(); // Dr.Yang
        double detatot = particle.eta() - etaJet1Raw;
        double deta = flip * (particle.eta() - etaJet1Raw);
        double dphi = particle.phi() - leadingJet.phi();
        dphi = RecoDecay::constrainAngle(dphi, -PIHalf);

        if (hpt > assoHadronPtMaxCut)
          continue;
        registry.fill(HIST("h_mixevent_stats_part"), 5);
        registry.fill(HIST("h_mixjeth_detatot_part"), detatot, weight);
        registry.fill(HIST("h_mixjeth_deta_part"), deta, weight);
        registry.fill(HIST("h_mixjeth_dphi_part"), dphi, weight);
        if (doDijetRaa)
          registry.fill(HIST("thn_mixljeth_correlations_part"), ptLeadingCorr, ptSubleadingCorr, kNaN, multEta1Eta2, deltaEtaJetsNoflip, kNaN, kNaN, kNaN, weight);
        else
          registry.fill(HIST("thn_mixljeth_correlations_part"), ptLeadingCorr, ptSubleadingCorr, hpt, multEta1Eta2, deltaEtaJets, heta, deta, dphi, weight);
        if (hpt < assoHadronPtCut) {
          registry.fill(HIST("h2_mixjeth_deta_dphi_part"), deta, dphi);
          if (std::abs(deltaEtaJets) >= etaGapup)
            registry.fill(HIST("h2_mixjeth_physicalcutsup_deta_dphi_part"), deta, dphi, weight);
          if (std::abs(deltaEtaJets) >= etaGapdw && std::abs(deltaEtaJets) < etaGapup)
            registry.fill(HIST("h2_mixjeth_physicalcutsmd_deta_dphi_part"), deta, dphi, weight);
          if (std::abs(deltaEtaJets) < etaGapdw)
            registry.fill(HIST("h2_mixjeth_physicalcutsdw_deta_dphi_part"), deta, dphi, weight);
        }
      }
    }
  }

  // ==========================================================
  //.............process staring...............................
  // ==========================================================
  void processCollisionsQCData(aod::JetCollision const& collision,
                               FilterJetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents))
      return;
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    if (!isGoodCollision(collision))
      return;
    registry.fill(HIST("h_collisions"), 3.5);
    registry.fill(HIST("h2_centrality_occupancy"), getCentrality(collision), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ());
    registry.fill(HIST("h_collisions_multFT0"), getMultiplicity(collision)); // collision.MultFT0M()

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;
      fillTrackHistograms(track);
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processCollisionsQCData, "QC of collisions and tracks for Data", true);

  void processSpectraAreaSubData(FilterCollision const& collision,
                                 CorrChargedJets const& jets,
                                 aod::JetTracks const&)
  {
    if (!isGoodCollision(collision))
      return;
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
      fillJetAreaSubHistograms(jet, collision.rho());
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processSpectraAreaSubData, "jet spectra without and with rho-area subtraction for Data", false);

  void processLeadingJetHadron(FilterCollision const& collision,
                               CorrChargedJets const& jets,
                               FilterJetTracks const& tracks)
  {
    if (!isGoodCollision(collision))
      return;
    fillLeadingJetHadronHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processLeadingJetHadron, "ame event subleading jet-h for Data", false);

  void processMixLeadingJetHadron(FilterCollisions const& collisions,
                                  CorrChargedJets const& jets,
                                  FilterJetTracks const& tracks)
  {
    if (collisions.size() == 0)
      return;
    fillMixLeadingJetHadronHistograms(collisions, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixLeadingJetHadron, "mixed event leading jet-h for Data", false);

  void processJetHadron(FilterCollision const& collision,
                        CorrChargedJets const& jets,
                        FilterJetTracks const& tracks)
  {
    if (!isGoodCollision(collision))
      return;
    fillJetHadronHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processJetHadron, "same event jet-h for Data", false);

  void processMixJetHadron(FilterCollisions const& collisions,
                           CorrChargedJets const& jets,
                           FilterJetTracks const& tracks)
  {
    if (collisions.size() == 0)
      return;
    fillMixJetHadronHistograms(collisions, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixJetHadron, "mixed event jet-h for Data", false);

  //...HF jet correlations....................
  void processHFJetCorrelation(FilterCollision const& collision,
                               CorrChargedJets const& jets,
                               aod::CandidatesD0Data const& candidates)
  {
    if (!isGoodCollision(collision))
      return;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      registry.fill(HIST("h_d0jet_pt"), jet.pt());
      registry.fill(HIST("h_d0jet_corrpt"), jet.pt() - collision.rho() * jet.area());
      registry.fill(HIST("h_d0jet_eta"), jet.eta());
      registry.fill(HIST("h_d0jet_phi"), jet.phi());
    }
    for (const auto& candidate : candidates) {
      registry.fill(HIST("h_d0_mass"), candidate.m());
      registry.fill(HIST("h_d0_pt"), candidate.pt());
      registry.fill(HIST("h_d0_eta"), candidate.eta());
      registry.fill(HIST("h_d0_phi"), candidate.phi());
      for (const auto& jet : jets) {
        double deltaeta = candidate.eta() - jet.eta();
        double deltaphi = candidate.phi() - jet.phi();
        deltaphi = RecoDecay::constrainAngle(deltaphi, -PIHalf);
        registry.fill(HIST("h2_d0jet_detadphi"), deltaeta, deltaphi);
      }
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processHFJetCorrelation, "D0-jet for Data", false);

  //........MCD..................................................
  void processCollisionsQCMCD(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                              FilterJetTracks const& tracks)
  {
    const float eventWeight = doEventWeighted ? collision.weight() : 1.f;
    if (!collision.has_mcCollision()) {
      registry.fill(HIST("h_fakecollisions"), 0.5);
    }
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents))
      return;
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange())
      return;
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);

    if (!isGoodCollision(collision))
      return;
    registry.fill(HIST("h_collisions"), 3.5);
    registry.fill(HIST("h_collisions_weighted"), 3.5, eventWeight);
    registry.fill(HIST("h2_centrality_occupancy"), getCentrality(collision), collision.trackOccupancyInTimeRange(), eventWeight);
    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ(), eventWeight);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;
      fillTrackHistograms(track, eventWeight);
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processCollisionsQCMCD, "QC of collisions and tracks for MCD", false);

  void processSpectraAreaSubMCD(FilterCollision const& collision,
                                CorrChargedMCDJets const& jets,
                                aod::JetTracks const&)
  {
    const float eventWeight = doEventWeighted ? collision.weight() : 1.f;
    if (!isGoodCollision(collision))
      return;
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet, eventWeight);
      fillJetAreaSubHistograms(jet, collision.rho(), eventWeight);
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processSpectraAreaSubMCD, "jet spectra with rho-area subtraction for MCD", false);

  void processLeadinJetHadronMCD(FilterCollision const& collision,
                                 CorrChargedMCDJets const& jets,
                                 FilterJetTracks const& tracks)
  {
    const float eventWeight = doEventWeighted ? collision.weight() : 1.f;
    if (!isGoodCollision(collision)) {
      return;
    }
    fillLeadingJetHadronHistograms(collision, jets, tracks, eventWeight);
  }
  PROCESS_SWITCH(ChargedJetHadron, processLeadinJetHadronMCD, "same event leading jet-hadron correlations for MCD", false);

  void processMixLeadinJetHadronMCD(FilterCollisions const& collisions,
                                    CorrChargedMCDJets const& jets,
                                    FilterJetTracks const& tracks)
  {
    if (collisions.size() == 0)
      return;
    fillMixLeadingJetHadronHistograms(collisions, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixLeadinJetHadronMCD, "mixed event leading jet-hadron correlations for MCD", false);

  void processJetHadronMCD(FilterCollision const& collision,
                           CorrChargedMCDJets const& jets,
                           FilterJetTracks const& tracks)
  {
    const float eventWeight = doEventWeighted ? collision.weight() : 1.f;
    if (!isGoodCollision(collision))
      return;
    fillJetHadronHistograms(collision, jets, tracks, eventWeight);
  }
  PROCESS_SWITCH(ChargedJetHadron, processJetHadronMCD, "same event jet-hadron correlations for MCD", false);

  void processMixJetHadronMCD(FilterCollisions const& collisions,
                              CorrChargedMCDJets const& jets,
                              FilterJetTracks const& tracks)
  {
    if (collisions.size() == 0)
      return;
    fillMixJetHadronHistograms(collisions, jets, tracks);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixJetHadronMCD, "mixed event jet-hadron correlations for MCD", false);

  //........MCP..................................................
  void processCollisionsQCMCP(aod::JetMcCollision const& mccollision,
                              soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                              soa::Filtered<aod::JetParticles> const& particles)
  {
    const float eventWeight = doEventWeighted ? mccollision.weight() : 1.f;
    registry.fill(HIST("h_mcColl_counts"), 0.5);

    if (std::abs(mccollision.posZ()) > vertexZCut)
      return;
    registry.fill(HIST("h_mcColl_counts"), 1.5);

    if (collisions.size() < 1)
      return;
    registry.fill(HIST("h_mcColl_counts"), 2.5);

    if (acceptSplitCollisions.value == NonSplitOnly && collisions.size() > 1)
      return;
    registry.fill(HIST("h_mcColl_counts"), 3.5);

    int nGood = 0;
    for (auto const& collision : collisions) {
      if (isGoodCollision(collision)) {
        registry.fill(HIST("h_mcdColl_mult"), getMultiplicity(collision), eventWeight);
        nGood++;
      }
    }
    if (nGood == 0)
      return;
    registry.fill(HIST("h_mcColl_counts"), 4.5);
    registry.fill(HIST("h_mcpColl_Zvertex"), mccollision.posZ(), eventWeight);
    registry.fill(HIST("h_mcpColl_centrality"), mccollision.centFT0M(), eventWeight);
    registry.fill(HIST("h_mcpColl_mult"), getMultiplicity(mccollision), eventWeight);
    for (auto const& particle : particles) {
      fillParticleHistograms(particle, eventWeight);
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processCollisionsQCMCP, "QC of collisions and particles for MCP", false);

  void processSpectraAreaSubMCP(FilterMcpCollision const& mccollision,
                                soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                CorrChargedMCPJets const& jets,
                                soa::Filtered<aod::JetParticles> const&)
  {
    bool mcLevelIsParticleLevel = true;
    const float eventWeight = doEventWeighted ? mccollision.weight() : 1.f;
    if (!applyMCCollisionCuts(mccollision, collisions))
      return;
    registry.fill(HIST("h_mcColl_rho"), mccollision.rho(), eventWeight);
    if (doEventWeighted) {
      registry.fill(HIST("h_mcColl_counts_weight"), 1.5);
      registry.fill(HIST("h_mcColl_counts_weight"), 2.5, mccollision.xsectGen());
      registry.fill(HIST("h_mcColl_counts_weight"), 3.5, eventWeight);
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
        continue;
      }
      fillMCPHistograms(jet, eventWeight);
      fillMCPAreaSubHistograms(jet, mccollision.rho(), eventWeight);
    }
  }
  PROCESS_SWITCH(ChargedJetHadron, processSpectraAreaSubMCP, "jet spectra without and with UE subtraction of area-based for MCP", false);

  void processLeadingJetHadronMCP(FilterMcpCollision const& mccollision,
                                  soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                  CorrChargedMCPJets const& jets,
                                  soa::Filtered<aod::JetParticles> const& particles)
  {
    const float eventWeight = doEventWeighted ? mccollision.weight() : 1.f;
    if (!applyMCCollisionCuts(mccollision, collisions))
      return;

    fillMCPLeadingJetHadronHistograms(mccollision, jets, particles, eventWeight);
  }
  PROCESS_SWITCH(ChargedJetHadron, processLeadingJetHadronMCP, "same event leading jet-hadron for MCP", false);

  void processMixLeadingJetHadronMCP(FilterMcpCollisions const& mccollisions,
                                     soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                     CorrChargedMCPJets const& jets,
                                     soa::Filtered<aod::JetParticles> const& particles)
  {
    if (mccollisions.size() < 1 || collisions.size() < 1)
      return;
    fillMCPMixLeadingJetHadronHistograms(mccollisions, collisions, jets, particles);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixLeadingJetHadronMCP, "mixed event leading jet-hadron for MCP", false);

  void processJetHadronMCP(FilterMcpCollision const& mccollision,
                           soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                           CorrChargedMCPJets const& jets,
                           soa::Filtered<aod::JetParticles> const& particles)
  {
    const float eventWeight = doEventWeighted ? mccollision.weight() : 1.f;
    if (!applyMCCollisionCuts(mccollision, collisions))
      return;

    fillMCPJetHadronHistograms(mccollision, jets, particles, eventWeight);
  }
  PROCESS_SWITCH(ChargedJetHadron, processJetHadronMCP, "same event jet-hadron for MCP", false);

  void processMixJetHadronMCP(FilterMcpCollisions const& mccollisions,
                              soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                              CorrChargedMCPJets const& jets,
                              soa::Filtered<aod::JetParticles> const& particles)
  {
    if (mccollisions.size() < 1 || collisions.size() < 1)
      return;
    fillMCPMixJetHadronHistograms(mccollisions, collisions, jets, particles);
  }
  PROCESS_SWITCH(ChargedJetHadron, processMixJetHadronMCP, "mixed event jet-hadron for MCP", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ChargedJetHadron>(cfgc)};
}
