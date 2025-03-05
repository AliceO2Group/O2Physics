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

#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH3F.h>
#include <cmath>
#include <TRandom3.h>
#include <TF1.h>
#include <algorithm>
// o2Physics includes.
#include "CommonConstants/MathConstants.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "Framework/runDataProcessing.h"

#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/ctpRateFetcher.h"

//< evt pln .h >//
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
//< evt pln .h | end >//

// o2 includes.
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetChargedV2 {
  HistogramRegistry registry;
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
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
  Configurable<float> jetRadius{"jetRadius", 0.2, "jet resolution parameters"};

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
  //=====================< evt pln | end >=====================//

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  double evtnum = 0;
  double accptTrack = 0;
  double fitTrack = 0;
  TH1F* hPtsumSumptFit = nullptr;

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

    auto jetPtTemp = 0.0;
    jetPtBins.push_back(jetPtTemp);
    jetPtBinsRhoAreaSub.push_back(jetPtTemp);
    while (jetPtTemp < jetPtMax) {
      if (jetPtTemp < 100.0) {
        jetPtTemp += 1.0;
        jetPtBins.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(jetPtTemp);
        jetPtBinsRhoAreaSub.push_back(-jetPtTemp);
      } else if (jetPtTemp < 200.0) {
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

    AxisSpec jetPtAxis = {jetPtBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {jetPtBinsRhoAreaSub, "#it{p}_{T} (GeV/#it{c})"};

    AxisSpec axisPt = {40, 0.0, 4.0};
    AxisSpec axisEta = {32, -0.8, 0.8};
    AxisSpec axixCent = {20, 0, 100};
    AxisSpec axisChID = {220, 0, 220};

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    //< \sigma p_T at local rho test plot >
    registry.add("h_accept_Track", "all and accept track;Track;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
    registry.add("h_accept_Track_init", "all and accept track;Track;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
    registry.add("h_accept_Track_Fit", "all and accept track;Track;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});

    registry.add("h_ptsum_collnum", "ptsum collnum;collnum;entries", {HistType::kTH1F, {{40, 0.0, 40}}});
    registry.add("h_ptsum_sumpt", "jet sumpt;sum p_{T};entries", {HistType::kTH1F, {{40, 0., o2::constants::math::TwoPI}}});
    registry.add("h2_phi_track_eta", "phi vs track eta; #eta (GeV/#it{c}); #varphi", {HistType::kTH2F, {{100, -1.0, 1.0}, {40, 0., o2::constants::math::TwoPI}}});
    registry.add("h2_phi_track_pt", "phi vs track pT; #it{p}_{T,track} (GeV/#it{c}); #varphi", {HistType::kTH2F, {{200, 0., 200.}, {40, 0., o2::constants::math::TwoPI}}});
    registry.add("h2_centrality_phi_w_pt", "centrality vs jet #varphi; centrality; entries", {HistType::kTH2F, {{100, 0.0, 100.0}, {40, 0., o2::constants::math::TwoPI}}});
    registry.add("h2_evtnum_phi_w_pt", "eventNumber vs jet #varphi; #eventNumber; entries", {HistType::kTH2F, {{1000, 0.0, 1000}, {40, 0., o2::constants::math::TwoPI}}});

    //< fit quality >//
    registry.add("h_PvalueCDF_CombinFit", "cDF #chi^{2}; entries", {HistType::kTH1F, {{50, 0, 1}}});
    registry.add("h2_PvalueCDFCent_CombinFit", "p-value cDF vs centrality; centrality; p-value", {HistType::kTH2F, {{100, 0, 100}, {40, 0, 1}}});
    registry.add("h2_Chi2Cent_CombinFit", "Chi2 vs centrality; centrality; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 5}}});
    registry.add("h2_PChi2_CombinFit", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});

    registry.add("Thn_PChi2_CombinFitCent", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTHnSparseF, {{100, 0.0, 100.0}, {100, 0, 1}, {100, 0, 5}}});
    registry.add("h2_PChi2_CombinFitA", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});
    registry.add("h2_PChi2_CombinFitB", "p-value vs #tilde{#chi^{2}}; p-value; #tilde{#chi^{2}}", {HistType::kTH2F, {{100, 0, 1}, {100, 0, 5}}});

    registry.add("h_evtnum_centrlity", "eventNumber vs centrality ; #eventNumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
    registry.add("h_evtnum_NTrk", "eventNumber vs Number of Track ; #eventNumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

    registry.add("Thn_evtnum_phi_centrality", "eventNumber vs jet #varphi; #eventNumber; entries", {HistType::kTHnSparseF, {{1000, 0.0, 1000}, {40, 0., o2::constants::math::TwoPI}, {100, 0.0, 100.0}}});

    registry.add("h2_evt_fitpara", "event vs fit parameter; evtnum; parameter", {HistType::kTH2F, {cfgAxisEvtfit, {5, 0., 5}}});

    registry.add("h_v2obs_centrality", "fitparameter v2obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});
    registry.add("h_v3obs_centrality", "fitparameter v3obs vs centrality ; #centrality", {HistType::kTProfile, {cfgAxisVnCent}});

    registry.add("h_fitparaRho_evtnum", "fitparameter #rho_{0} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
    registry.add("h_fitparaPsi2_evtnum", "fitparameter #Psi_{2} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
    registry.add("h_fitparaPsi3_evtnum", "fitparameter #Psi_{3} vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
    registry.add("h_fitparav2obs_evtnum", "fitparameter v2obs vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});
    registry.add("h_fitparav3obs_evtnum", "fitparameter v3obs vs evtnum ; #eventnumber", {HistType::kTH1F, {{1000, 0.0, 1000}}});

    registry.add("h2_phi_rhophi", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho_{ch}(#varphi) ", {HistType::kTH2F, {{40, 0., o2::constants::math::TwoPI}, {210, -10.0, 200.0}}});
    registry.add("h2_phi_rholocal", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho(#varphi) ", {HistType::kTH2F, {{40, 0., o2::constants::math::TwoPI}, {210, -10.0, 200.0}}});
    registry.add("h2_phi_rhomedian", "#varphi vs #rho(#varphi); #varphi - #Psi_{EP,2};  #rho_{median} ", {HistType::kTH2F, {{40, 0., o2::constants::math::TwoPI}, {210, -10.0, 200.0}}});

    registry.add("h2_fitParaZero_cent", "#varphi vs #rho(#varphi); #cent;  #fitParameter[0] ", {HistType::kTH2F, {{100, 0., 100}, {210, -10.0, 200.0}}});
    registry.add("h2_phi_rholocal_cent", "#varphi vs #rho(#varphi); #cent;  #rho(#varphi) ", {HistType::kTH2F, {{100, 0., 100}, {210, -10.0, 200.0}}});
    registry.add("h2_phi_rhomedian_cent", "#varphi vs #rho(#varphi); #cent;  #rho_{median} ", {HistType::kTH2F, {{100, 0., 100}, {210, -10.0, 200.0}}});

    registry.add("h3_centrality_rhovsphi_phi", "centrality; #rho(#varphi); #Delta#varphi_{jet}", {HistType::kTH3F, {{120, -10.0, 110.0}, {200, 0.0, 200.0}, {40, 0., o2::constants::math::TwoPI}}});
    //< \sigma p_T at local rho test plot | end >

    registry.add("h2_centrality_track_pt", "centrality vs track pT; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
    registry.add("h2_centrality_track_eta", "centrality vs track #eta; centrality; #eta_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
    registry.add("h2_centrality_track_phi", "centrality vs track #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});

    registry.add("h_jet_pt_rhoareasubtracted", "jet pT rhoareasubtracted;#it{p}_{T,jet} (GeV/#it{c}); entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_rholocal", "jet pT rholocal;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

    registry.add("h2_centrality_jet_pt_rhoareasubtracted", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, jetPtAxisRhoAreaSub}});

    registry.add("leadJetPt", "leadJet Pt ", {HistType::kTH1F, {{200, 0., 200.0}}});
    registry.add("leadJetPhi", "leadJet constituent #phi ", {HistType::kTH1F, {{80, -1.0, 7.}}});
    registry.add("leadJetEta", "leadJet constituent #eta ", {HistType::kTH1F, {{100, -1.0, 1.0}}});

    //< RC test plots >//
    registry.add("h3_centrality_deltapT_RandomCornPhi_rhorandomconewithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
    registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphi", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
    registry.add("h3_centrality_deltapT_RandomCornPhi_RCprocess_rhorandomconewithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});

    registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}(#varphi); #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
    registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutoneleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}(#varphi); #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
    registry.add("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithouttwoleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}(#varphi); #Delta#varphi_{jet}", {HistType::kTH3F, {{100, 0.0, 100.0}, {400, -200.0, 200.0}, {100, 0., o2::constants::math::TwoPI}}});
    //< bkg sub plot | end >//
    //< median rho >//
    registry.add("h_jet_pt_in_plane_v2", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v2", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_in_plane_v3", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v3", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

    registry.add("h2_centrality_jet_pt_in_plane_v2", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v2", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_in_plane_v3", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v3", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    //< rho(phi) >//
    registry.add("h_jet_pt_inclusive_v2_rho", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_in_plane_v2_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v2_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_in_plane_v3_rho", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v3_rho", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

    registry.add("h2_centrality_jet_pt_in_plane_v2_rho", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v2_rho", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_in_plane_v3_rho", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v3_rho", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});

    //=====================< evt pln plot >=====================//
    AxisSpec axisCent{cfgAxisCent, "centrality"};
    AxisSpec axisQvec{cfgAxisQvec, "Q"};
    AxisSpec axisQvecF{cfgAxisQvecF, "Q"};

    AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("histCent", "Centrality TrkProcess", HistType::kTH1F, {axisCent});

    for (uint i = 0; i < cfgnMods->size(); i++) {
      histosQA.add(Form("histQvecUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecRectrV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecTwistV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

      histosQA.add(Form("histEvtPlUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlRectrV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlTwistV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    }
    //=====================< evt pln plot | end >=====================//
  }

  Preslice<aod::ChargedJets> jetsPerJCollision = o2::aod::jet::collisionId;
  Preslice<aod::JetTracks> tracksPerJCollision = o2::aod::jtrack::collisionId;

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
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

      for (const auto& constituent : jet.template tracks_as<T>()) {
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

  template <typename T, typename U>
  bool trackIsInJet(T const& track, U const& jet)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  void fillLeadingJetQA(double leadingJetPt, double leadingJetPhi, double leadingJetEta)
  {
    registry.fill(HIST("leadJetPt"), leadingJetPt);
    registry.fill(HIST("leadJetPhi"), leadingJetPhi);
    registry.fill(HIST("leadJetEta"), leadingJetEta);
  }

  double chiSquareCDF(int nDF, double x)
  {
    return TMath::Gamma(nDF / 2., x / 2.);
  }

  void processInOutJetV2(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                         aod::JetTracks const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    histosQA.fill(HIST("histCent"), collision.cent());
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == 2) {
        if (collision.qvecAmp()[detId] > 1e-8) {
          histosQA.fill(HIST("histQvecUncorV2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
          histosQA.fill(HIST("histQvecRectrV2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
          histosQA.fill(HIST("histQvecTwistV2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
          histosQA.fill(HIST("histQvecFinalV2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

          histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());
        }
      } else if (nmode == 3) {
        histosQA.fill(HIST("histQvecUncorV3"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.cent());
        histosQA.fill(HIST("histQvecRectrV3"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.cent());
        histosQA.fill(HIST("histQvecTwistV3"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.cent());
        histosQA.fill(HIST("histQvecFinalV3"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.cent());

        histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlRectrV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlTwistV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.cent());
        histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.cent());
      }

      if (nmode == 2) {
        double phiMinusPsi2;
        if (collision.qvecAmp()[detId] < 1e-8) {
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
          registry.fill(HIST("h2_centrality_jet_pt_rhoareasubtracted"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);

          phiMinusPsi2 = jet.phi() - ep2;

          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= 7 * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= 3 * o2::constants::math::PIQuarter && phiMinusPsi2 < 5 * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      } else if (nmode == 3) {
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

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= 7 * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= 3 * o2::constants::math::PIQuarter && phiMinusPsi3 < 5 * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v3"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v3"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetChargedV2, processInOutJetV2, "Jet V2 in and out of plane", true);

  void processSigmaPt(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                      soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                      aod::JetTracks const& tracks)
  {
    // for (const auto& collision : collisions) {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    for (const auto& jet : jets) {
      if (jet.pt() > leadingJetPt) {
        leadingJetPt = jet.pt();
        leadingJetEta = jet.eta();
        leadingJetPhi = jet.phi();
      }
    }
    fillLeadingJetQA(leadingJetPt, leadingJetPhi, leadingJetEta);

    int nTrk = 0;
    if (jets.size() > 0) {
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= 0.2 && track.pt() <= 5.) {
          registry.fill(HIST("h_accept_Track"), 2.5);
          nTrk += 1;
        }
      }
      registry.fill(HIST("h_evtnum_NTrk"), evtnum, nTrk);
    }

    if (nTrk <= 0) {
      return;
    }

    hPtsumSumptFit = new TH1F("h_ptsum_sumpt_fit", "h_ptsum_sumpt fit use", TMath::CeilNint(std::sqrt(nTrk)), 0., o2::constants::math::TwoPI);

    if (jets.size() > 0) {
      for (auto const& trackfit : tracks) {
        registry.fill(HIST("h_accept_Track"), 0.5);
        if (jetderiveddatautilities::selectTrack(trackfit, trackSelection) && (std::fabs(trackfit.eta() - leadingJetEta) > jetRadius) && trackfit.pt() >= 0.2 && trackfit.pt() <= 5.) {
          registry.fill(HIST("h_accept_Track_Fit"), 0.5);
          fitTrack += 1;
        }
      }

      for (auto const& track : tracks) {
        registry.fill(HIST("h_accept_Track"), 1.5);
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && (std::fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= 0.2 && track.pt() <= 5.) {
          accptTrack += 1;
          registry.fill(HIST("h_accept_Track"), 2.5);
          hPtsumSumptFit->Fill(track.phi(), track.pt());
          registry.fill(HIST("h2_phi_track_pt"), track.pt(), track.phi());
          registry.fill(HIST("h2_phi_track_eta"), track.eta(), track.phi());
          registry.fill(HIST("h_ptsum_sumpt"), track.phi(), track.pt());
          registry.fill(HIST("h2_centrality_phi_w_pt"), collision.centrality(), track.phi(), track.pt());
          registry.fill(HIST("h2_evtnum_phi_w_pt"), evtnum, track.phi(), track.pt());
          registry.fill(HIST("Thn_evtnum_phi_centrality"), evtnum, track.phi(), collision.centrality());
          registry.fill(HIST("h_accept_Track_init"), accptTrack);
          registry.fill(HIST("h_accept_Track_Fit"), 1.5);
        }
      }
    }
    registry.fill(HIST("h_ptsum_collnum"), 0.5);

    double ep2 = 0.;
    double ep3 = 0.;
    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == 2) {
        if (collision.qvecAmp()[detId] > 1e-8) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      } else if (nmode == 3) {
        if (collision.qvecAmp()[detId] > 1e-8) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
        }
      }
    }

    TF1* fFitModulationV2v3 = 0x0;
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
    for (int i = 1; i <= 5; i++) {
      registry.fill(HIST("h2_evt_fitpara"), evtnum, i - 0.5, temppara[i - 1]);
    }

    registry.fill(HIST("h_fitparaRho_evtnum"), evtnum, temppara[0]);
    registry.fill(HIST("h_fitparav2obs_evtnum"), evtnum, temppara[1]);
    registry.fill(HIST("h_fitparaPsi2_evtnum"), evtnum, temppara[2]);
    registry.fill(HIST("h_fitparav3obs_evtnum"), evtnum, temppara[3]);
    registry.fill(HIST("h_fitparaPsi3_evtnum"), evtnum, temppara[4]);

    registry.fill(HIST("h_v2obs_centrality"), collision.centrality(), temppara[1]);
    registry.fill(HIST("h_v3obs_centrality"), collision.centrality(), temppara[3]);
    registry.fill(HIST("h_evtnum_centrlity"), evtnum, collision.centrality());

    if (temppara[0] == 0) {
      return;
    }
    registry.fill(HIST("h2_fitParaZero_cent"), collision.centrality(), temppara[0], 1.0);

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

    registry.fill(HIST("h_PvalueCDF_CombinFit"), cDF);
    registry.fill(HIST("h2_PvalueCDFCent_CombinFit"), collision.centrality(), cDF);
    registry.fill(HIST("h2_Chi2Cent_CombinFit"), collision.centrality(), chiSqr / (static_cast<float>(nDF)));
    registry.fill(HIST("h2_PChi2_CombinFit"), cDF, chiSqr / (static_cast<float>(nDF)));
    registry.fill(HIST("Thn_PChi2_CombinFitCent"), collision.centrality(), cDF, chiSqr / (static_cast<float>(nDF)));
    double evtcent = collision.centrality();
    if (evtcent >= 0 && evtcent <= 5) {
      registry.fill(HIST("h2_PChi2_CombinFitA"), cDF, chiSqr / (static_cast<float>(nDF)));
    } else if (evtcent >= 30 && evtcent <= 50) {
      registry.fill(HIST("h2_PChi2_CombinFitB"), cDF, chiSqr / (static_cast<float>(nDF)));
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
        registry.fill(HIST("h2_phi_rholocal_cent"), collision.centrality(), rholocal, 1.0);
        registry.fill(HIST("h2_phi_rhomedian_cent"), collision.centrality(), collision.rho(), 1.0);

        if (nmode == 2) {
          registry.fill(HIST("h_jet_pt_rholocal"), jet.pt() - (rholocal * jet.area()), 1.0);

          double phiMinusPsi2;
          if (collision.qvecAmp()[detId] < 1e-8) {
            continue;
          }
          phiMinusPsi2 = jet.phi() - ep2;

          registry.fill(HIST("h2_phi_rholocal"), jet.phi() - ep2, rholocal, 1.0);
          registry.fill(HIST("h_jet_pt_inclusive_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);

          if ((phiMinusPsi2 < o2::constants::math::PIQuarter) || (phiMinusPsi2 >= 7 * o2::constants::math::PIQuarter) || (phiMinusPsi2 >= 3 * o2::constants::math::PIQuarter && phiMinusPsi2 < 5 * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2_rho"), collision.centrality(), jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v2_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2_rho"), collision.centrality(), jet.pt() - (rholocal * jet.area()), 1.0);
          }
        } else if (nmode == 3) {
          double phiMinusPsi3;
          if (collision.qvecAmp()[detId] < 1e-8) {
            continue;
          }
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
          phiMinusPsi3 = jet.phi() - ep3;

          if ((phiMinusPsi3 < o2::constants::math::PIQuarter) || (phiMinusPsi3 >= 7 * o2::constants::math::PIQuarter) || (phiMinusPsi3 >= 3 * o2::constants::math::PIQuarter && phiMinusPsi3 < 5 * o2::constants::math::PIQuarter)) {
            registry.fill(HIST("h_jet_pt_in_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_in_plane_v3_rho"), collision.centrality(), jet.pt() - (rholocal * jet.area()), 1.0);
          } else {
            registry.fill(HIST("h_jet_pt_out_of_plane_v3_rho"), jet.pt() - (rholocal * jet.area()), 1.0);
            registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v3_rho"), collision.centrality(), jet.pt() - (rholocal * jet.area()), 1.0);
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
      if (nmode == 2) {
        double rcPhiPsi2;
        rcPhiPsi2 = randomConePhi - ep2;

        for (auto const& jet : jets) {
          registry.fill(HIST("h2_phi_rhophi"), jet.phi() - ep2, fFitModulationV2v3->Eval(jet.phi()), 1.0);
          registry.fill(HIST("h3_centrality_rhovsphi_phi"), collision.centrality(), fFitModulationV2v3->Eval(jet.phi()), jet.phi() - ep2);
        }

        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphi"), collision.centrality(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);

        // removing the leading jet from the random cone
        if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
          float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
          float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

          bool jetWasInCone = false;
          while (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR) {
            jetWasInCone = true;
            randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
            randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
            dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
          }
          if (jetWasInCone) {
            randomConePt = 0.0;
            for (auto const& track : tracks) {
              if (jetderiveddatautilities::selectTrack(track, trackSelection)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
                float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
                float dEta = track.eta() - randomConeEta;
                if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
                  randomConePt += track.pt();
                }
              }
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutleadingjet"), collision.centrality(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), rcPhiPsi2, 1.0);

        // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles, removing tracks from 2 leading jets
        double randomConePtWithoutOneLeadJet = 0;
        double randomConePtWithoutTwoLeadJet = 0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              if (!trackIsInJet(track, jets.iteratorAt(0))) {
                randomConePtWithoutOneLeadJet += track.pt();
                if (!trackIsInJet(track, jets.iteratorAt(1))) {
                  randomConePtWithoutTwoLeadJet += track.pt();
                }
              }
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithoutoneleadingjet"), collision.centrality(), randomConePtWithoutOneLeadJet - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_localrhovsphiwithouttwoleadingjet"), collision.centrality(), randomConePtWithoutTwoLeadJet - o2::constants::math::PI * randomConeR * randomConeR * rholocalRC, rcPhiPsi2, 1.0);
      } else if (nmode == 3) {
        continue;
      }
    }
    delete hPtsumSumptFit;
    evtnum += 1;
  }
  PROCESS_SWITCH(JetChargedV2, processSigmaPt, "Sigma pT and bkg as fcn of phi", true);

  void processRandomConeDataV2(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                               soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                               soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    for (uint i = 0; i < cfgnMods->size(); i++) {
      TRandom3 randomNumber(0);
      float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
      float randomConePt = 0;

      int nmode = cfgnMods->at(i);
      int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      if (nmode == 2) {
        double rcPhiPsi2;
        float evtPl2 = helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode);
        rcPhiPsi2 = randomConePhi - evtPl2;

        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
        // removing the leading jet from the random cone
        if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
          float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
          float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

          bool jetWasInCone = false;
          while (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR) {
            jetWasInCone = true;
            randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
            randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
            dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
            dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
          }
          if (jetWasInCone) {
            randomConePt = 0.0;
            for (auto const& track : tracks) {
              if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
                float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-o2::constants::math::PI));
                float dEta = track.eta() - randomConeEta;
                if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
                  randomConePt += track.pt();
                }
              }
            }
          }
        }
        registry.fill(HIST("h3_centrality_deltapT_RandomCornPhi_RCprocess_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), rcPhiPsi2, 1.0);
      } else if (nmode == 3) {
        continue;
      }
    }
  }
  PROCESS_SWITCH(JetChargedV2, processRandomConeDataV2, "QA for random cone estimation of background fluctuations in data", true);

  void processTracksQA(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                       soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h2_centrality_track_pt"), collision.centrality(), track.pt());
      registry.fill(HIST("h2_centrality_track_eta"), collision.centrality(), track.eta());
      registry.fill(HIST("h2_centrality_track_phi"), collision.centrality(), track.phi());
    }
  }
  PROCESS_SWITCH(JetChargedV2, processTracksQA, "QA for charged tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetChargedV2>(cfgc)};
}
