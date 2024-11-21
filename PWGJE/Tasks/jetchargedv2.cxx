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

// jet v2 task
/// \author Yubiao Wang <yubiao.wang@cern.ch>
// C++/ROOT includes.

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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//=====================< evt pln >=====================//
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
//=====================< evt pln | end >=====================//

struct Jetchargedv2Task {
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
  Configurable<float> jetPtMin{"jetPtMin", 0.15, "minimum pT acceptance for jets"};
  Configurable<float> jetPtMax{"jetPtMax", 200.0, "maximum pT acceptance for jets"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.9, "minimum eta acceptance for jets"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.9, "maximum eta acceptance for jets"};
  Configurable<float> jetRadius{"jetRadius", 0.2, "jet resolution parameters"};

  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};

  //=====================< evt pln >=====================//
  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "event selection"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0M", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};

  EventPlaneHelper helperEP;
  int DetId;
  int RefAId;
  int RefBId;

  template <typename T>
  int GetDetId(const T& name)
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

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  int eventSelection = -1;
  int trackSelection = -1;
  double evtnum = 1; // evt sum for local rho test

  void init(o2::framework::InitContext&)
  {
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);
    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
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
    //< bkg sub | end >//

    AxisSpec axisPt = {40, 0.0, 4.0};
    AxisSpec axisEta = {32, -0.8, 0.8};
    AxisSpec axixCent = {20, 0, 100};
    AxisSpec axisChID = {220, 0, 220};

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    //< \sigma p_T at local rho test plot >
    registry.add("h_ptsum_collnum", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{40, 0.0, 40}}});
    registry.add("h_ptsum_sumpt", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., TMath::TwoPi()}}});
    registry.add("h2_phi_track_eta", "phi vs track eta; #varphi; #eta (GeV/#it{c})", {HistType::kTH2F, {{100, -1.0, 1.0}, {160, 0., TMath::TwoPi()}}});
    registry.add("h2_phi_track_pt", "phi vs track pT; #varphi; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {160, 0., TMath::TwoPi()}}});
    registry.add("h2_centrality_phi_w_pt", "centrality vs jet #varphi; #varphi_{jet}; entries", {HistType::kTH2F, {{100, 0.0, 100.0}, {160, 0., TMath::TwoPi()}}});
    registry.add("h2_evtnum_phi_w_pt", "eventNumber vs jet #varphi; #eventNumber; entries", {HistType::kTH2F, {{100000, 0.0, 100000}, {160, 0., TMath::TwoPi()}}});

    //< \sigma p_T at local rho test plot | end >
    registry.add("h2_centrality_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
    registry.add("h2_centrality_track_pt", "centrality vs track pT; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
    registry.add("h2_centrality_track_eta", "centrality vs track #eta; centrality; #eta_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
    registry.add("h2_centrality_track_phi", "centrality vs track #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});

    registry.add("h_recoil_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_recoil_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
    registry.add("h_recoil_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
    registry.add("h_recoil_jet_dphi", "hadron-jet #Delta#phi;#Delta#phi_{jet,trigger hadron};entries", {HistType::kTH1F, {{40, -2.0, 2.0}}});

    registry.add("leadJetPt", "leadJet Pt ", {HistType::kTH1F, {{200, 0., 200.0}}});
    registry.add("leadJetPhi", "leadJet constituent #phi ", {HistType::kTH1F, {{80, -1.0, 7.}}});
    registry.add("leadJetEta", "leadJet constituent #eta ", {HistType::kTH1F, {{100, -1.0, 1.0}}});

    //< RC test plots >//
    registry.add("h3_centrality_RCpt_RandomCornPhi_rhorandomcone", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{120, -10.0, 110.0}, {800, -400.0, 400.0}, {160, 0., TMath::TwoPi()}}});
    registry.add("h3_centrality_RCpt_RandomCornPhi_rhorandomconewithoutleadingjet", "centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho}; #Delta#varphi_{jet}", {HistType::kTH3F, {{120, -10.0, 110.0}, {800, -400.0, 400.0}, {160, 0., TMath::TwoPi()}}});
    //< bkg sub plot | end >//

    registry.add("h_jet_pt_in_plane_v2", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v2", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_in_plane_v3", "jet pT;#it{p}^{in-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
    registry.add("h_jet_pt_out_of_plane_v3", "jet pT;#it{p}^{out-of-plane}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});

    registry.add("h2_centrality_jet_pt_in_plane_v2", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v2", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_in_plane_v3", "centrality vs #it{p}^{in-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    registry.add("h2_centrality_jet_pt_out_of_plane_v3", "centrality vs #it{p}^{out-of-plane}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{120, -10.0, 110.0}, jetPtAxisRhoAreaSub}});
    //< bkg sub DeltaPhi plot | end >//

    //=====================< evt pln plot >=====================//
    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};

    AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("histCentFull", "Centrality distribution for valid events", HistType::kTH1F, {axisCent});
    for (std::size_t i = 0; i < cfgnMods->size(); i++) {
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

  Preslice<aod::ChargedJets> JetsPerJCollision = aod::jet::collisionId;
  Preslice<aod::JetTracks> tracksPerJCollision = o2::aod::jtrack::collisionId;

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > -98.0) {
      bool isMinleadingConstituent = false;
      for (auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }
      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }

  //=====================< q-vector & evtpln check >=====================//
  template <typename T>
  void fillHistosQvec(const T& vec, int nmode)
  {
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if (nmode == 2) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV2"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histQvecRectrV2"), vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], vec.cent());
        histosQA.fill(HIST("histQvecTwistV2"), vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], vec.cent());
        histosQA.fill(HIST("histQvecFinalV2"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
      }
    }
  }
  //=====================< q-vector & evtpln check | end >=====================//
  void fillLeadingJetQA(double leadingJetPt, double leadingJetPhi, double leadingJetEta)
  {
    registry.fill(HIST("leadJetPt"), leadingJetPt);
    registry.fill(HIST("leadJetPhi"), leadingJetPhi);
    registry.fill(HIST("leadJetEta"), leadingJetEta);
  } // end of fillLeadingJetQA template

  void processjetQA(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                    soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, aod::JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double leadingTrackpT = 0.0;
    double leadingTrackPhi = 0.0;
    for (auto const& track : tracks) {
      if (track.pt() > 6.0 && track.pt() < 10.0) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
          leadingTrackPhi = track.phi();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return;
    for (auto& jet : jets) {
      if (TMath::Abs(RecoDecay::constrainAngle(RecoDecay::constrainAngle(jet.phi(), -o2::constants::math::PIHalf) - RecoDecay::constrainAngle(leadingTrackPhi, -o2::constants::math::PIHalf), -o2::constants::math::PIHalf) > 0.6)) {
        registry.fill(HIST("h_recoil_jet_pt"), jet.pt());
        registry.fill(HIST("h_recoil_jet_eta"), jet.eta());
        registry.fill(HIST("h_recoil_jet_phi"), jet.phi());
        registry.fill(HIST("h_recoil_jet_dphi"), jet.phi() - leadingTrackPhi);
      }
    }
  }
  PROCESS_SWITCH(Jetchargedv2Task, processjetQA, "jet rho v2 jet QA", true);

  void processSigmaPt(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>> const& collisions,
                      soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                      aod::JetTracks const& tracks)
  {
    // double collnum = 1;
    for (const auto& collision : collisions) {
      double leadingJetPt = -1;
      double leadingJetPhi = -1;
      double leadingJetEta = -1;
      for (auto& jet : jets) {
        if (jet.pt() > leadingJetPt) {
          leadingJetPt = jet.pt();
          leadingJetEta = jet.eta();
          leadingJetPhi = jet.phi();
        }
      }
      fillLeadingJetQA(leadingJetPt, leadingJetPhi, leadingJetEta);

      if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
        return;
      }

      //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
      for (std::size_t i = 0; i < cfgnMods->size(); i++) {
        int nmode = cfgnMods->at(i);
        int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
        if (nmode == 2) {
          if (collision.qvecAmp()[DetId] > 1e-8) {
            histosQA.fill(HIST("histQvecUncorV2"), collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], collision.cent());
            histosQA.fill(HIST("histQvecRectrV2"), collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], collision.cent());
            histosQA.fill(HIST("histQvecTwistV2"), collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], collision.cent());
            histosQA.fill(HIST("histQvecFinalV2"), collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], collision.cent());

            histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode), collision.cent());
            histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], nmode), collision.cent());
            histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], nmode), collision.cent());
            histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), collision.cent());
          }
        } else if (nmode == 3) {
          histosQA.fill(HIST("histQvecUncorV3"), collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], collision.cent());
          histosQA.fill(HIST("histQvecRectrV3"), collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], collision.cent());
          histosQA.fill(HIST("histQvecTwistV3"), collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], collision.cent());
          histosQA.fill(HIST("histQvecFinalV3"), collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], collision.cent());

          histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlRectrV3"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlTwistV3"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], nmode), collision.cent());
          histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), collision.cent());
        }
        //< Psi_EP,2, JetPtCorr = Jet_pT-<rho>*A in-plane and out-of-plane >//
        auto collJets = jets.sliceBy(JetsPerJCollision, collision.globalIndex()); // select the jet in collisions
        if (nmode == 2) {
          Double_t phiMinusPsi2;
          if (collision.qvecAmp()[DetId] < 1e-8) {
            continue;
          }
          float evtPl2 = helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode);
          for (auto const& jet : collJets) {
            if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              continue;
            }
            if (!isAcceptedJet<aod::JetTracks>(jet)) {
              continue;
            }
            phiMinusPsi2 = jet.phi() - evtPl2;
            // Double_t jetPtCorr = 0.0;
            // jetPtCorr = jet.pt() - collision.rho() * jet.area();

            if ((phiMinusPsi2 < TMath::Pi() / 4) || (phiMinusPsi2 >= 7 * TMath::Pi() / 4) || (phiMinusPsi2 >= 3 * TMath::Pi() / 4 && phiMinusPsi2 < 5 * TMath::Pi() / 4)) {
              registry.fill(HIST("h_jet_pt_in_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
              registry.fill(HIST("h2_centrality_jet_pt_in_plane_v2"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
            } else {
              registry.fill(HIST("h_jet_pt_out_of_plane_v2"), jet.pt() - (collision.rho() * jet.area()), 1.0);
              registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v2"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
            }
          }
          //< JetPtCorr = Jet_pT-<rho>*A in-plane and out-of-plane | end >//
        } else if (nmode == 3) {
          Double_t phiMinusPsi3;
          float evtPl3 = helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode);
          for (auto const& jet : collJets) {
            if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              continue;
            }
            if (!isAcceptedJet<aod::JetTracks>(jet)) {
              continue;
            }
            phiMinusPsi3 = jet.phi() - evtPl3;
            // Double_t jetPtCorr = 0.0;
            // jetPtCorr = jet.pt() - collision.rho() * jet.area();

            if ((phiMinusPsi3 < TMath::Pi() / 4) || (phiMinusPsi3 >= 7 * TMath::Pi() / 4) || (phiMinusPsi3 >= 3 * TMath::Pi() / 4 && phiMinusPsi3 < 5 * TMath::Pi() / 4)) {
              registry.fill(HIST("h_jet_pt_in_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
              registry.fill(HIST("h2_centrality_jet_pt_in_plane_v3"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
            } else {
              registry.fill(HIST("h_jet_pt_out_of_plane_v3"), jet.pt() - (collision.rho() * jet.area()), 1.0);
              registry.fill(HIST("h2_centrality_jet_pt_out_of_plane_v3"), collision.centrality(), jet.pt() - (collision.rho() * jet.area()), 1.0);
            }
          }
        }
      }
      //=====================< evt pln | end >=====================//
      auto collTracks = tracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      if (jets.size() > 0) {
        for (auto const& track : collTracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection) && (fabs(track.eta() - leadingJetEta) > jetRadius) && track.pt() >= 0.2 && track.pt() <= 5.) {
            registry.fill(HIST("h2_phi_track_pt"), track.pt(), track.phi());
            registry.fill(HIST("h2_phi_track_eta"), track.eta(), track.phi());
            registry.fill(HIST("h_ptsum_sumpt"), track.phi(), track.pt());                                  // \sigma p_T distribution test
            registry.fill(HIST("h2_centrality_phi_w_pt"), collision.centrality(), track.phi(), track.pt()); // \sigma track.pt() distribution with centrality test
            registry.fill(HIST("h2_evtnum_phi_w_pt"), evtnum, track.phi(), track.pt());
          }
        }
      }
      registry.fill(HIST("h_ptsum_collnum"), 0.5);
      evtnum += 1;
    }
  }
  PROCESS_SWITCH(Jetchargedv2Task, processSigmaPt, "QA for charged tracks", true);

  void processRandomConeDataV2(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                               soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                               soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    for (std::size_t i = 0; i < cfgnMods->size(); i++) {
      TRandom3 randomNumber(0);
      float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      float randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
      float randomConePt = 0;

      int nmode = cfgnMods->at(i);
      int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

      Double_t RcPhiPsi2;
      float evtPl2 = helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode);
      RcPhiPsi2 = randomConePhi - evtPl2;

      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
          float dEta = track.eta() - randomConeEta;
          if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            randomConePt += track.pt();
          }
        }
      }
      registry.fill(HIST("h3_centrality_RCpt_RandomCornPhi_rhorandomcone"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho(), RcPhiPsi2, 1.0);
      // removing the leading jet from the random cone
      if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
        float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
        float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

        bool jetWasInCone = false;
        while (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR) {
          jetWasInCone = true;
          randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
          randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
          dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
          dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
        }
        if (jetWasInCone) {
          randomConePt = 0.0;
          for (auto const& track : tracks) {
            if (jetderiveddatautilities::selectTrack(track, trackSelection)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
              float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
              float dEta = track.eta() - randomConeEta;
              if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
                randomConePt += track.pt();
              }
            }
          }
        }
      }
      registry.fill(HIST("h3_centrality_RCpt_RandomCornPhi_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho(), RcPhiPsi2, 1.0);
    }
  }
  PROCESS_SWITCH(Jetchargedv2Task, processRandomConeDataV2, "QA for random cone estimation of background fluctuations in data", true);

  void processTracksQA(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                       soa::Filtered<aod::JetTracks> const& tracks)
  {
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 1.5);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h2_centrality_track_pt"), collision.centrality(), track.pt());
      registry.fill(HIST("h2_centrality_track_eta"), collision.centrality(), track.eta());
      registry.fill(HIST("h2_centrality_track_phi"), collision.centrality(), track.phi());
    }
  }
  PROCESS_SWITCH(Jetchargedv2Task, processTracksQA, "QA for charged tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Jetchargedv2Task>(cfgc, TaskName{"jet-charged-v2"})};
}
