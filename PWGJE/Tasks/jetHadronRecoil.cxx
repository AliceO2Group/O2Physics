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
/// \file jetHadronRecoil.cxx
/// \brief Task for analysing hadron triggered events.
/// \author Daniel Jones <djones22@liverpool.ac.uk>

#include <vector>
#include <string>

#include "TRandom3.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

#include "PWGJE/DataModel/JetSubstructure.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetHadronRecoil {

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> maxLeadingTrackPt{"maxLeadingTrackPt", 1000.0, "maximum acceptance for leading track in jets"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> ptTTrefMin{"ptTTrefMin", 5, "reference minimum trigger track pt"};
  Configurable<float> ptTTrefMax{"ptTTrefMax", 7, "reference maximum trigger track pt"};
  Configurable<float> ptTTsigMin{"ptTTsigMin", 20, "signal minimum trigger track pt"};
  Configurable<float> ptTTsigMax{"ptTTsigMax", 50, "signal maximum trigger track pt"};
  Configurable<float> fracSig{"fracSig", 0.9, "fraction of events to use for signal"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatTrackMaxMCD{"pTHatTrackMaxMCD", 999.0, "maximum fraction of hard scattering for track acceptance in detector MC"};
  Configurable<float> pTHatTrackMaxMCP{"pTHatTrackMaxMCP", 999.0, "maximum fraction of hard scattering for track acceptance in particle MC"};
  Configurable<float> rhoReferenceShift{"rhoReferenceShift", 0.0, "shift in rho calculated in reference events for consistency with signal events"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection applied at the jet finder level, here rejection is applied for collision and track process functions"};
  Configurable<int> wtaMethod{"wtaMethod", 1, "method for WTA axis definition: 0 = matching closest WTA jet (incorrect), 1 = recluster original jet"};

  Preslice<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> partJetsPerCollision = aod::jet::mcCollisionId;

  TRandom3* rand = new TRandom3(0);

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter particleCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  std::vector<double> ptBinningPart = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                       15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
                                       65.0, 70.0, 75.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0,
                                       140.0, 150.0, 160.0, 180.0, 200.0};
  std::vector<double> ptBinningDet = {-100.0, -70.0, -60.0, -50.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0,
                                      0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                      15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
                                      65.0, 70.0, 75.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0,
                                      140.0, 150.0, 160.0, 180.0, 200.0};
  std::vector<double> dRBinning = {0.0, 1.0e-9, 0.003, 0.006, 0.009, 0.012, 0.015, 0.018, 0.021, 0.024,
                                   0.027, 0.03, 0.033, 0.036, 0.039, 0.042, 0.045, 0.048, 0.051, 0.054,
                                   0.057, 0.06, 0.063, 0.066, 0.069, 0.072, 0.075, 0.078, 0.081, 0.084,
                                   0.087, 0.09, 0.093, 0.096, 0.099, 0.102, 0.105, 0.108, 0.111, 0.114,
                                   0.117, 0.12, 0.123, 0.126, 0.129, 0.132, 0.135, 0.138, 0.141, 0.144,
                                   0.147, 0.15, 0.153, 0.156, 0.159, 0.162, 0.165, 0.168, 0.171, 0.174,
                                   0.177, 0.18, 0.183, 0.186, 0.189, 0.192, 0.195, 0.198, 0.201, 0.204,
                                   0.207, 0.21, 0.213, 0.216, 0.219, 0.222, 0.225, 0.228, 0.231, 0.234,
                                   0.237, 0.24, 0.27, 0.30, 0.33, 0.36, 0.39, 0.42, 0.45, 0.48, 0.51, 0.54,
                                   0.57, 0.60};

  AxisSpec dRAxis = {dRBinning, "#Delta R"};
  AxisSpec ptAxisDet = {ptBinningDet, "#it{p}_{T,det} (GeV/c)"};
  AxisSpec ptAxisPart = {ptBinningPart, "#it{p}_{T,part} (GeV/c)"};
  AxisSpec phiAxisDet = {100, 0.0, o2::constants::math::TwoPI, "#phi_{det}"};
  AxisSpec phiAxisPart = {100, 0.0, o2::constants::math::TwoPI, "#phi_{part}"};
  AxisSpec dRAxisDet = {dRBinning, "#Delta R_{det}"};
  AxisSpec dRAxisPart = {dRBinning, "#Delta R_{part}"};

  HistogramRegistry registry{"registry",
                             {{"hNtrig", "number of triggers;trigger type;entries", {HistType::kTH1F, {{2, 0, 2}}}},
                              {"hSignalTriggersPtHard", "Signal triggers vs PtHard", {HistType::kTH1F, {{20, 0, 5}}}},
                              {"hReferenceTriggersPtHard", "Reference triggers vs PtHard", {HistType::kTH1F, {{20, 0, 5}}}},
                              {"hZvtxSelected", "Z vertex position;Z_{vtx};entries", {HistType::kTH1F, {{80, -20, 20}}}},
                              {"hPtTrack", "Track p_{T};p_{T};entries", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"hEtaTrack", "Track #eta;#eta;entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hPhiTrack", "Track #phi;#phi;entries", {HistType::kTH1F, {{100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hTrack3D", "3D tracks histogram;p_{T};#eta;#phi", {HistType::kTH3F, {{200, 0, 200}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hPtTrackPtHard", "Track p_{T} vs #hat{p};p_{T};#frac{p_{T}}{#hat{p}}", {HistType::kTH2F, {{200, 0, 200}, {20, 0, 5}}}},
                              {"hConstituents3D", "3D constituents histogram;p_{T};#eta;#phi", {HistType::kTH3F, {{200, 0, 200}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hReferencePtDPhi", "jet p_{T} vs DPhi;#Delta#phi;p_{T,jet}", {HistType::kTH2F, {{100, 0, o2::constants::math::TwoPI}, {500, -100, 400}}}},
                              {"hReferencePtDPhiShifts", "rho shifts;#Delta#phi;p_{T,jet};shifts", {HistType::kTH3F, {{100, 0, o2::constants::math::TwoPI}, {500, -100, 400}, {20, 0.0, 2.0}}}},
                              {"hSignalPtDPhi", "jet p_{T} vs DPhi;#Delta#phi;p_{T,jet}", {HistType::kTH2F, {{100, 0, o2::constants::math::TwoPI}, {500, -100, 400}}}},
                              {"hReferencePt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{500, -100, 400}}}},
                              {"hSignalPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{500, -100, 400}}}},
                              {"hSignalTriggers", "trigger p_{T};p_{T,trig};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hSignalPtHard", "jet p_{T} vs #hat{p};p_{T,jet};#frac{p_{T,trig}}{#hat{p}}", {HistType::kTH2F, {{500, -100, 400}, {20, 0, 5}}}},
                              {"hReferenceTriggers", "trigger p_{T};p_{T,trig};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hReferencePtHard", "jet p_{T} vs #hat{p};p_{T,jet};#frac{p_{T,trig}}{#hat{p}}", {HistType::kTH2F, {{500, -100, 400}, {20, 0, 5}}}},
                              {"hSigEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hRefEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{500, -100, 400}}}},
                              {"hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hJet3D", "3D jet distribution;p_{T};#eta;#phi", {HistType::kTH3F, {{500, -100, 400}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hPtPart", "Particle p_{T};p_{T};entries", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"hEtaPart", "Particle #eta;#eta;entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hPhiPart", "Particle #phi;#phi;entries", {HistType::kTH1F, {{100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hPart3D", "3D tracks histogram;p_{T};#eta;#phi", {HistType::kTH3F, {{200, 0, 200}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hPtPartPtHard", "Track p_{T} vs #hat{p};p_{T};#frac{p_{T}}{#hat{p}}", {HistType::kTH2F, {{200, 0, 200}, {20, 0, 5}}}},
                              {"hDeltaR", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRPart", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRpT", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{500, -100, 400}, dRAxis}}},
                              {"hDeltaRpTPart", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{400, 0, 400}, dRAxis}}},
                              {"hRhoSignal", "Signal Rho bkg;#rho;entries", {HistType::kTH1F, {{220, 0, 220}}}},
                              {"hRhoReference", "Reference Rho bkg;#rho;entries", {HistType::kTH1F, {{220, 0, 220}}}},
                              {"hRhoReferenceShift", "Testing reference shifts;#rho;shift", {HistType::kTH2F, {{220, 0, 220}, {20, 0.0, 2.0}}}},
                              {"hDeltaRSignal", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRSignalPart", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRpTSignal", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{500, -100, 400}, dRAxis}}},
                              {"hDeltaRpTSignalPart", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{400, 0, 400}, dRAxis}}},
                              {"hDeltaRpTDPhiSignal", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{500, -100, 400}, {100, 0, o2::constants::math::TwoPI}, dRAxis}}},
                              {"hDeltaRpTDPhiSignalPart", "Particle jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{400, 0, 400}, {100, 0, o2::constants::math::TwoPI}, dRAxis}}},
                              {"hDeltaRReference", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRPartReference", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {dRAxis}}},
                              {"hDeltaRpTReference", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{500, -100, 400}, dRAxis}}},
                              {"hDeltaRpTPartReference", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{400, 0, 400}, dRAxis}}},
                              {"hDeltaRpTDPhiReference", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{500, -100, 400}, {100, 0, o2::constants::math::TwoPI}, dRAxis}}},
                              {"hDeltaRpTDPhiReferenceShifts", "testing shifts;p_{T,jet};#Delta#phi;#DeltaR;shifts", {HistType::kTHnSparseD, {{500, -100, 400}, {100, 0, o2::constants::math::TwoPI}, dRAxis, {20, 0.0, 2.0}}}},
                              {"hDeltaRpTDPhiReferencePart", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{400, 0, 400}, {100, 0, o2::constants::math::TwoPI}, dRAxis}}},
                              {"hPtMatched", "p_{T} matching;p_{T,det};p_{T,part}", {HistType::kTH2F, {{500, -100, 400}, {400, 0, 400}}}},
                              {"hPhiMatched", "#phi matching;#phi_{det};#phi_{part}", {HistType::kTH2F, {{100, 0.0, o2::constants::math::TwoPI}, {100, 0.0, o2::constants::math::TwoPI}}}},
                              {"hDeltaRMatched", "#DeltaR matching;#DeltaR_{det};#DeltaR_{part}", {HistType::kTH2F, {dRAxisDet, dRAxisPart}}},
                              {"hPtMatched1d", "p_{T} matching 1d;p_{T,part}", {HistType::kTH1F, {{400, 0, 400}}}},
                              {"hDeltaRMatched1d", "#DeltaR matching 1d;#DeltaR_{part}", {HistType::kTH1F, {dRAxisPart}}},
                              {"hPtResolution", "p_{T} resolution;p_{T,part};Relative Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -5.0, 5.0}}}},
                              {"hPhiResolution", "#phi resolution;#p{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -7.0, 7.0}}}},
                              {"hDeltaRResolution", "#DeltaR Resolution;p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -0.15, 0.15}}}},
                              {"hFullMatching", "Full 6D matching;p_{T,det};p_{T,part};#phi_{det};#phi_{part};#DeltaR_{det};#DeltaR_{part}", {HistType::kTHnSparseD, {ptAxisDet, ptAxisPart, phiAxisDet, phiAxisPart, dRAxisDet, dRAxisPart}}}}};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  Service<o2::framework::O2DatabasePDG> pdg;

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

    Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);
    Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
    Filter particleCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);
    Filter eventTrackLevelCuts = nabs(aod::jcollision::posZ) < vertexZCut;

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = 2 * jetR;
    jetReclusterer.recombScheme = fastjet::WTA_pt_scheme;
  }

  template <typename T, typename U, typename W>
  void fillHistograms(T const& jets, W const& jetsWTA, U const& tracks, float weight = 1.0, float rho = 0.0)
  {
    bool isSigCol;
    std::vector<double> phiTTAr;
    std::vector<double> ptTTAr;
    double phiTT = 0;
    double ptTT = 0;
    int trigNumber = 0;
    int nTT = 0;
    double leadingPT = 0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    float rhoReference = rho + rhoReferenceShift;

    float dice = rand->Rndm();
    if (dice < fracSig)
      isSigCol = true;
    else
      isSigCol = false;

    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      if (track.pt() > pTHatTrackMaxMCD * pTHat) {
        return;
      }
      if (isSigCol && track.pt() < ptTTsigMax && track.pt() > ptTTsigMin) {
        phiTTAr.push_back(track.phi());
        ptTTAr.push_back(track.pt());
        registry.fill(HIST("hSignalTriggers"), track.pt(), weight);
        nTT++;
      }
      if (!isSigCol && track.pt() < ptTTrefMax && track.pt() > ptTTrefMin) {
        phiTTAr.push_back(track.phi());
        ptTTAr.push_back(track.pt());
        registry.fill(HIST("hReferenceTriggers"), track.pt(), weight);
        nTT++;
      }
      registry.fill(HIST("hPtTrack"), track.pt(), weight);
      registry.fill(HIST("hEtaTrack"), track.eta(), weight);
      registry.fill(HIST("hPhiTrack"), track.phi(), weight);
      registry.fill(HIST("hTrack3D"), track.pt(), track.eta(), track.phi(), weight);
      registry.fill(HIST("hPtTrackPtHard"), track.pt(), track.pt() / pTHat, weight);
    }

    if (nTT > 0) {
      trigNumber = rand->Integer(nTT);
      phiTT = phiTTAr[trigNumber];
      ptTT = ptTTAr[trigNumber];
      if (isSigCol) {
        registry.fill(HIST("hNtrig"), 1.5, weight);
        registry.fill(HIST("hSigEventTriggers"), nTT, weight);
        registry.fill(HIST("hRhoSignal"), rho, weight);
        registry.fill(HIST("hSignalTriggersPtHard"), ptTT / pTHat, weight);
      }
      if (!isSigCol) {
        registry.fill(HIST("hNtrig"), 0.5, weight);
        registry.fill(HIST("hRefEventTriggers"), nTT, weight);
        registry.fill(HIST("hRhoReference"), rhoReference, weight);
        for (double shift = 0.0; shift <= 2.0; shift += 0.1) {
          registry.fill(HIST("hRhoReferenceShift"), rho + shift, shift, weight);
        }
        registry.fill(HIST("hReferenceTriggersPtHard"), ptTT / pTHat, weight);
      }
    }

    for (const auto& jet : jets) {
      if (jet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
      for (const auto& constituent : jet.template tracks_as<U>()) {
        if (constituent.pt() > leadingPT) {
          leadingPT = constituent.pt();
        }
        registry.fill(HIST("hConstituents3D"), constituent.pt(), constituent.eta(), constituent.phi());
      }
      if (leadingPT > maxLeadingTrackPt) {
        continue;
      }
      registry.fill(HIST("hJetPt"), jet.pt() - (rho * jet.area()), weight);
      registry.fill(HIST("hJetEta"), jet.eta(), weight);
      registry.fill(HIST("hJetPhi"), jet.phi(), weight);
      registry.fill(HIST("hJet3D"), jet.pt() - (rho * jet.area()), jet.eta(), jet.phi(), weight);

      double dR = getWTAaxisDifference(jet, jetsWTA, tracks);

      registry.fill(HIST("hDeltaR"), dR, weight);
      registry.fill(HIST("hDeltaRpT"), jet.pt() - (rho * jet.area()), dR, weight);
      // try with fjcontrib

      if (nTT > 0) {
        float dphi = RecoDecay::constrainAngle(jet.phi() - phiTT);
        if (isSigCol) {
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hDeltaRpTSignal"), jet.pt() - (rho * jet.area()), dR, weight);
            registry.fill(HIST("hDeltaRSignal"), dR, weight);
          }
          registry.fill(HIST("hDeltaRpTDPhiSignal"), jet.pt() - (rho * jet.area()), dphi, dR, weight);
          registry.fill(HIST("hSignalPtDPhi"), dphi, jet.pt() - (rho * jet.area()), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hSignalPt"), jet.pt() - (rho * jet.area()), weight);
            registry.fill(HIST("hSignalPtHard"), jet.pt() - (rho * jet.area()), ptTT / pTHat, weight);
          }
        }
        if (!isSigCol) {
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hDeltaRpTReference"), jet.pt() - (rhoReference * jet.area()), dR, weight);
            registry.fill(HIST("hDeltaRReference"), dR, weight);
          }
          registry.fill(HIST("hDeltaRpTDPhiReference"), jet.pt() - (rhoReference * jet.area()), dphi, dR, weight);
          for (double shift = 0.0; shift <= 2.0; shift += 0.1) {
            registry.fill(HIST("hDeltaRpTDPhiReferenceShifts"), jet.pt() - ((rho + shift) * jet.area()), dphi, dR, shift, weight);
          }
          registry.fill(HIST("hReferencePtDPhi"), dphi, jet.pt() - (rhoReference * jet.area()), weight);
          for (double shift = 0.0; shift <= 2.0; shift += 0.1) {
            registry.fill(HIST("hReferencePtDPhiShifts"), dphi, jet.pt() - ((rho + shift) * jet.area()), shift, weight);
          }
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hReferencePt"), jet.pt() - (rhoReference * jet.area()), weight);
            registry.fill(HIST("hReferencePtHard"), jet.pt() - (rhoReference * jet.area()), ptTT / pTHat, weight);
          }
        }
      }
    }
  }

  template <typename T, typename W, typename U>
  void fillMCPHistograms(T const& jets, W const& jetsWTA, U const& particles, float weight = 1.0)
  {
    bool isSigCol;
    std::vector<double> phiTTAr;
    std::vector<double> ptTTAr;
    double phiTT = 0;
    double ptTT = 0;
    int trigNumber = 0;
    int nTT = 0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    float dice = rand->Rndm();
    if (dice < fracSig)
      isSigCol = true;
    else
      isSigCol = false;

    for (const auto& particle : particles) {
      if (particle.pt() > pTHatTrackMaxMCD * pTHat) {
        return;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle) {
        continue;
      }
      if ((pdgParticle->Charge() == 0.0) || (!particle.isPhysicalPrimary())) {
        continue;
      }
      if (isSigCol && particle.pt() < ptTTsigMax && particle.pt() > ptTTsigMin) {
        phiTTAr.push_back(particle.phi());
        ptTTAr.push_back(particle.pt());
        nTT++;
      }
      if (!isSigCol && particle.pt() < ptTTrefMax && particle.pt() > ptTTrefMin) {
        phiTTAr.push_back(particle.phi());
        ptTTAr.push_back(particle.pt());
        nTT++;
      }
      registry.fill(HIST("hPtPart"), particle.pt(), weight);
      registry.fill(HIST("hEtaPart"), particle.eta(), weight);
      registry.fill(HIST("hPhiPart"), particle.phi(), weight);
      registry.fill(HIST("hPart3D"), particle.pt(), particle.eta(), particle.phi(), weight);
      registry.fill(HIST("hPtPartPtHard"), particle.pt(), particle.pt() / pTHat, weight);
    }

    if (nTT > 0) {
      trigNumber = rand->Integer(nTT);
      phiTT = phiTTAr[trigNumber];
      ptTT = ptTTAr[trigNumber];
      if (isSigCol) {
        registry.fill(HIST("hNtrig"), 1.5, weight);
        registry.fill(HIST("hSigEventTriggers"), nTT, weight);
        registry.fill(HIST("hSignalTriggersPtHard"), ptTT / pTHat, weight);
      }
      if (!isSigCol) {
        registry.fill(HIST("hNtrig"), 0.5, weight);
        registry.fill(HIST("hRefEventTriggers"), nTT, weight);
        registry.fill(HIST("hReferenceTriggersPtHard"), ptTT / pTHat, weight);
      }
    }

    for (const auto& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat) {
        return;
      }
      for (const auto& constituent : jet.template tracks_as<U>()) {
        registry.fill(HIST("hConstituents3D"), constituent.pt(), constituent.eta(), constituent.phi());
      }
      registry.fill(HIST("hJetPt"), jet.pt(), weight);
      registry.fill(HIST("hJetEta"), jet.eta(), weight);
      registry.fill(HIST("hJetPhi"), jet.phi(), weight);
      registry.fill(HIST("hJet3D"), jet.pt(), jet.eta(), jet.phi(), weight);

      double dR = getWTAaxisDifference(jet, jetsWTA, particles);

      registry.fill(HIST("hDeltaRPart"), dR, weight);
      registry.fill(HIST("hDeltaRpTPart"), jet.pt(), dR, weight);
      if (nTT > 0) {
        float dphi = RecoDecay::constrainAngle(jet.phi() - phiTT);
        if (isSigCol) {
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hDeltaRpTSignalPart"), jet.pt(), dR, weight);
            registry.fill(HIST("hDeltaRSignalPart"), dR, weight);
          }
          registry.fill(HIST("hDeltaRpTDPhiSignalPart"), jet.pt(), dphi, dR, weight);
          registry.fill(HIST("hSignalPtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hSignalPt"), jet.pt(), weight);
            registry.fill(HIST("hSignalPtHard"), jet.pt(), ptTT / pTHat, weight);
          }
        }
        if (!isSigCol) {
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hDeltaRpTPartReference"), jet.pt(), dR, weight);
            registry.fill(HIST("hDeltaRPartReference"), dR, weight);
          }
          registry.fill(HIST("hDeltaRpTDPhiReferencePart"), jet.pt(), dphi, dR, weight);
          registry.fill(HIST("hReferencePtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hReferencePt"), jet.pt(), weight);
            registry.fill(HIST("hReferencePtHard"), jet.pt(), ptTT / pTHat, weight);
          }
        }
      }
    }
  }

  template <typename T, typename V, typename W, typename U, typename X, typename Y>
  void fillMatchedHistograms(T const& jetsBase, V const& mcdjetsWTA, W const& mcpjetsWTA, U const&, X const& tracks, Y const& particles, float weight = 1.0, float rho = 0.0)
  {
    for (const auto& jetBase : jetsBase) {
      double dR = 0;
      double dRp = 0;

      float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
      if (jetBase.pt() > pTHatMaxMCD * pTHat) {
        return;
      }

      dR = getWTAaxisDifference(jetBase, mcdjetsWTA, tracks, true);

      if (jetBase.has_matchedJetGeo()) {
        for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
          if (jetTag.pt() > pTHatMaxMCP * pTHat) {
            return;
          }

          dRp = getWTAaxisDifference(jetTag, mcpjetsWTA, particles, true);

          registry.fill(HIST("hPtMatched"), jetBase.pt() - (rho * jetBase.area()), jetTag.pt(), weight);
          registry.fill(HIST("hPhiMatched"), jetBase.phi(), jetTag.phi(), weight);
          registry.fill(HIST("hPtResolution"), jetTag.pt(), (jetTag.pt() - (jetBase.pt() - (rho * jetBase.area()))) / jetTag.pt(), weight);
          registry.fill(HIST("hPhiResolution"), jetTag.pt(), jetTag.phi() - jetBase.phi(), weight);
          registry.fill(HIST("hDeltaRMatched"), dR, dRp, weight);
          registry.fill(HIST("hDeltaRResolution"), jetTag.pt(), dRp - dR, weight);
          registry.fill(HIST("hFullMatching"), jetBase.pt() - (rho * jetBase.area()), jetTag.pt(), jetBase.phi(), jetTag.phi(), dR, dRp, weight);
          registry.fill(HIST("hPtMatched1d"), jetTag.pt(), weight);
          registry.fill(HIST("hDeltaRMatched1d"), dRp, weight);
        }
      }
    }
  }

  template <typename T, typename V, typename W, typename U, typename X, typename Y>
  void fillRecoilJetMatchedHistograms(T const& jetsBase, V const& mcdjetsWTA, W const& mcpjetsWTA, U const&, X const& tracks, Y const& particles, float weight = 1.0, float rho = 0.0)
  {
    std::vector<double> phiTTAr;
    std::vector<double> phiTTArPart;
    double phiTT = 0;
    double phiTTPart = 0;
    int trigNumber = 0;
    int nTT = 0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      if (track.pt() > pTHatTrackMaxMCD * pTHat) {
        return;
      }
      if (track.pt() < ptTTsigMax && track.pt() > ptTTsigMin) {
        phiTTAr.push_back(track.phi());
        nTT++;
        auto particle = track.template mcParticle_as<Y>();
        phiTTArPart.push_back(particle.phi());
      }
    }

    if (nTT > 0) {
      trigNumber = rand->Integer(nTT);
      phiTT = phiTTAr[trigNumber];
      phiTTPart = phiTTAr[trigNumber];
    } else {
      return;
    }

    for (const auto& jetBase : jetsBase) {
      double dR = 0;
      double dRp = 0;

      if (jetBase.pt() > pTHatMaxMCD * pTHat) {
        return;
      }

      float dphi = RecoDecay::constrainAngle(jetBase.phi() - phiTT);
      dR = getWTAaxisDifference(jetBase, mcdjetsWTA, tracks, true);

      if (jetBase.has_matchedJetGeo()) {
        for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
          if (jetTag.pt() > pTHatMaxMCP * pTHat) {
            return;
          }

          float dphip = RecoDecay::constrainAngle(jetTag.phi() - phiTTPart);
          dRp = getWTAaxisDifference(jetTag, mcpjetsWTA, particles, true);

          registry.fill(HIST("hPtMatched"), jetBase.pt() - (rho * jetBase.area()), jetTag.pt(), weight);
          registry.fill(HIST("hPhiMatched"), dphi, dphip, weight);
          registry.fill(HIST("hPtResolution"), jetTag.pt(), (jetTag.pt() - (jetBase.pt() - (rho * jetBase.area()))) / jetTag.pt(), weight);
          registry.fill(HIST("hPhiResolution"), dphip, dphip - dphi, weight);
          registry.fill(HIST("hDeltaRMatched"), dR, dRp, weight);
          registry.fill(HIST("hDeltaRResolution"), jetTag.pt(), dRp - dR, weight);
          registry.fill(HIST("hFullMatching"), jetBase.pt() - (rho * jetBase.area()), jetTag.pt(), dphi, dphip, dR, dRp, weight);
          registry.fill(HIST("hPtMatched1d"), jetTag.pt(), weight);
          registry.fill(HIST("hDeltaRMatched1d"), dRp, weight);
        }
      }
    }
  }

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToCharged1Jets>> const& jets,
                   soa::Filtered<soa::Join<aod::Charged1Jets, aod::Charged1JetConstituents, aod::Charged1JetsMatchedToChargedJets>> const& jetsWTA,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks);
  }
  PROCESS_SWITCH(JetHadronRecoil, processData, "process data", true);

  void processDataWithRhoSubtraction(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                     soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToCharged1Jets>> const& jets,
                                     soa::Filtered<soa::Join<aod::Charged1Jets, aod::Charged1JetConstituents, aod::Charged1JetsMatchedToChargedJets>> const& jetsWTA,
                                     soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks, 1.0, collision.rho());
  }
  PROCESS_SWITCH(JetHadronRecoil, processDataWithRhoSubtraction, "process data with rho subtraction", false);

  void processMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                  soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                  soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                  soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks);
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCD, "process MC detector level", false);

  void processMCDWithRhoSubtraction(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                    soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                                    soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                                    soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks, 1.0, collision.rho());
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCDWithRhoSubtraction, "process MC detector level with rho subtraction", false);

  void processMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                          aod::JMcCollisions const&,
                          soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                          soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                          soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ(), collision.mcCollision().weight());
    fillHistograms(jets, jetsWTA, tracks, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCDWeighted, "process MC detector level with event weights", false);

  void processMCDWeightedWithRhoSubtraction(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::BkgChargedRhos>>::iterator const& collision,
                                            aod::JMcCollisions const&,
                                            soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                                            soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                                            soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ(), collision.mcCollision().weight());
    fillHistograms(jets, jetsWTA, tracks, collision.mcCollision().weight(), collision.rho());
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCDWeightedWithRhoSubtraction, "process MC detector level with event weights and rho subtraction", false);

  void processMCP(aod::JetMcCollision const& collision,
                  soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets>> const& jets,
                  soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents, aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets>> const& jetsWTA,
                  soa::Filtered<aod::JetParticles> const& particles)
  {
    if (std::abs(collision.posZ()) > vertexZCut) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillMCPHistograms(jets, jetsWTA, particles);
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCP, "process MC particle level", false);

  void processMCPWeighted(aod::JetMcCollision const& collision,
                          soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets>> const& jets,
                          soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents, aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets>> const& jetsWTA,
                          soa::Filtered<aod::JetParticles> const& particles)
  {
    if (std::abs(collision.posZ()) > vertexZCut) {
      return;
    }
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ(), collision.weight());
    fillMCPHistograms(jets, jetsWTA, particles, collision.weight());
  }
  PROCESS_SWITCH(JetHadronRecoil, processMCPWeighted, "process MC particle level with event weights", false);

  void processJetsMCPMCDMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                aod::JetTracks const& tracks,
                                aod::JetParticles const& particles,
                                aod::JetMcCollisions const&,
                                soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles);
  }
  PROCESS_SWITCH(JetHadronRecoil, processJetsMCPMCDMatched, "process MC matched (inc jets)", false);

  void processJetsMCPMCDMatchedWithRhoSubtraction(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
                                                  soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                                  soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                                  soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                                  aod::JetTracks const& tracks,
                                                  aod::JetParticles const& particles,
                                                  aod::JetMcCollisions const&,
                                                  soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles);
  }
  PROCESS_SWITCH(JetHadronRecoil, processJetsMCPMCDMatchedWithRhoSubtraction, "process MC matched (inc jets) with rho subtraction", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                        soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                        soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                        soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                        aod::JetTracks const& tracks,
                                        aod::JetParticles const& particles,
                                        aod::JetMcCollisions const&,
                                        soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetHadronRecoil, processJetsMCPMCDMatchedWeighted, "process MC matched with event weights (inc jets)", false);

  void processJetsMCPMCDMatchedWeightedWithRhoSubtraction(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
                                                          soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                                          soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                                          soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                                          aod::JetTracks const& tracks,
                                                          aod::JetParticles const& particles,
                                                          aod::JetMcCollisions const&,
                                                          soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles, collision.mcCollision().weight(), collision.rho());
  }
  PROCESS_SWITCH(JetHadronRecoil, processJetsMCPMCDMatchedWeightedWithRhoSubtraction, "process MC matched with event weights (inc jets) and rho subtraction", false);

  void processRecoilJetsMCPMCDMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                      soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                      soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                      soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                      soa::Filtered<aod::JetTracksMCD> const& tracks,
                                      soa::Filtered<aod::JetParticles> const& particles,
                                      aod::JetMcCollisions const&,
                                      soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillRecoilJetMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles);
  }
  PROCESS_SWITCH(JetHadronRecoil, processRecoilJetsMCPMCDMatched, "process MC matched (recoil jets)", false);

  void processRecoilJetsMCPMCDMatchedWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                              soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                              soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                              soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                              soa::Filtered<aod::JetTracksMCD> const& tracks,
                                              soa::Filtered<aod::JetParticles> const& particles,
                                              aod::JetMcCollisions const&,
                                              soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillRecoilJetMatchedHistograms(mcdjets, mcdjetsWTA, mcpjetsWTACut, mcpjets, tracks, particles, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetHadronRecoil, processRecoilJetsMCPMCDMatchedWeighted, "process MC matched with event weights (recoil jets)", false);

  template <typename T, typename U, typename X>
  double getWTAaxisDifference(T const& jet, U const& jetsWTA, X const& /*tracks or particles*/, bool isMatched = false)
  {
    double deltaPhi = -1;
    double deltaEta = -1;
    double deltaY = -1;
    double dR = -1;
    if (wtaMethod == 0) {
      // get WTA matched jet - should just be one jet matched geometrically
      if (isMatched) {
        // response - requires alternative method
        for (const auto& jetWTA : jetsWTA) {
          double djetp = RecoDecay::sqrtSumOfSquares(RecoDecay::constrainAngle(jet.phi() - jetWTA.phi(), -o2::constants::math::PI), jet.eta() - jetWTA.eta());
          if (djetp < 0.6 * jetR) {
            dR = djetp;
            break;
          }
        }
      } else {
        // MCP or MCD
        for (const auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<U>>()) {
          deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
          deltaEta = jetWTA.eta() - jet.eta();
          dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
        }
      }
    } else if (wtaMethod == 1) {
      // recluster jet
      jetConstituents.clear();
      for (auto& jetConstituent : jet.template tracks_as<X>()) {
        fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
      }
      jetReclustered.clear();
      fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
      jetReclustered = sorted_by_pt(jetReclustered);

      deltaPhi = RecoDecay::constrainAngle(jet.phi() - jetReclustered[0].phi(), -o2::constants::math::PI);
      deltaY = jet.y() - jetReclustered[0].rap();
      dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaY);
      LOG(debug) << "orig. jet n const = " << jet.tracksIds().size() << " pt = " << jet.pt() << " eta = " << jet.eta() << " phi = " << jet.phi();
      LOG(debug) << "recl. jet n const = " << clusterSeq.constituents(jetReclustered[0]).size() << " pt = " << jetReclustered[0].pt() << " eta = " << jetReclustered[0].eta() << " phi = " << jetReclustered[0].phi();
      LOG(debug) << "distance = " << dR;
    }
    return dR;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetHadronRecoil>(cfgc)}; }
