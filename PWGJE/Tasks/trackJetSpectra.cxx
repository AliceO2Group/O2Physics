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

/// \file trackJetSpectra.cxx
/// \brief Analysis task for charged jet and track spectra studies
/// \author Janis Jan

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include "TRandom3.h"
#include <Math/Vector4D.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TVector2.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FilteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using FilteredEventMultiplicity = soa::Filtered<soa::Join<aod::JetCollisions, aod::ZDCMults>>::iterator;

using FilteredJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>>;
using FilteredJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredCollDetLevelGetWeight = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator;

using FilteredCollPartLevel = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>>::iterator;

using CollPartLevelEfficiencyWeighted = soa::Join<aod::JetMcCollisions, aod::JMcCollisionOutliers>::iterator;

using FilteredCollPartLevelMB = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator;

using CollPartLevelMBEfficiency = aod::JetMcCollisions::iterator;

using CollisionMCD = soa::SmallGroups<aod::JetCollisionsMCD>;

using FilteredParticles = soa::Filtered<aod::JetParticles>;

using FilteredTracksMatch = soa::Filtered<soa::Join<aod::JetTracks, aod::McTrackLabels>>;

using ParticlesFullTable = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

using TracksMCDetFullTable = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;

using FilteredTracks = soa::Filtered<aod::JetTracks>;

using FilteredMatchedJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using FilteredMatchedJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

struct TrackJetSpectra {

  // List of configurable parameters
  Configurable<std::string> evSel{"evSel", "sel8", "Choose event selection"};
  Configurable<std::string> trkSel{"trkSel", "globalTracks", "Set track selection"};
  Configurable<bool> setSumw2{"setSumw2", false, "Enable Sum2w() calculation for histograms"};
  Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};
  Configurable<float> maxJetConstituentPt{"maxJetConstituentPt", 100., "Remove jets with constituent above this pT cut"};

  Configurable<float> trkPtMin{"trkPtMin", 0.15, "Minimum pT of acceptanced tracks"};
  Configurable<float> trkPtMax{"trkPtMax", 100., "Maximum pT of acceptanced tracks"};

  Configurable<float> trkEtaCut{"trkEtaCut", 0.9, "Eta acceptance of TPC"};
  Configurable<float> jetR{"jetR", 0.4, "Jet cone radius"};

  Configurable<float> signalTriggerMin{"signalTriggerMin", 15.f, "Minimal signal trigger track pT"};

  Configurable<float> signalTriggerMax{"signalTriggerMax", 50.f, "Maximal signal trigger track pT"};

  Configurable<float> referenceTriggerMin{"referenceTriggerMin", 5.f, "Minimal reference trigger track pT"};

  Configurable<float> referenceTriggerMax{"referenceTriggerMax", 7.f, "Maximal reference trigger track pT"};

  Configurable<float> sigToRefFraction{"sigToRefFraction", 0.1f, "Signal-to-reference trigger fraction"};
  Configurable<float> dPhiCut{"dPhiCut", 0.6, ""};
  // List of configurable parameters for histograms
  Configurable<uint16_t> histJetPt{"histJetPt", 150, "Maximum value of jet pT shown in histograms"};

  Configurable<float> meanFT0A{"meanFT0A", -1., "Mean value of FT0A signal"};
  Configurable<float> meanFT0C{"meanFT0C", -1., "Mean value of FT0C signal"};

  Configurable<float> pionMass{"pionMass", o2::constants::physics::MassPiPlus, "Mass of a pion"};

  Configurable<bool> isMCJJProd{"isMCJJProd", false, "Flag to select MC production: MB (false) or JJ(true)"},
    skipMBGapEvents{"skipMBGapEvents", false,
                    "Flag to choose to reject min. bias gap events; jet-level rejection "
                    "applied at the jet finder level, here rejection is applied for "
                    "collision and track process functions"};

  //------------------------------------------------------------

  TRandom3* rand = new TRandom3(0);
  // Declare filter on collision Z vertex
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilterMC = nabs(aod::jmccollision::posZ) < vertexZCut;

  // Declare filters on accepted tracks and MC particles (settings for jet reco are provided in the jet finder wagon)
  Filter trackFilter = aod::jtrack::pt > trkPtMin&& aod::jtrack::pt < trkPtMax&& nabs(aod::jtrack::eta) < trkEtaCut;
  Filter partFilter = nabs(aod::jmcparticle::eta) < trkEtaCut;

  // Declare filter on jets
  Filter jetRadiusFilter = aod::jet::r == nround(jetR.node() * 100.);

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<FilteredMatchedJetsPartLevel> partJetsPerCollision = aod::jet::mcCollisionId;

  Preslice<aod::JetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;

  void init(InitContext const&)
  {
    // Initialize histogram axes
    AxisSpec pT{histJetPt + 20, -20.0, histJetPt * 1., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTRes{100, -5.0, 5.0, "#it{p}_{T} Resolution"};
    AxisSpec pTTrack{100, 0.0, 100.0, "#it{p}_{T,Track} (GeV/#it{c})"};
    AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
    AxisSpec dphiAngle{400, -1.0, 1.0, "#it{#varphi} (rad)"};
    AxisSpec etaTracks{100, -1.0, 1.0, "#it{#eta}_{trk}"};
    AxisSpec detaTracks{400, -0.9, 0.9, "#it{#eta}_{trk}"};
    AxisSpec etaJets{100, -0.6, 0.6, "#it{#eta}_{jet}"};

    AxisSpec jetArea{50, 0.0, 3., "Area_{jet}"};
    AxisSpec rho{50, 0.0, 50., "#it{#rho}"};
    AxisSpec rhoArea{60, 0.0, 60., "#it{#rho} #times Area_{jet}"};

    AxisSpec centrality{10, 0.0, 100., "centrality"};

    // Convert configurable strings to std::string
    std::string evSelToString = static_cast<std::string>(evSel);
    std::string trkSelToString = static_cast<std::string>(trkSel);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(evSelToString);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trkSelToString);

    if (doprocessFilteredCollisions) {
      // Event selection
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{3, 0.0, 3.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of events w/o cuts");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "# of events after z cut"); //?

      // Z coordinate of collision vertex
      spectra.add("hVertexZ_NoCut", "z vertex of collisions w/o cut", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
      spectra.add("hVertexZ_Cut", "z vertex of collisions w. cut", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
      spectra.add("hVertexZ_EventFiltering", "z vertex of collisions w. event filtering", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
    }

    // Disitribution of tracks
    if (doprocessJetsMCDet || doprocessJetsMCDetWeighted || doprocessJets) {
      spectra.add("hPocetTracku", "Pocet", kTH1F, {{1000, 0., 1000., "Pocet Eventu"}});
      spectra.add("hPocetSignalTriggeru", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hPocetReferenceTriggeru", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hRecoilJetRefPt", "#it{p}_{T} distribution of recoil jets", kTH2F, {{200, 0., 200., "Pocet Eventu"}, centrality});
      spectra.add("hRecoilJetSigPt", "#it{p}_{T} distribution of recoil jets", kTH2F, {{200, 0., 200., "Pocet Eventu"}, centrality});
      spectra.add("hRecoilJetRefCorrPt", "Corrected #it{p}_{T} distribution of recoil jets", kTH2F, {{200, 0., 200., "Pocet Eventu"}, centrality});
      spectra.add("hRecoilJetSigCorrPt", "Corrected #it{p}_{T} distribution of recoil jets", kTH2F, {{200, 0., 200., "Pocet Eventu"}, centrality});
      spectra.add("hTrackPt", "#it{p}_{T} distribution of tracks", kTH2F, {pT, centrality});
      spectra.add("hTrackPhi", "#varphi distribution of tracks", kTH1F, {phiAngle});
      spectra.add("hTrackEta", "#eta distribution of tracks", kTH1F, {etaTracks});

      spectra.add("hTrackPtPhi", "#it{p}_{T} vs. #varphi distribution of tracks", kTH2F, {pT, phiAngle});
      spectra.add("hTrackPtEta", "#it{p}_{T} vs. #eta distribution of tracks", kTH2F, {pT, etaTracks});

      spectra.add("hTrackPtPhiEta", "#it{p}_{T} vs. #varphi vs. #eta distribution of tracks", kTH3F, {pT, phiAngle, etaTracks});

      // Distribution of jets
      spectra.add("hJetPt", "#it{p}_{T} distribution of jets", kTH1F, {pT});
      spectra.add("hJetPhi", "#varphi distribution of jets", kTH1F, {phiAngle});
      spectra.add("hJetEta", "#eta distribution of jets", kTH1F, {etaJets});

      spectra.add("hJetPtPhi", "#it{p}_{T} vs. #varphi distribution of jets", kTH2F, {pT, phiAngle});
      spectra.add("hJetPtEta", "#it{p}_{T} vs. #eta distribution of jets", kTH2F, {pT, etaJets});

      spectra.add("hJetPtPhiEta", "#it{p}_{T} vs. #varphi vs. #eta distribution of jets", kTH3F, {pT, phiAngle, etaJets});

      spectra.add("hMultFT0A", "Mult. signal from FTOA", kTH2F, {{2000, 0.0, 40000., "FT0A"}, centrality});
      spectra.add("hMultFT0C", "Mult. signal from FTOC", kTH2F, {{2000, 0.0, 40000., "FT0C"}, centrality});

      spectra.add("hJetPtCorr", "#it{p}_{T} distribution of jets", kTH1F, {pT});

      spectra.add("hRho", "distribution of rho", kTH1F, {rho});
      spectra.add("hjetArea", "distribution of jet area", kTH1F, {jetArea});

      spectra.add("JetAreavsPt", "#it{p}_{T} vs. jet area distribution of jets", kTH2F, {pT, jetArea});
      spectra.add("hRhovsNTracks", "Rho vs number of tracks", kTH2F, {{2000, 0.0, 40000., "NTracks"}, rho});
      spectra.add("hRhovsCentrality", "Rho vs number of tracks", kTH2F, {centrality, rho});

      spectra.add("hSigCorrConePt", "distribution of Signal Corrected perpendicular Cone Pt", kTH1F, {pT});
      spectra.add("hRefCorrConePt", "distribution of Reference Corrected perpendicular Cone Pt", kTH1F, {pT});

      spectra.add("hSigRandCorrConePt", "distribution of Signal Corrected random Cone Pt", kTH2F, {pT, centrality});
      spectra.add("hRefRandCorrConePt", "distribution of Reference Corrected random Cone Pt", kTH2F, {pT, centrality});

      spectra.add("hTTCountData", "Pocet", kTH2F, {{2, 0., 2., "Pocet TT Data"}, centrality});
    }

    if (doprocessJetsMatched) {
      spectra.add("hDetLevelInclusiveJetsPt", "All reconstructed inclusive jets", kTH1F, {{200, 0.0, 200.}}, setSumw2);
      spectra.add("hFakeInclusiveJetsPt", "Det. level inclusive jets w/o matched pair", kTH1F, {{200, 0.0, 200.}}, setSumw2);

      spectra.add("hMissedJets_pT", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}}, setSumw2);
      spectra.add("hMissedJets_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}}, setSumw2);

      spectra.add("hJetPt_resolution", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT}, setSumw2);
      spectra.add("hJetPhi_resolution", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT}, setSumw2);

      spectra.add("hJetPt_DetLevel_vs_PartLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}}, setSumw2);
      spectra.add("hJetPt_resolution_RecoilJets", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT}, setSumw2);
      spectra.add("hJetPhi_resolution_RecoilJets", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT}, setSumw2);
    }

    if (doprocessMCPartLevel || doprocessMCPartLevelWeighted) {
      spectra.add("hRecoilJetRefPtMCP", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "p_{T,ch jet} (GeV/c)"}});
      spectra.add("hRecoilJetSigPtMCP", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "p_{T,ch jet} (GeV/c)"}});
      spectra.add("hRecoilJetRefCorrPtMCP", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "p_{T,ch jet} (GeV/c)"}});
      spectra.add("hRecoilJetSigCorrPtMCP", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "p_{T,ch jet} (GeV/c)"}});

      spectra.add("hEventSelectionCountPartLevel", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hMCPTrackPt", "#it{p}_{T} distribution of tracks", kTH1F, {pT});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, etaJets, phiAngle, rho, jetArea});
      spectra.add("hJetPtMCP", "#it{p}_{T} distribution of jets", kTH1F, {{200, 0., 200., "p_{T,ch jet} (GeV/c)"}});

      spectra.add("hTTCountPartLevel", "Pocet", kTH1F, {{2, 0., 2., "Pocet TT MCP"}});
    }

    if (doprocessTrackEfficiency || doprocessTrackEfficiencyMB) {
      spectra.add("hTrackingEfficiencyDenominatorPt", "#it{p}_{T} distribution", kTH2F, {{200, 0., 200., "p_{T,ch} (GeV/c)"}, centrality});
      spectra.add("hNonprimaryParticlesPt", "#it{p}_{T} distribution", kTH1F, {{200, 0., 200., "p_{T,ch} (GeV/c)"}});
      spectra.add("hTrackingEfficiencyNumeratorPt", "#it{p}_{T} distribution", kTH2F, {{200, 0., 200., "p_{T,ch} (GeV/c)"}, centrality});

      spectra.add("hTrackingPhiSmearingPt", "#it{p}_{T} distribution", kTH2F, {pT, dphiAngle});
      spectra.add("hTrackingEtaSmearingPt", "#it{p}_{T} distribution", kTH2F, {pT, detaTracks});
      spectra.add("hTrackingPtResolutionPt", "#it{p}_{T} distribution", kTH2F, {pTTrack, pTRes});

      spectra.add("hNTracksPerCollision", "Pocet", kTH1F, {{1000, 0., 1000., "Pocet Eventu"}});

      spectra.add("CentralityC", "#it{p}_{T} distribution", kTH1F, {{100, 0., 100., "p_{T,ch} (GeV/c)"}});
      spectra.add("CentralityM", "#it{p}_{T} distribution", kTH1F, {{100, 0., 100., "p_{T,ch} (GeV/c)"}});
      spectra.add("AllTracksPtCent", "#it{p}_{T} distribution", kTH2F, {{200, 0., 200., "p_{T,ch} (GeV/c)"}, centrality});
    }
  }

  //------------------------------------------------------------------------------

  template <typename Collision, typename Jets, typename Tracks>
  void fillHistograms(Collision const& collision, Jets const& jets,
                      Tracks const& tracks, float weight = 1.)

  {

    const float randomConeDistanceFactor = 2.;
    const int maxRandomConeIterations = 1000;

    std::vector<float> sigTrackPhi;
    std::vector<float> sigTrackEta;
    std::vector<float> refTrackPhi;
    std::vector<float> refTrackEta;

    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float centralityM = collision.centFT0M();

    spectra.fill(HIST("hMultFT0A"), multFT0A, centralityM, weight);
    spectra.fill(HIST("hMultFT0C"), multFT0C, centralityM, weight);

    float jetR2 = jetR * jetR;
    float rho = collision.rho();

    spectra.fill(HIST("hRho"), rho, weight);

    int i = 0;
    for (const auto& track : tracks) {
      // check whether track passes the selection flags
      if (skipTrack(track))
        continue;
      i++;
      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();

      spectra.fill(HIST("hTrackPt"), trackPt, centralityM, weight);
      spectra.fill(HIST("hTrackPhi"), trackPhi, weight);
      spectra.fill(HIST("hTrackEta"), trackEta, weight);

      spectra.fill(HIST("hTrackPtPhi"), trackPt, trackPhi, weight);
      spectra.fill(HIST("hTrackPtEta"), trackPt, trackEta, weight);

      spectra.fill(HIST("hTrackPtPhiEta"), trackPt, trackPhi, trackEta, weight);

      if ((signalTriggerMin <= trackPt) && (trackPt <= signalTriggerMax)) {
        sigTrackPhi.push_back(trackPhi);
        sigTrackEta.push_back(trackEta);
      }

      if ((referenceTriggerMin <= trackPt) && (trackPt <= referenceTriggerMax)) {
        refTrackPhi.push_back(trackPhi);
        refTrackEta.push_back(trackEta);
      }
    }
    spectra.fill(HIST("hRhovsNTracks"), i, rho, weight);
    spectra.fill(HIST("hRhovsCentrality"), centralityM, rho, weight);

    float perpConePhi = -99;
    float phiTT = -999; // this will trigger track phi
    float etaTT = -999; // trigger track eta
    float rnd = rand->Rndm();
    bool analyzeSignal = 0; // 0= reference TT  and 1= signal TT

    if (rnd < sigToRefFraction) {
      analyzeSignal = 1;
    }

    if (analyzeSignal == 0 && refTrackPhi.size() > 0) {

      int ii = rand->Integer(refTrackPhi.size());
      spectra.fill(HIST("hTTCountData"), 0.5, centralityM, weight);
      spectra.fill(HIST("hPocetReferenceTriggeru"), static_cast<float>(refTrackPhi.size()), weight);
      phiTT = refTrackPhi[ii];
      etaTT = refTrackEta[ii];

    } else if (analyzeSignal == 1 && sigTrackPhi.size() > 0) {

      int ii = rand->Integer(sigTrackPhi.size());
      spectra.fill(HIST("hTTCountData"), 1.5, centralityM, weight);
      spectra.fill(HIST("hPocetSignalTriggeru"), static_cast<float>(sigTrackPhi.size()), weight);
      phiTT = sigTrackPhi[ii];
      etaTT = sigTrackEta[ii];
    }

    spectra.fill(HIST("hPocetTracku"), i, weight);

    for (const auto& jet : jets) {
      // jet spectra before TT selection
      auto jetPt = jet.pt();

      auto jetPhi = jet.phi();
      auto jetEta = jet.eta();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtCorr"), jetPtCorr, weight);
      spectra.fill(HIST("hJetPt"), jetPt, weight);
      spectra.fill(HIST("hJetPhi"), jetPhi, weight);
      spectra.fill(HIST("hJetEta"), jetEta, weight);

      spectra.fill(HIST("hJetPtPhi"), jetPt, jetPhi, weight);
      spectra.fill(HIST("hJetPtEta"), jetPt, jetEta, weight);

      spectra.fill(HIST("hJetPtPhiEta"), jetPt, jetPhi, jetEta, weight);

      spectra.fill(HIST("hjetArea"), jetArea, weight);
      spectra.fill(HIST("JetAreavsPt"), jetPt, jetArea, weight);
    }

    const int minPhiTT = -100; // phi for TT not found
    if (phiTT < minPhiTT) {
      return; // trigger track was not found
    }

    perpConePhi = TVector2::Phi_mpi_pi(phiTT - constants::math::PIHalf);

    float conePt = 0.;
    ROOT::Math::PtEtaPhiMVector sumlv(0., 0., 0., 0.);

    for (const auto& track : tracks) {
      // check whether track passes the selection flags
      if (skipTrack(track))
        continue;

      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();

      float dPhiPerp = TVector2::Phi_mpi_pi(trackPhi - perpConePhi);
      float dEtaPerp = trackEta - etaTT;

      if ((dPhiPerp) * (dPhiPerp) + (dEtaPerp) * (dEtaPerp) < jetR2) {
        ROOT::Math::PtEtaPhiMVector lv(trackPt, trackEta, trackPhi, pionMass);
        sumlv += lv;
      }
    }
    conePt = static_cast<float>(sumlv.Pt());
    // calculate delta pT from perpendicular cone
    float corrConeDeltaPt = conePt - constants::math::PI * jetR2 * rho;
    if (analyzeSignal == 0) {
      spectra.fill(HIST("hRefCorrConePt"), corrConeDeltaPt, weight);
    } else {
      spectra.fill(HIST("hSigCorrConePt"), corrConeDeltaPt, weight);
    }

    int jj = 0;
    const int jet1 = 1;
    const int jet2 = 2;
    float jet1Phi = 0.;
    float jet1Eta = 0.;
    float jet2Phi = 0.;
    float jet2Eta = 0.;

    for (const auto& jet : jets) {

      jj++;

      auto jetPhi = jet.phi();
      auto jetEta = jet.eta();

      if (jj == jet1) {
        jet1Phi = jetPhi;
        jet1Eta = jetEta;
      }

      if (jj == jet2) {
        jet2Phi = jetPhi;
        jet2Eta = jetEta;
      }
    }

    bool bCloseJet = true;
    float randJetPhi = -999.;
    float randJetEta = -999.;
    int iii = 0;

    while (bCloseJet) {
      iii++;
      if (iii > maxRandomConeIterations) {
        break;
      }

      randJetPhi = rand->Uniform(-constants::math::PI, constants::math::PI);
      randJetEta = rand->Uniform(-0.5, 0.5);

      float vphi1 = TVector2::Phi_mpi_pi(randJetPhi - jet1Phi);
      float vphi2 = TVector2::Phi_mpi_pi(randJetPhi - jet2Phi);
      float vphiTT = TVector2::Phi_mpi_pi(randJetPhi - phiTT);
      float veta1 = randJetEta - jet1Eta;
      float veta2 = randJetEta - jet2Eta;
      float vetaTT = randJetEta - etaTT;
      float dist1 = std::sqrt(vphi1 * vphi1 + veta1 * veta1);
      float dist2 = std::sqrt(vphi2 * vphi2 + veta2 * veta2);
      float distTT = std::sqrt(vphiTT * vphiTT + vetaTT * vetaTT);

      if (dist1 > randomConeDistanceFactor * jetR &&
          dist2 > randomConeDistanceFactor * jetR &&
          distTT > randomConeDistanceFactor * jetR)
        bCloseJet = false;
    }

    float randConePt = 0.;
    sumlv = ROOT::Math::PtEtaPhiMVector(0., 0., 0., 0.);
    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();
      float dPhiRand = TVector2::Phi_mpi_pi(trackPhi - randJetPhi);
      float dEtaRand = trackEta - randJetEta;

      if ((dPhiRand) * (dPhiRand) + (dEtaRand) * (dEtaRand) < jetR2) {
        ROOT::Math::PtEtaPhiMVector lv(trackPt, trackEta, trackPhi, pionMass);
        sumlv += lv;
      }
    }
    randConePt = static_cast<float>(sumlv.Pt());

    float corrRandConeDeltaPt = randConePt - constants::math::PI * jetR2 * rho;

    if (analyzeSignal == 0) {
      spectra.fill(HIST("hRefRandCorrConePt"), corrRandConeDeltaPt, centralityM, weight);
    } else {
      spectra.fill(HIST("hSigRandCorrConePt"), corrRandConeDeltaPt, centralityM, weight);
    }

    // Recoil jet spectra
    for (const auto& jet : jets) {
      auto jetPt = jet.pt();
      auto jetPhi = jet.phi();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      float dPhi = -100;
      dPhi = std::fabs(TVector2::Phi_mpi_pi(jetPhi - phiTT));

      if (dPhi > constants::math::PI - dPhiCut) {
        if (analyzeSignal == 0) {
          spectra.fill(HIST("hRecoilJetRefPt"), jetPt, centralityM, weight);
          spectra.fill(HIST("hRecoilJetRefCorrPt"), jetPtCorr, centralityM, weight);
        } else {
          spectra.fill(HIST("hRecoilJetSigPt"), jetPt, centralityM, weight);
          spectra.fill(HIST("hRecoilJetSigCorrPt"), jetPtCorr, centralityM, weight);
        }
      }
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // template <typename Collision, typename Jets, typename Particles>
  template <typename JTracksTable, typename JetsBase, typename JetsTag>
  void fillMatchedHistograms(JTracksTable const& tracks,
                             JetsBase const& jetsBase,
                             JetsTag const& jetsTag,
                             float weight = 1.)
  {

    std::vector<double> sigTrackPhi;
    double phiTT = 0.0;

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;
      if ((signalTriggerMin <= track.pt()) && (track.pt() <= signalTriggerMax))
        sigTrackPhi.push_back(track.phi());
    }

    bool bIsThereTTSig = sigTrackPhi.size() > 0;

    if (bIsThereTTSig)
      phiTT = sigTrackPhi[rand->Integer(sigTrackPhi.size())];

    for (const auto& jetBase : jetsBase) {
      bool bIsBaseJetRecoil = false; ///////////////////////////////////////////////////////////////////////////////////////////////////

      if (bIsThereTTSig) {
        float dPhi = std::fabs(TVector2::Phi_mpi_pi(jetBase.phi() - phiTT));
        if (dPhi > constants::math::PI - dPhiCut)
          bIsBaseJetRecoil = true;
      }

      dataForUnfolding(jetBase, jetsTag, bIsBaseJetRecoil, tracks, weight);
    }
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processFilteredCollisions(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hVertexZ_EventFiltering"), collision.posZ()); // compare with "hVertexZ_Cut" histogram
  }

  PROCESS_SWITCH(TrackJetSpectra, processFilteredCollisions, "process filtered collisions", false);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processJets(FilteredColl const& collision,
                   soa::Filtered<aod::ChargedJets> const& jets,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (skipEvent(collision)) {
      return;
    }

    // check whether event contains the required flag for event selection (sel8 in this case)

    fillHistograms(collision, jets, tracks);
  }

  PROCESS_SWITCH(TrackJetSpectra, processJets, "process inclusive jets", false);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processJetsMCDet(FilteredColl const& collision,
                        FilteredJetsDetLevel const& jets,
                        soa::Filtered<aod::JetTracks> const& tracks)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision)) {
      return;
    }
    fillHistograms(collision, jets, tracks);
  }

  PROCESS_SWITCH(TrackJetSpectra, processJetsMCDet, "process MCDet jets", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processJetsMCDetWeighted(FilteredCollDetLevelGetWeight const& collision,
                                FilteredJetsDetLevel const& jets,
                                soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (skipEvent(collision) || collision.isOutlier()) {
      return;
    }
    // check whether event contains the required flag for event selection (sel8 in this case)
    auto weight = collision.weight();
    fillHistograms(collision, jets, tracks, weight);
  }

  PROCESS_SWITCH(TrackJetSpectra, processJetsMCDetWeighted, "process MCDet jets", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processMCPartLevel(FilteredCollPartLevelMB const& collision,
                          FilteredParticles const& particles,
                          FilteredJetsPartLevel const& jets)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision))
      return;

    fillMCPHistograms(collision, jets, particles);
  }

  PROCESS_SWITCH(TrackJetSpectra, processMCPartLevel, "process MCDet jets", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processMCPartLevelWeighted(FilteredCollPartLevel const& collision,
                                  FilteredParticles const& particles,
                                  FilteredJetsPartLevel const& jets)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision) || collision.isOutlier())
      return;

    auto weight = collision.weight();
    fillMCPHistograms(collision, jets, particles, weight);
  }

  PROCESS_SWITCH(TrackJetSpectra, processMCPartLevelWeighted, "process MCDet jets", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processTrackEfficiency(CollPartLevelEfficiencyWeighted const& collisionMC,
                              CollisionMCD const& collision,
                              TracksMCDetFullTable const& tracks,
                              ParticlesFullTable const& mcParticles)
  {
    if (collisionMC.isOutlier())
      return;
    auto weight = collisionMC.weight();

    fillEfficiency(collisionMC,
                   collision,
                   tracks,
                   mcParticles,
                   weight);
  }

  PROCESS_SWITCH(TrackJetSpectra, processTrackEfficiency, "process efficiency of tracks", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processTrackEfficiencyMB(CollPartLevelMBEfficiency const& collisionMC,
                                CollisionMCD const& collision,
                                TracksMCDetFullTable const& tracks,
                                ParticlesFullTable const& mcParticles)
  {
    // if (skipEvent(collision)) return;

    fillEfficiency(collisionMC,
                   collision,
                   tracks,
                   mcParticles);
  }
  PROCESS_SWITCH(TrackJetSpectra, processTrackEfficiencyMB, "process efficiency of tracks", false);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename JCollisionMC, typename JCollision, typename JTracks, typename JMcParticles>
  void fillEfficiency(JCollisionMC const& collisionMC,
                      JCollision const& collisions,
                      JTracks const& tracks,
                      JMcParticles const& mcParticles,
                      float weight = 1.)

  {

    if (!(std::fabs(collisionMC.posZ()) < vertexZCut)) {
      return;
    }

    if (collisions.size() < 1) {
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
      return;
    }

    float centralityC = collisions.begin().centFT0C();

    float centralityM = collisions.begin().centFT0M();

    spectra.fill(HIST("CentralityC"), centralityC, weight);

    spectra.fill(HIST("CentralityM"), centralityM, weight);

    for (auto const& particle : mcParticles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      if (static_cast<int8_t>(pdgParticle->Charge()) == 0)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;

      if (std::fabs(particle.eta()) > trkEtaCut)
        continue;
      if (particle.pt() < trkPtMin)
        continue;

      spectra.fill(HIST("hTrackingEfficiencyDenominatorPt"), particle.pt(), centralityM, weight);
    }

    for (auto const& collision : collisions) {
      if (skipEvent(collision))
        continue;
      if (!(std::fabs(collision.posZ()) < vertexZCut))
        continue;

      auto collTracks = tracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      spectra.fill(HIST("hNTracksPerCollision"), collTracks.size());
      // for (auto const& particle : particles){

      for (auto const& track : collTracks) {
        if (skipTrack(track))
          continue;

        if (std::fabs(track.eta()) > trkEtaCut)
          continue;
        if (track.pt() < trkPtMin)
          continue;

        spectra.fill(HIST("AllTracksPtCent"), track.pt(), centralityM, weight);

        if (!track.has_mcParticle())
          continue;

        auto mcPart = track.template mcParticle_as<ParticlesFullTable>();

        auto pdgParticle = pdg->GetParticle(mcPart.pdgCode());
        if (static_cast<int8_t>(pdgParticle->Charge()) == 0)
          continue;

        if (std::fabs(mcPart.eta()) > trkEtaCut)
          continue;
        if (mcPart.pt() < trkPtMin)
          continue;

        if (!mcPart.isPhysicalPrimary()) {
          spectra.fill(HIST("hNonprimaryParticlesPt"), track.pt(), weight);
          continue;
        }
        spectra.fill(HIST("hTrackingEfficiencyNumeratorPt"), mcPart.pt(), centralityM, weight);

        float dPhi = TVector2::Phi_mpi_pi(mcPart.phi() - track.phi());
        float dEta = mcPart.eta() - track.eta();
        float dPtResolution = (mcPart.pt() - track.pt()) / mcPart.pt();

        spectra.fill(HIST("hTrackingPhiSmearingPt"), mcPart.pt(), dPhi, weight);
        spectra.fill(HIST("hTrackingEtaSmearingPt"), mcPart.pt(), dEta, weight);
        spectra.fill(HIST("hTrackingPtResolutionPt"), mcPart.pt(), dPtResolution, weight);
      }
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename Collision, typename Jets, typename Particles>
  void fillMCPHistograms(Collision const& collision, Jets const& jets,
                         Particles const& particles, float weight = 1.)

  {

    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5, weight);

    std::vector<float> vPhiOfTT;

    float rnd = rand->Rndm();
    bool analyzeSignal = 0;
    float rho = collision.rho();

    if (rnd < sigToRefFraction)
      analyzeSignal = 1;

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      // Need charge and physical primary particles
      bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
      if (bParticleNeutral || !particle.isPhysicalPrimary())
        continue;

      float particlePt = particle.pt();
      float particlePhi = particle.phi();

      if (analyzeSignal && (particlePt > signalTriggerMin && particlePt < signalTriggerMax)) {
        vPhiOfTT.push_back(particlePhi);
      }

      if (!analyzeSignal && (particlePt > referenceTriggerMin && particlePt < referenceTriggerMax)) {
        vPhiOfTT.push_back(particlePhi);
      }
    }

    int ii = -1; // index of selected particle level trigger track
    if (vPhiOfTT.size() > 0) {
      ii = rand->Integer(vPhiOfTT.size());
      if (analyzeSignal == 0) {
        spectra.fill(HIST("hTTCountPartLevel"), 0.5, weight);
      } else {
        spectra.fill(HIST("hTTCountPartLevel"), 1.5, weight);
      }
    }

    for (const auto& jet : jets) {
      float jetPt = jet.pt();
      float jetPhi = jet.phi();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea_Part"), jetPt, jet.eta(), jet.phi(), rho, jetArea, weight);
      spectra.fill(HIST("hJetPtMCP"), jetPt, weight);

      if (ii > -1) { // only events with TT particle level

        float phiTT = vPhiOfTT[ii];

        float deltaPhi = std::fabs(TVector2::Phi_mpi_pi(phiTT - jetPhi));
        if (deltaPhi > constants::math::PI - dPhiCut) {
          if (analyzeSignal == 0) {
            spectra.fill(HIST("hRecoilJetRefPtMCP"), jetPt, weight);
            spectra.fill(HIST("hRecoilJetRefCorrPtMCP"), jetPtCorr, weight);
          } else {
            spectra.fill(HIST("hRecoilJetSigPtMCP"), jetPt, weight);
            spectra.fill(HIST("hRecoilJetSigCorrPtMCP"), jetPtCorr, weight);
          }
        }
      }
    }
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void processJetsMatched(FilteredCollDetLevelGetWeight const& collision,
                          aod::JetMcCollisions const&,
                          FilteredTracks const& tracks,
                          FilteredMatchedJetsDetLevel const& mcdjets,
                          FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision) || collision.isOutlier())
      return;

    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets);
  }
  PROCESS_SWITCH(TrackJetSpectra, processJetsMatched, "process matching of MC jets (no weight)", false);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Auxiliary functions

  template <typename Collision>
  bool skipEvent(const Collision& coll)
  {
    /// \brief: trigger cut is needed for pp data
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits, skipMBGapEvents);
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename Collision>
  bool skipMCEvent(const Collision& coll)
  {
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits);
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename Track>
  bool skipTrack(const Track& track)
  {

    // check track quality select global tracks
    return !jetderiveddatautilities::selectTrack(track, trackSelection);
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename Jet, typename Tracks>
  bool isJetWithHighPtConstituent(Jet const& jet, Tracks const&)
  {
    // kill all jets which have a constituent track with pT > 100 GeV
    bool bIsJetWithHighPtConstituent = false;
    for (const auto& jetConstituent : jet.template tracks_as<Tracks>()) {
      if (jetConstituent.pt() > maxJetConstituentPt) {
        bIsJetWithHighPtConstituent = true;
        break;
      }
    }
    return bIsJetWithHighPtConstituent;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <typename PartJet, typename DetJet, typename TracksTable>
  void dataForUnfolding(PartJet const& partJet, DetJet const& detJets, bool bIsBaseJetRecoil, TracksTable const& tracks, float weight = 1.)
  {
    // function which matches  particle jet with detector level jet
    float partJetPt = partJet.pt();
    bool bIsThereMatchedJet = partJet.has_matchedJetGeo();

    if (bIsThereMatchedJet) {
      const auto& jetsMatched = partJet.template matchedJetGeo_as<std::decay_t<DetJet>>();

      for (const auto& jetMatched : jetsMatched) {
        // skip matches where detector level jets have a constituent with pT above specified cut
        bool skipMatchedDetJet = isJetWithHighPtConstituent(jetMatched, tracks);

        if (skipMatchedDetJet) {
          // Miss jets
          spectra.fill(HIST("hMissedJets_pT"), partJetPt, weight);

          if (bIsBaseJetRecoil)
            spectra.fill(HIST("hMissedJets_pT_RecoilJets"), partJetPt, weight);

        } else {
          float detJetPt = jetMatched.pt();

          spectra.fill(HIST("hJetPt_resolution"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
          spectra.fill(HIST("hJetPhi_resolution"), partJet.phi() - jetMatched.phi(), partJetPt, weight);

          if (bIsBaseJetRecoil) {
            spectra.fill(HIST("hJetPt_DetLevel_vs_PartLevel_RecoilJets"), detJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPt_resolution_RecoilJets"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPhi_resolution_RecoilJets"), partJet.phi() - jetMatched.phi(), partJetPt, weight);
          }
        }
      }
    }

    for (const auto& jetDet : detJets) {
       if(isJetWithHighPtConstituent(jetDet, tracks))
         continue;

       auto detJetPt = jetDet.pt();
       spectra.fill(HIST("hDetLevelInclusiveJetsPt"), detJetPt, weight);
       if (!jetDet.has_matchedJetGeo()) {
         spectra.fill(HIST("hFakeInclusiveJetsPt"), detJetPt, weight);
       }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackJetSpectra>(cfgc)};
}
