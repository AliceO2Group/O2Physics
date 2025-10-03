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

// QA task for MC-based outliers
//
/// \author Jaime Norman <jaime.norman@cern.ch>

#include "RecoDecay.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>
#include <TMathBase.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetOutlierQATask {

  HistogramRegistry registry;

  Preslice<aod::JetTracks> perCol = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> perColMC = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>> perColJets = aod::jet::collisionId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>> perColJetsMatched = aod::jet::collisionId;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};
  Configurable<double> jetPtMax{"jetPtMax", 200., "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "number of bins for eta axes"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  
  Configurable<float> splitCollisionsDeltaZ{"splitCollisionsDeltaZ", 0.2, "threshold in delta z to assign as split collision"};
  Configurable<float> splitCollisionsDeltaZPart{"splitCollisionsDeltaZPart", 0.2, "threshold in delta z to assign as split collision particle level"};
  Configurable<int> splitCollisionsDeltaBC{"splitCollisionsDeltaBC", 5, "threshold in BC to assign as split collision"};
  Configurable<int> mergeCollisionsDeltaMin{"mergeCollisionsDeltaMin", -10, "number of prior collisions to search for close Z position"};
  Configurable<int> mergeCollisionsDeltaMax{"mergeCollisionsDeltaMax", 10, "number of following collisions to search for close Z position"};

  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices

  std::vector<bool> filledJetR_Both;
  std::vector<bool> filledJetR_Low;
  std::vector<bool> filledJetR_High;
  std::vector<double> jetRadiiValues;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    jetRadiiValues = (std::vector<double>)jetRadii;

    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR_Both.push_back(0.0);
      filledJetR_Low.push_back(0.0);
      filledJetR_High.push_back(0.0);
    }
    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
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

    AxisSpec jetEtaAxis = {nBinsEta, jetEtaMin, jetEtaMax, "#eta"};
    AxisSpec trackEtaAxis = {nBinsEta, trackEtaMin, trackEtaMax, "#eta"};

    if (doprocessJetsAmbiguous ) {
      registry.add("h3_jet_pt_track_pt_pt_hat_ambiguous", "ambiguous;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {100,0,100}, {200,0,600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_unambiguous", "matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {100,0,100}, {200,0,600}}});
      registry.add("h3_jet_pt_frac_pt_ambiguous_pt_hat", "fraction pT;#it{p}_{T,jet} (GeV/#it{c});fraction of #it{p}_{T,track} unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {40,0,1.1}, {200,0,600}}});
      registry.add("h3_jet_pt_frac_constituents_ambiguous_pt_hat", "fraction const;#it{p}_{T,jet} (GeV/#it{c});fraction of constituents matched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {40,0,1.1}, {200,0,600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_no_particle", "no matching particle;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {100,0,100}, {200,0,600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_with_particle", "with matching particle;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {100,0,100}, {200,0,600}}});
      registry.add("h3_jet_pt_frac_pt_unmatched_particle_pt_hat", "fraction pT;#it{p}_{T,jet} (GeV/#it{c});fraction of #it{p}_{T,track} unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {40,0,1.1}, {200,0,600}}});
      registry.add("h3_jet_pt_frac_constituents_unmatched_particle_pt_hat", "fraction const;#it{p}_{T,jet} (GeV/#it{c});fraction of constituents unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150,0,300}, {40,0,1.1}, {200,0,600}}});
    }
    if(doprocessCollisionsBC) {
      registry.add("h_DeltaZ_InBunch", "Delta Z between two events in bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_InBunch", "Delta Z between two events in bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_InBunch", "Z1 vs Z2 between two events in bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_OutOfBunch", "Delta Z between two events out of bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_OutOfBunch", "Delta Z between two events out of bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_OutOfBunch", "Z1 vs Z2 between two events out of bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_InBunch_JJ", "Delta Z between two events in bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_InBunch_JJ", "Delta Z between two events in bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_InBunch_JJ", "Z1 vs Z2 between two events in bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_OutOfBunch_JJ", "Delta Z between two events out of bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_OutOfBunch_JJ", "Delta Z between two events out of bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_OutOfBunch_JJ", "Z1 vs Z2 between two events out of bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_Z", "Delta Z between two events", {HistType::kTH1F, {{400, -12, 12}}});
    }
    if(doprocessTracksBC) {
      registry.add("h_Z_resolution", "Z resolution", {HistType::kTH1F, {{200, -0.05, 0.05}}});
      registry.add("h_Z_resolution_wide", "Z resolution", {HistType::kTH1F, {{200, -1, 1}}});
      registry.add("h_Z_reco_rejected", "Z position reconstructed rejected", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_true_rejected", "Z position particle rejected", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_reco_accepted", "Z position reconstructed accepted", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_true_accepted", "Z position particle accepted", {HistType::kTH1F, {{400, -12, 12}}});

      registry.add("h_track_pt", "track pt;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0, 300}}});
      registry.add("h_track_eta", "track eta;#eta_{track};entries", {HistType::kTH1F, {{100, -5, 5}}});
      registry.add("h_track_phi", "track phi;#varphi_{track} (rad);entries", {HistType::kTH1F, {{160, -1.0, 7.0}}});
      registry.add("h_track_pt_eta", "track pt vs eta;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{300, 0, 300}, {100, -5, 5}}});
      registry.add("h_track_pt_phi", "track pt vs phi;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{300, 0, 300}, {160, -1.0, 7.0}}});

      registry.add("h_track_pt_accepted", "track pt accepted;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0, 300}}});
      registry.add("h_track_eta_accepted", "track eta accepted;#eta_{track};entries", {HistType::kTH1F, {{100, -5, 5}}});
      registry.add("h_track_phi_accepted", "track phi accepted;#varphi_{track} (rad);entries", {HistType::kTH1F, {{160, -1.0, 7.0}}});
      registry.add("h_track_pt_eta_accepted", "track pt vs eta accepted;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{300, 0, 300}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_accepted", "track pt vs phi accepted;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{300, 0, 300}, {160, -1.0, 7.0}}});

      registry.add("h_track_pt_no_collision", "track pt no collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0, 300}}});
      registry.add("h_track_pt_collision", "track pt collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0, 300}}});
      registry.add("h2_track_pt_pt_hat_no_particle", "track pt vs pt hat no particle;p_{T,track} (GeV/#it{c});#hat{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{300, 0, 300}, {600, 0, 600}}});
      registry.add("h2_track_pt_pt_hat_particle", "track pt vs pt hat particle;p_{T,track} (GeV/#it{c});#hat{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{300, 0, 300}, {600, 0, 600}}});

      registry.add("h_track_pt_outlier", "weight track pt", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h2_pt_hat_track_pt", "track; #hat{#it{p}_{T}} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{600, 0, 600}, {150, 0, 300}}});
      registry.add("h2_pt_hat_track_pt_outlier", "track; #hat{#it{p}_{T}} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{600, 0, 600}, {150, 0, 300}}});
      registry.add("h2_neighbour_pt_hat_outlier", "neighbour; distance from collision; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {600, 0, 600}}});
      registry.add("h2_neighbour_track_pt_outlier", "neighbour; distance from collision; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {200, 0, 100}}});
      registry.add("h2_neighbour_pt_hat_all", "neighbour; distance from collision; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {600, 0, 600}}});
      registry.add("h2_neighbour_track_pt_all", "neighbour; distance from collision; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {200, 0, 100}}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);
  
  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
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
 
  void fillHistogramsAmbiguous(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>::iterator const& jet, float weight, aod::AmbiguousTracks const& tracksAmbiguous)
  {
    // take method from qaEventTrack.cxx

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
     // registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
    }

    auto iterAmbiguous = tracksAmbiguous.begin();
    bool searchAmbiguous = true;
    int nAmbTracks = 0;
    int nUnmatchedTracks = 0;
    double pt_total = 0;
    double pt_amb = 0;
    double pt_unmatched = 0;
    //LOG(info) << "===== constituent.size()=" << jet.template tracks_as<aod::JetTracksMCD>().size() << ", tracksAmbiguous.size()=" << tracksAmbiguous.size();
    for (auto& constituent : jet.template tracks_as<aod::JetTracksMCD>()) {
      pt_total += constituent.pt();
      if(!constituent.has_collision()) {
        //LOG(info) << "NO COLLISION : track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex();
      }
      bool has_MCparticle = constituent.has_mcParticle();
      if(!has_MCparticle) {
        //LOG(info) << "constituent NO MC PARTICLE: track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex();
        registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_no_particle"), jet.pt(), constituent.pt(), pTHat, weight);
        pt_unmatched += constituent.pt();
        nUnmatchedTracks++;
      }
      else {
        //LOG(info) << "constituent HAS MC PARTICLE: track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex();
        registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_with_particle"), jet.pt(), constituent.pt(), pTHat, weight);
      }

      bool goFillHisto = (iterAmbiguous != tracksAmbiguous.end());
      if (goFillHisto) {
        while (constituent.globalIndex() > iterAmbiguous.trackId()) {
          //LOG(info) << "track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex() << " ||| iterAmbiguous.trackId()=" << iterAmbiguous.trackId() << ", iterAmbiguous.globalIndex()=" << iterAmbiguous.globalIndex();
          iterAmbiguous++;
          if (iterAmbiguous == tracksAmbiguous.end()) { /// all ambiguous tracks found
            searchAmbiguous = false;
            goFillHisto = false;
            break;
          }
        }
      }
      if (goFillHisto) {
        if (constituent.globalIndex() == iterAmbiguous.trackId()) {
          nAmbTracks++;
          //LOG(info) << "   >>> FOUND AMBIGUOUS TRACK! track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex() << " ||| iterAmbiguous.trackId()=" << iterAmbiguous.trackId() << ", iterAmbiguous.globalIndex()=" << iterAmbiguous.globalIndex();
          registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_ambiguous"), jet.pt(), constituent.pt(), pTHat, weight);
          pt_amb += constituent.pt();
          nAmbTracks++;
        }
        else {
          registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_unambiguous"), jet.pt(), constituent.pt(), pTHat, weight);
        }
      }
    }
    registry.fill(HIST("h3_jet_pt_frac_pt_ambiguous_pt_hat"), jet.pt(), pt_amb/pt_total, pTHat, weight);
    registry.fill(HIST("h3_jet_pt_frac_constituents_ambiguous_pt_hat"), jet.pt(), double(nAmbTracks)/double(jet.template tracks_as<aod::JetTracksMCD>().size()), pTHat, weight);

    registry.fill(HIST("h3_jet_pt_frac_pt_unmatched_particle_pt_hat"), jet.pt(), pt_unmatched/pt_total, pTHat, weight);
    registry.fill(HIST("h3_jet_pt_frac_constituents_unmatched_particle_pt_hat"), jet.pt(), double(nUnmatchedTracks)/double(jet.template tracks_as<aod::JetTracksMCD>().size()), pTHat, weight);
  }

void processJetsAmbiguous(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets, aod::JetTracksMCD const&, const aod::AmbiguousTracks& tracksAmbiguous)
  {
    // 
    // jet-based outlier checks based on ambiguous tracks
    //
    if (collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramsAmbiguous(jet, jet.eventWeight(), tracksAmbiguous);
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processJetsAmbiguous, "jet finder QA mcd with weighted events", false);

  void processCollisionsBC(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionBCs> const& collisions, 
      aod::JetMcCollisions const&,
      aod::JetTracks const& tracks
      )
  {
    //
    // collision-based outlier checks based on BC and z position
    // based on 2-event correlation checks in PWGDQ/Tasks/tableReader_withAssoc.cxx
    //

    fBCCollMap.clear();
    for (auto const& collision : collisions) {
      // Fill the BC map of events
      if (fBCCollMap.find(collision.bcId()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {collision.globalIndex()};
        fBCCollMap[collision.bcId()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[collision.bcId()];
        evIndices.push_back(collision.globalIndex());
      }
      registry.fill(HIST("h_Z"), collision.posZ());
    }

    // Create a map for collisions which are candidate of being split
    // key: event global index, value: whether pileup event is a possible splitting
    // (not used ATM, but could be used to flag events in the future)
    std::map<int64_t, bool> collisionSplittingMap;

    // loop over the BC map, get the collision vectors and make in-bunch and out of bunch 2-event correlations
    for (auto bc1It = fBCCollMap.begin(); bc1It != fBCCollMap.end(); ++bc1It) {
      uint64_t bc1 = bc1It->first;
      auto bc1Events = bc1It->second;

      // same bunch event correlations, if more than 1 collisions in this bunch
      if (bc1Events.size() > 1) {
        for (auto ev1It = bc1Events.begin(); ev1It != bc1Events.end(); ++ev1It) {
          auto ev1 = collisions.rawIteratorAt(*ev1It);
          for (auto ev2It = std::next(ev1It); ev2It != bc1Events.end(); ++ev2It) {
            auto ev2 = collisions.rawIteratorAt(*ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            float deltaZ = ev1.posZ() - ev2.posZ();
            if (TMath::Abs(deltaZ) < splitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[*ev1It] = true;
              collisionSplittingMap[*ev2It] = true;
            }
            registry.fill(HIST("h_DeltaZ_InBunch"), deltaZ);
            registry.fill(HIST("h_DeltaZ_Z1_InBunch"), deltaZ, ev1.posZ());
            registry.fill(HIST("h_Z1_Z2_InBunch"), ev1.posZ(), ev2.posZ());
            if (ev1.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap &&
                ev2.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) { // both are non-gap events
              registry.fill(HIST("h_DeltaZ_InBunch_JJ"), deltaZ);
              registry.fill(HIST("h_DeltaZ_Z1_InBunch_JJ"), deltaZ, ev1.posZ());
              registry.fill(HIST("h_Z1_Z2_InBunch_JJ"), ev1.posZ(), ev2.posZ());
            }
          } // end second event loop
        } // end first event loop
      } // end if BC1 events > 1

      // loop over the following BCs in the TF
      for (auto bc2It = std::next(bc1It); bc2It != fBCCollMap.end(); ++bc2It) {
        uint64_t bc2 = bc2It->first;
        if ((bc2 > bc1 ? bc2 - bc1 : bc1 - bc2) > splitCollisionsDeltaBC) {
          break;
        }
        auto bc2Events = bc2It->second;

        // loop over events in the first BC
        for (auto ev1It : bc1Events) {
          auto ev1 = collisions.rawIteratorAt(ev1It);
          // loop over events in the second BC
          for (auto ev2It : bc2Events) {
            auto ev2 = collisions.rawIteratorAt(ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            float deltaZ = ev1.posZ() - ev2.posZ();
            if (TMath::Abs(deltaZ) < splitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[ev1It] = true;
              collisionSplittingMap[ev2It] = true;
            }
            registry.fill(HIST("h_DeltaZ_OutOfBunch"), deltaZ);
            registry.fill(HIST("h_DeltaZ_Z1_OutOfBunch"), deltaZ, ev1.posZ());
            registry.fill(HIST("h_Z1_Z2_OutOfBunch"), ev1.posZ(), ev2.posZ());
            if (ev1.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap &&
                ev2.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) { // both are non-gap events
              registry.fill(HIST("h_DeltaZ_OutOfBunch_JJ"), deltaZ);
              registry.fill(HIST("h_DeltaZ_Z1_OutOfBunch_JJ"), deltaZ, ev1.posZ());
              registry.fill(HIST("h_Z1_Z2_OutOfBunch_JJ"), ev1.posZ(), ev2.posZ());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processCollisionsBC, "jet finder QA outliers", false);

  void processTracksBC(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionBCs>> const& collisions, 
      aod::JetMcCollisions const& collisionsMC,
      aod::JetTracksMCD const& tracks
      )
  {
    //
    // track-based outlier checks
    //
    
    // first check for collisions occuring close by in time and z in MC
    std::set<int> closeByCollisionIDs;  
    for (auto const& collisionMC : collisionsMC) {
      if (collisionMC.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        return;
      }
      float posZtrue = collisionMC.posZ();
      for (auto const& collisionCloseMC : collisionsMC) { // check for closeby collisions in MC
        int diffColl = collisionCloseMC.globalIndex() - collisionMC.globalIndex();
        if( diffColl >= mergeCollisionsDeltaMin && diffColl <= mergeCollisionsDeltaMax) { // check if n collisions prior or after
          if (collisionCloseMC.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
            continue;
          }
          if( diffColl == 0) {
            continue;
          }
          if (TMath::Abs(collisionCloseMC.posZ() - posZtrue) < splitCollisionsDeltaZPart) {
            closeByCollisionIDs.insert(collisionMC.globalIndex()); // Save the ID of the close-by collision
            break; // closeby collision in MC, don't use this event
          }
        }
      }
    }
    // now make reconstructed-level checks
    for (auto const& collision : collisions) { // loop over reconstructed collisions
      if (collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        return;
      }
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        continue;
      }
      if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
        continue;
      }
      float weight = collision.weight();
      float pTHat = collision.mcCollision().ptHard();
      const auto tracksColl = tracks.sliceBy(perCol, collision.globalIndex());

      // fill track histograms for all collisions
      for (auto const& track : tracksColl) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
        registry.fill(HIST("h_track_pt_eta"), track.pt(), track.eta(), weight);
        registry.fill(HIST("h_track_pt_phi"), track.pt(), track.phi(), weight);
        // checks on track distributions with/without collision and with/without MC particle
        if(!track.has_collision() || track.collisionId() != collision.globalIndex()) {
          registry.fill(HIST("h_track_pt_no_collision"),  track.pt(), weight);
        }
        else {
          registry.fill(HIST("h_track_pt_collision"),  track.pt(), weight);
        }
        bool has_MCparticle = track.has_mcParticle();
        if(!has_MCparticle) {
          registry.fill(HIST("h2_track_pt_pt_hat_no_particle"), track.pt(), collision.mcCollision().ptHard(), weight);
        }
        else {
          registry.fill(HIST("h2_track_pt_pt_hat_particle"), track.pt(), collision.mcCollision().ptHard(), weight);
        }

        // check outlier tracks and neighbouring collisions
        registry.fill(HIST("h2_pt_hat_track_pt"), pTHat, track.pt());
        if(track.pt() > 1.5 * pTHat) { // high weight outlier track
          registry.fill(HIST("h_track_pt_outlier"), track.pt());
          registry.fill(HIST("h2_pt_hat_track_pt_outlier"),pTHat, track.pt());
          for (auto const& collisionOutlier : collisions) { // find collisions closeby
            int diffColl = collision.globalIndex() - collisionOutlier.globalIndex();
            if(abs(diffColl) < 6) {
              float eventWeightOutlier = collisionOutlier.mcCollision().weight();
              double pTHatOutlier = collisionOutlier.mcCollision().ptHard();
              registry.fill(HIST("h2_neighbour_pt_hat_outlier"), float(diffColl+0.1), pTHatOutlier, eventWeightOutlier);
              registry.fill(HIST("h2_neighbour_track_pt_outlier"), float(diffColl+0.1), track.pt(), eventWeightOutlier);
            }
          }
        }
        // all
        for (auto const& collisionOutlier : collisions) { // find collisions closeby
          float eventWeightOutlier = collisionOutlier.mcCollision().weight();
          double pTHatOutlier = collisionOutlier.mcCollision().ptHard();
          int diffColl = collision.globalIndex() - collisionOutlier.globalIndex();

          if(abs(diffColl) < 6) {
            //LOG(info) << "pThat = " << pTHat << "pThat neighbour = "<<pTHatOutlier;
            registry.fill(HIST("h2_neighbour_pt_hat_all"), float(diffColl+0.1), pTHatOutlier, eventWeightOutlier);
            registry.fill(HIST("h2_neighbour_track_pt_all"), float(diffColl+0.1), track.pt(), eventWeightOutlier);
          }
        }
      }

      // now check for close collisions in Z and fill track histograms excluding these events
      float posZtrue = collision.mcCollision().posZ();
      float posZrec = collision.posZ();
      registry.fill(HIST("h_Z_resolution"), posZtrue - posZrec);
      registry.fill(HIST("h_Z_resolution_wide"), posZtrue - posZrec);
      int mcCollisionID = collision.mcCollisionId(); // Get the corresponding MC collision ID
      if (closeByCollisionIDs.find(mcCollisionID) != closeByCollisionIDs.end()) {
        registry.fill(HIST("h_Z_reco_rejected"), posZrec);
        registry.fill(HIST("h_Z_true_rejected"), posZtrue);
        continue; // Skip this reconstructed collision as it corresponds to a close-by MC collision
      }
      registry.fill(HIST("h_Z_reco_accepted"), posZrec);
      registry.fill(HIST("h_Z_true_accepted"), posZtrue);

      // fill track histograms for accepted collisions
      for (auto const& track : tracksColl) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
        registry.fill(HIST("h_track_pt_accepted"), track.pt(), weight);
        registry.fill(HIST("h_track_eta_accepted"), track.eta(), weight);
        registry.fill(HIST("h_track_phi_accepted"), track.phi(), weight);
        registry.fill(HIST("h_track_pt_eta_accepted"), track.pt(), track.eta(), weight);
        registry.fill(HIST("h_track_pt_phi_accepted"), track.pt(), track.phi(), weight);
      }
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processTracksBC, "jet finder QA outliers", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetOutlierQATask>(cfgc, TaskName{"jet-outlier-qa"})}; // o2-linter: disable=name/o2-task,name/workflow-file
}
