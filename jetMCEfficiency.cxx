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

//jet MC particle level and detector level tracks and jets task
//
// \author Wenhui Feng 

#include <cmath>
#include <TRandom3.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

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
using namespace o2::aod::track;
using namespace o2::aod::evsel;


struct JetMCEffTask {

    HistogramRegistry registry;
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
    Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
    Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
    Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
    Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
    Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};

    int eventSelection = -1;
    int trackSelection = -1;
    std::vector<double> jetRadiiValues;
    Service<o2::framework::O2DatabasePDG> pdg;



    void init(o2::framework::InitContext&)
    {
        eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
        trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
        jetRadiiValues = (std::vector<double>)jetRadii;
        auto jetRadiiBins = (std::vector<double>)jetRadii;
        if (jetRadiiBins.size() > 1) {
            jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
        } else {
            jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
        }

        if (doprocessMCPAndMCPtracks) {
            registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
            registry.add("h2_centrality_collisions", "event status vs. centrality;number of event; centrality ;entries", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
            registry.add("h2_MC_Particle_pt_eta", "MC particle pT & Eta; #it{p}_{T}; #eta", {HistType::kTH2F, {{200, 0., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_MC_Particle_pt_eta_phi", "MC particle pt vs eta vs phi;pT (GeV/c); #eta; #varphi", {HistType::kTH3F, {{200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
            registry.add("h_MCParticle_counts", "MC particle couts after each cuts; counts", {HistType::kTH1F, {{5, 0.0, 5.0}}});
            registry.add("h_mcd_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
            registry.add("h2_centrality_mcd_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
            registry.add("h2_track_pt_eta", "track pT & eta; #it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{200, 0., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_track_pt_eta_phi", "MC detector level track pt vs eta vs phi;pT (GeV/c); #eta; #varphi", {HistType::kTH3F, {{200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
            registry.add("h_MCD_tracks_counts", "MC detector level track couts after each cuts; counts", {HistType::kTH1F, {{5, 0.0, 5.0}}});
        }
    
        if (doprocessMCRho) {
            registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
            registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
            registry.add("h2_ntracks_rhom", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}});
            registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
            registry.add("h2_centrality_rhom", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}});
        }
        if (doprocessJetsMCD) {
            registry.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
            registry.add("h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
            registry.add("h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
            registry.add("h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
            registry.add("h2_centrality_jet_pt", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
            registry.add("h2_centrality_jet_eta", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
            registry.add("h2_centrality_jet_phi", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
            registry.add("h2_centrality_jet_ntracks", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
            registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
            registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
        }
        if (doprocessJetsMCP) {
            registry.add("h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
            registry.add("h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
            registry.add("h_jet_phi_part", "jet #varphi;#varphi_{jet}^{part};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
            registry.add("h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
            registry.add("h3_jet_r_part_jet_pt_part_jet_eta_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_jet_r_part_jet_pt_part_jet_phi_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
            registry.add("h3_jet_r_part_jet_eta_part_jet_phi_part", ";#it{R}_{jet}^{part};#eta_{jet}^{part};#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
            registry.add("h3_jet_r_part_jet_pt_part_jet_ntracks_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
            registry.add("h3_jet_r_part_jet_pt_part_track_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet tracks}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
            registry.add("h3_jet_r_part_jet_pt_part_track_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_jet_r_part_jet_pt_part_track_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
            // registry.add("h_jet_phat_part", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0, 350}}});
            // registry.add("h_jet_ptcut_part", "p_{T} cut;p_{T,jet}^{part} (GeV/#it{c});N;entries", {HistType::kTH2F, {{200, 0, 200}, {20, 0, 5}}});
            // registry.add("h_jet_phat_part_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0, 350}}});
        }
        if (doprocessJetsRhoAreaSubMCD) {

            registry.add("h_jet_pt_mcd_rhoareasubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{400, -200., 200.}}});
            registry.add("h_jet_eta_mcd_rhoareasubtracted", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
            registry.add("h_jet_phi_mcd_rhoareasubtracted", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
            registry.add("h_jet_ntracks_mcd_rhoareasubtracted", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
            registry.add("h2_centrality_jet_pt_mcd_rhoareasubtracted", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {400, -200., 200.}}});
            registry.add("h2_centrality_jet_eta_mcd_rhoareasubtracted", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
            registry.add("h2_centrality_jet_phi_mcd_rhoareasubtracted", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
            registry.add("h2_centrality_jet_ntracks_mcd_rhoareasubtracted", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
            registry.add("h3_jet_r_jet_pt_centrality_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {1200, -10.0, 110.0}}});
            registry.add("h3_jet_r_jet_pt_jet_eta_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_jet_r_jet_eta_jet_phi_mcd_rhoareasubtracted", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
            registry.add("h3_jet_r_jet_pt_jet_ntracks_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {200, -0.5, 199.5}}});
            registry.add("h3_jet_r_jet_pt_jet_area_mcd_rhoareasubtracted", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {300, 0., 3.}}});
            registry.add("h3_jet_r_jet_pt_track_pt_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {200, 0., 200.}}});
            registry.add("h3_jet_r_jet_pt_track_eta_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
            registry.add("h3_jet_r_jet_pt_track_phi_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {160, -1.0, 7.0}}});
            registry.add("h3_jet_r_jet_pt_jet_pt_mcd_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {400, -200., 200.}}});
        }

    }

    // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
    // Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax);
    // Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

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

    void processMCPAndMCPtracks(soa::Join<JetCollisions, aod::JMcCollisionLbs, aod::EvSels>::iterator const& collision,
                                soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& jmcparticles,
                                soa::Join<JetTracksMCD, aod::JTrackExtras> const& mcdtracks)
    {
        // registry.fill(HIST("h_collisions"), 0.5);
        // registry.fill(HIST("h2_centrality_collisions"), mcCollision.centrality(), 0.5);
        // if (!jetderiveddatautilities::selectCollision(mcCollision, eventSelection)) {
        //     return;
        // }
        // registry.fill(HIST("h_collisions"), 1.5);
        // registry.fill(HIST("h2_centrality_collisions"), mcCollision.centrality(), 1.5);
        // if (!(abs(mcCollision.posZ()) < vertexZCut)) {
        //     return;
        // }   
        // registry.fill(HIST("h_collisions"), 1.5);
        // registry.fill(HIST("h2_centrality_collisions"), mcCollision.centrality(), 2.5);
        registry.fill(HIST("h_mcd_collisions"), 0.5);
        registry.fill(HIST("h2_centrality_mcd_collisions"), collision.centrality(), 0.5);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
            return;
        }
        registry.fill(HIST("h_mcd_collisions"), 1.5);
        registry.fill(HIST("h2_centrality_mcd_collisions"), collision.centrality(), 1.5);
        if (!(abs(collision.posZ()) < vertexZCut)) {
            return;
        }   
        registry.fill(HIST("h_mcd_collisions"), 2.5);
        registry.fill(HIST("h2_centrality_mcd_collisions"), collision.centrality(), 2.5);
        //MC particle level tracks loop
        for (auto& jmcparticle : jmcparticles) {
            registry.fill(HIST("h_MCParticle_counts"), 0.5);
            if (!jmcparticle.isPhysicalPrimary()) {
                continue;
            }
            registry.fill(HIST("h_MCParticle_counts"), 1.5);
            auto pdgParticle = pdg->GetParticle(jmcparticle.pdgCode());
            if (pdgParticle == nullptr) {
                continue;
            }
            
            if (std::abs(pdgParticle->Charge()) < 3) {
                continue;
            }
            registry.fill(HIST("h_MCParticle_counts"), 2.5);
            if (jmcparticle.eta() < -0.9 || jmcparticle.eta() > 0.9) {
                continue;
            }
            registry.fill(HIST("h_MCParticle_counts"), 3.5);
            registry.fill(HIST("h3_MC_Particle_pt_eta_phi"), jmcparticle.pt(), jmcparticle.eta(), jmcparticle.phi());
            registry.fill(HIST("h2_MC_Particle_pt_eta"), jmcparticle.pt(), jmcparticle.eta());
           
        }

        for (auto const& track : mcdtracks) {
            registry.fill(HIST("h_MCD_tracks_counts"), 0.5);
            if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
                continue;
            }
            registry.fill(HIST("h_MCD_tracks_counts"), 1.5);
            if (!track.has_mcParticle()) {
                continue;
            }
            registry.fill(HIST("h_MCD_tracks_counts"), 2.5);
            if (track.eta() < -0.9 || track.eta() > 0.9) {
                continue;
            }
            registry.fill(HIST("h_MCD_tracks_counts"), 3.5);
            registry.fill(HIST("h3_track_pt_eta_phi"), track.pt(), track.eta(), track.phi());
            registry.fill(HIST("h2_track_pt_eta"), track.pt(), track.eta());
        }

    }
    PROCESS_SWITCH(JetMCEffTask, processMCPAndMCPtracks, "process track in MC particle level track histos", true);
    
    
    void processMCDtracks(JetCollision const& collision,
                          soa::Join<JetTracksMCD, aod::JTrackExtras> const& mcdtracks)
    {
        
        //MC detector level tracks loop
        
         
    }
    PROCESS_SWITCH(JetMCEffTask, processMCDtracks, "process track in MC detector level track histos", true);
    
    void processMCRho(soa::Join<JetCollisions, aod::BkgChargedRhos>::iterator const& mcCollision,
                      soa::Join<JetTracksMCD, aod::JTrackExtras> const& mcdtracks)
    {
        if (!jetderiveddatautilities::selectCollision(mcCollision, eventSelection)) {
            return;
        }
        if (!(abs(mcCollision.posZ()) < vertexZCut)) {
            return;
        } 
        int nTracks = 0;
        for (auto const& track : mcdtracks) {
            if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
                nTracks++;
            }
        }
        registry.fill(HIST("h2_centrality_ntracks"), mcCollision.centrality(), nTracks);
        registry.fill(HIST("h2_ntracks_rho"), nTracks, mcCollision.rho());
        registry.fill(HIST("h2_ntracks_rhom"), nTracks, mcCollision.rhoM());
        registry.fill(HIST("h2_centrality_rho"), mcCollision.centrality(), mcCollision.rho());
        registry.fill(HIST("h2_centrality_rhom"), mcCollision.centrality(), mcCollision.rhoM());
           
        }
    PROCESS_SWITCH(JetMCEffTask, processMCRho, "process MC collisions rho histos", true);
    
    void processJetsMCD(JetCollision const& mccollisions, 
                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, 
                        JetTracks const&)
    {
        if (!jetderiveddatautilities::selectCollision(mccollisions, eventSelection)) {
            return;
        }
        if (!(abs(mccollisions.posZ()) < vertexZCut)) {
            return;
        } 
        for (auto const& jet : jets) {
            if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
                continue;
            }
            if (!isAcceptedJet<JetTracks>(jet)) {
                continue;
            }
            if (jet.r() == round(selectedJetsRadius * 100.0f)) {
                registry.fill(HIST("h_jet_pt"), jet.pt());
                registry.fill(HIST("h_jet_eta"), jet.eta());
                registry.fill(HIST("h_jet_phi"), jet.phi());
                registry.fill(HIST("h_jet_ntracks"), jet.tracksIds().size());
                registry.fill(HIST("h2_centrality_jet_pt"), mccollisions.centrality(), jet.pt());
                registry.fill(HIST("h2_centrality_jet_eta"), mccollisions.centrality(), jet.eta());
                registry.fill(HIST("h2_centrality_jet_phi"), mccollisions.centrality(), jet.phi());
                registry.fill(HIST("h2_centrality_jet_ntracks"), mccollisions.centrality(), jet.tracksIds().size());
            }
            for (auto& constituent : jet.template tracks_as<JetTracks>()) {
                registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), constituent.pt());
                registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), constituent.eta());
                registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), constituent.phi());
            }
        }
        
    }
    PROCESS_SWITCH(JetMCEffTask, processJetsMCD, "jet finder QA mcd", true);
    
    void processJetsRhoAreaSubMCD(soa::Join<JetCollisions, aod::BkgChargedRhos>::iterator const& mccollision, 
                                  soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, 
                                  JetTracks const&)
    {
        if (!jetderiveddatautilities::selectCollision(mccollision, eventSelection)) {
            return;
        }
        if (!(abs(mccollision.posZ()) < vertexZCut)) {
            return;
        } 
        for (auto const& jet : jets) {
            if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
                continue;
            }
            if (!isAcceptedJet<JetTracks>(jet)) {
                continue;
            }
            if (jet.r() == round(selectedJetsRadius * 100.0f)) {
                registry.fill(HIST("h_jet_pt_mcd_rhoareasubtracted"), jet.pt() - (mccollision.rho() * jet.area()));
                registry.fill(HIST("h_jet_eta_mcd_rhoareasubtracted"), jet.eta());
                registry.fill(HIST("h_jet_phi_mcd_rhoareasubtracted"), jet.phi());
                registry.fill(HIST("h_jet_ntracks_mcd_rhoareasubtracted"), jet.tracksIds().size());
                registry.fill(HIST("h2_centrality_jet_pt_mcd_rhoareasubtracted"), mccollision.centrality(), jet.pt() - (mccollision.rho() * jet.area()));
                registry.fill(HIST("h2_centrality_jet_eta_mcd_rhoareasubtracted"), mccollision.centrality(), jet.eta());
                registry.fill(HIST("h2_centrality_jet_phi_mcd_rhoareasubtracted"), mccollision.centrality(), jet.phi());
                registry.fill(HIST("h2_centrality_jet_ntracks_mcd_rhoareasubtracted"), mccollision.centrality(), jet.tracksIds().size());
            }

            registry.fill(HIST("h3_jet_r_jet_pt_centrality_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho()  * jet.area()), mccollision.centrality());
            registry.fill(HIST("h3_jet_r_jet_pt_jet_eta_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho()  * jet.area()), jet.eta());
            registry.fill(HIST("h3_jet_r_jet_eta_jet_phi_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.eta(), jet.phi());
            registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho()  * jet.area()), jet.tracksIds().size());
            registry.fill(HIST("h3_jet_r_jet_pt_jet_area_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho()  * jet.area()), jet.area());
            registry.fill(HIST("h3_jet_r_jet_pt_jet_pt_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt(), jet.pt() - (mccollision.rho()  * jet.area()));

            for (auto& constituent : jet.template tracks_as<JetTracks>()) {

                registry.fill(HIST("h3_jet_r_jet_pt_track_pt_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho() * jet.area()), constituent.pt());
                registry.fill(HIST("h3_jet_r_jet_pt_track_eta_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho() * jet.area()), constituent.eta());
                registry.fill(HIST("h3_jet_r_jet_pt_track_phi_mcd_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (mccollision.rho() * jet.area()), constituent.phi());
            }
        }
        
    }
    PROCESS_SWITCH(JetMCEffTask, processJetsRhoAreaSubMCD, "jet mcd rho area-base subtraction", true);

    void processJetsMCP(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet, 
                        JetParticles const&)
    {
        // for (auto const& jet : jets) {
            if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
                return;
            }
            if (!isAcceptedJet<JetParticles>(jet)) {
                return;
            }
            if (jet.r() == round(selectedJetsRadius * 100.0f)) {
                registry.fill(HIST("h_jet_pt_part"), jet.pt());
                registry.fill(HIST("h_jet_eta_part"), jet.eta());
                registry.fill(HIST("h_jet_phi_part"), jet.phi());
                registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size());
            }

            registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta());
            registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi());
            registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi());
            registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size());

            // for (auto& constituent : jet.template tracks_as<JetParticles>()) {
            //     registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt());
            //     registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta());
            //     registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi());
            // }
        // }
    }
    PROCESS_SWITCH(JetMCEffTask, processJetsMCP, "jet finder QA mcp", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetMCEffTask>(cfgc, TaskName{"jet-mcp-efficiency"})}; }
