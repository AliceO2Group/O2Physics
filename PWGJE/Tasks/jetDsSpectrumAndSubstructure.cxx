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

// Jet substructure and spectrum task for D_s mesons
//
// This task is used to reconstruct and analyse jets containing charged D_s
// mesons
//
/// \author Monalisa Melo <monalisa.melo@cern.ch>
//


#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGHF/Core/DecayChannels.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include <Framework/AnalysisTask.h>
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>



#include "TVector3.h"
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jet_distance
{
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);
} // namespace jet_distance

DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                    jet_distance::JetHfDist,
                    jet_distance::JetPt,
                    jet_distance::JetEta,
                    jet_distance::JetPhi,
                    jet_distance::JetNConst,
                    jet_distance::HfPt,
                    jet_distance::HfEta,
                    jet_distance::HfPhi,
                    jet_distance::HfMass,
                    jet_distance::HfY,
                    jet_distance::HfMlScore0,
                    jet_distance::HfMlScore1,
                    jet_distance::HfMlScore2);
}

struct JetDsSpecSubs {
    HistogramRegistry registry{"registry",
                            {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}},
                            {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                            {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                            {"h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                            {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                            {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                            {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                            {"h_collision_counter", "# of collisions;", {HistType::kTH1F, {{200, 0., 200.}}}},
                            {"h_jet_counter", ";# of D_{S} jets;", {HistType::kTH1F, {{6, 0., 3.0}}}},
                            {"h_ds_jet_projection", ";z^{D_{S},jet}_{||};dN/dz^{D_{S},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}}},
                            {"h_ds_jet_distance_vs_projection", ";#DeltaR_{D_{S},jet};z^{D_{S},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}}},
                            {"h_ds_jet_distance", ";#DeltaR_{D_{S},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}}},
                            {"h_ds_jet_pt", ";p_{T,D_{S} jet};dN/dp_{T,D_{S} jet}", {HistType::kTH1F, {{200, 0., 10.}}}},
                            {"h_ds_jet_eta", ";#eta_{T,D_{S} jet};dN/d#eta_{D_{S} jet}", {HistType::kTH1F, {{250, -5., 5.}}}},
                            {"h_ds_jet_phi", ";#phi_{T,D_{S} jet};dN/d#phi_{D_{S} jet}", {HistType::kTH1F, {{250, -10., 10.}}}},
                            {"h_ds_mass", ";m_{D_{S}} (GeV/c^{2});dN/dm_{D_{S}}", {HistType::kTH1F, {{1000, 0., 10.}}}},
                            {"h_ds_eta", ";#eta_{D_{S}} (GeV/c^{2});dN/d#eta_{D_{S}}", {HistType::kTH1F, {{250, -5., 5.}}}},
                            {"h_ds_phi", ";#phi_{D_{S}} (GeV/c^{2});dN/d#phi_{D_{S}}", {HistType::kTH1F, {{250, -10., 10.}}}}}
                        };
    
    Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

    Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
    Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

    Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
    Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
    
    std::vector<int> eventSelectionBits;
    int trackSelection = -1;

    Produces<aod::JetDistanceTable> distJetTable;
    
    
    void init(o2::framework::InitContext&)
    {
        eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
        trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    }


    Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
    Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;


    void processCollisions(aod::JetCollision const& collision, aod::JetTracks const& tracks)
    {

        registry.fill(HIST("h_collisions"), 0.5);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        return;
        }
        registry.fill(HIST("h_collisions"), 1.5);
        for (auto const& track : tracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
            continue;
        }
        registry.fill(HIST("h_track_pt"), track.pt());
        registry.fill(HIST("h_track_eta"), track.eta());
        registry.fill(HIST("h_track_phi"), track.phi());
        }
    }
    PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "process JE collisions", false);

    void processDataCharged(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets)
    {
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        return;
        }
        //jets -> charged jet 
        for (auto& jet : jets) {
        registry.fill(HIST("h_jet_pt"), jet.pt());
        registry.fill(HIST("h_jet_eta"), jet.eta());
        registry.fill(HIST("h_jet_phi"), jet.phi());
        }
    }
    PROCESS_SWITCH(JetDsSpecSubs, processDataCharged, "charged jets in data", false);

    void processDataChargedSubstructure(aod::JetCollision const& collision,
                                        soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents> const& jets,
                                        aod::CandidatesDsData const&,
                                        aod::JetTracks const&)
    {
        // apply event selection and fill histograms for sanity check
        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
        return;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);
        // jets -> charged jet with Ds
        for (const auto& jet : jets) {
            //number of charged jets with Ds
            registry.fill(HIST("h_jet_counter"), 0.5);
            // obtaining jet 3-vector
            TVector3 jetVector(jet.px(), jet.py(), jet.pz());

            for (const auto& dsCandidate : jet.candidates_as<aod::CandidatesDsData>()) {

                // obtaining jet 3-vector
                TVector3 dsVector(dsCandidate.px(), dsCandidate.py(), dsCandidate.pz());

                // calculating fraction of the jet momentum carried by the Ds along the direction of the jet axis
                double zParallel = (jetVector * dsVector) / (jetVector * jetVector);

                // calculating angular distance in eta-phi plane
                double axisDistance = jetutilities::deltaR(jet, dsCandidate);

                // filling histograms
                registry.fill(HIST("h_ds_jet_projection"), zParallel);
                registry.fill(HIST("h_ds_jet_distance_vs_projection"), axisDistance, zParallel);
                registry.fill(HIST("h_ds_jet_distance"), axisDistance);
                registry.fill(HIST("h_ds_jet_pt"), jet.pt());
                registry.fill(HIST("h_ds_jet_eta"), jet.eta());
                registry.fill(HIST("h_ds_jet_phi"), jet.phi());
                registry.fill(HIST("h_ds_mass"), dsCandidate.m());
                registry.fill(HIST("h_ds_eta"), dsCandidate.eta());
                registry.fill(HIST("h_ds_phi"), dsCandidate.phi());

                // filling table
                distJetTable(axisDistance,
                            jet.pt(), jet.eta(), jet.phi(), jet.tracks_as<aod::JetTracks>().size(),
                            dsCandidate.pt(), dsCandidate.eta(), dsCandidate.phi(), dsCandidate.m(), dsCandidate.y(), dsCandidate.mlScores()[0], dsCandidate.mlScores()[1], dsCandidate.mlScores()[2]);

                break; // get out of candidates' loop after first HF particle is found in jet
            } // end of DS candidates loop

        } // end of jets loop

    } // end of process function
    PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "charged HF jet substructure", false);



};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc, TaskName{"jet-ds-spectrum-subs"})}; }