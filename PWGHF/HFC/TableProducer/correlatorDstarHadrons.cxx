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

/// \file correlatorDstarHadron.cxx
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

// O2
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// O2Physics
#include "Common/DataModel/Multiplicity.h"

// PWGHF
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// // flaging a collision if D* meson is found.
// struct HfCollisionSelector{
//     Produces<aod::DmesonSelection> collisionWDstar;
//     Configurable<bool> selectionFlagDstar{"selectionFlagDstar",true,"selection flag for Dstar"};
//     Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
//     Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
//     using DstarCandidates = soa::Join<aod::HfCandDstar,aod::HfSelDstarToD0Pi>;
//     SliceCache cache;
//     // candidates who passed the slection criteria defined in "CandidateSelectionTables.h"
//     Partition<DstarCandidates> selectedDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar; 
//     void processCollisionSelWDstar(aod::Collision const& collision,
//                                     DstarCandidates const& candidates){
//         bool isDstarFound = false;
//         if(selectedDstarCandidates.size() > 0){
//             auto selectedDstarCandidatesGrouped = selectedDstarCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
//             for(const auto & selectedCandidate: selectedDstarCandidatesGrouped){
//                 auto yDstar = selectedCandidate.y(constants::physics::MassDStar);
//                 auto pTDstar = selectedCandidate.pt();
//                 if(yCandMax >= 0 && yDstar > yCandMax){
//                     continue;
//                 }
//                 if(ptCandMin >= 0 && pTDstar < ptCandMin){
//                     continue;
//                 }
//                 isDstarFound = true;
//                 break;
//             }
//         }
// LOG(info)<<"processCollisionSelWDstar: isDstarFound = "<< isDstarFound;
//         collisionWDstar(isDstarFound); // compatible with collision table (filled collision by collision)
//     }
//     PROCESS_SWITCH(HfCollisionSelector,processCollisionSelWDstar,"process only data for dstar hadron correlation", true);
// };

struct HfCorrelatorDstarHadrons{
    Produces<aod::DstarHadronPair> rowsDstarHadronPair;
    Produces<aod::DmesonSelection> collisionWDstar;

    // Configurable<bool> selectOnlyCollisionWDstar{"selectOnlyCollisionWDstar",true," select on collisions which have atleast a Dstar candidate"};

    // Dstar candidate related configurable
    Configurable<bool> selectionFlagDstar{"selectionFlagDstar",true,"selection flag for Dstar"};
    Configurable<float> pTMinDstar{"pTMinDstar",1.5,"min pT of dstar candidate"};
    Configurable<float> pTMaxDstar{"pTMaxDstar",50,"max pT of dstar Candidate"};
    // Configurable<float> etaAbsMaxDstar{"etaAbsMaxDstar",1.0,"max Abs(eta) cut on Dstar candidate"};
    Configurable<float> yMaxDstar{"yMaxDstar", 0.8, "max. cand. rapidity"};
    // track related configurable
    Configurable<float> etaAbsMaxAssoTrack{"etaAbsMaxAssoTrack",1.0,"max Abs(eta) cut on Associated Track"};
    Configurable<float> dcaxyMinAssoTrack{"dcaxyMinAssoTrack",0.0,"min DCAxy of Associated Track"};
    Configurable<float> dcaxyMaxAssoTrack{"dcaxyMaxAssoTrack",0.003,"max DCAxy of Associated Track"};
    Configurable<float> dcazMinAssoTrack{"dcazMinAssoTrack",0.0, "min DCAz of Associated Track"};
    Configurable<float> dcazMaxAssoTrack{"dcazMaxAssoTrack",0.003,"max DCAz of Associated Track"};
    Configurable<float> pTMinAssoTrack{"pTMinAssoTrack",0.5,"min Pt of Associated Track"};
    Configurable<float> pTMaxAssoTrack{"pTMaxAssoTrack",50.0,"max pT of Associated Track"};

    ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};

    ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>> binningScheme{{binsZVtx, binsMultiplicity},true};
    // ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>> binningScheme;

    using CollisionsWMult = soa::Join<aod::Collisions, aod::Mults/*, aod::DmesonSelection*/>;
    using DstarCandidates = soa::Join<aod::HfCandDstar, aod::HfSelDstarToD0Pi>; // if we add two extra columns of prong0Id and prong1d then we need not to join HfD0FromDstar here.

    // Filter CollisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == selectOnlyCollisionWDstar;

    SliceCache cache;
    Preslice<DstarCandidates> perColCandidates = aod::hf_cand::collisionId;
    Preslice<aod::TracksWDca> perColTracks = aod::track::collisionId;

    // cabdidate partition
    Partition<DstarCandidates> selectedCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar;

    // track partition 
    Partition<aod::TracksWDca> selectedTracks = nabs(aod::track::eta) < etaAbsMaxAssoTrack && aod::track::pt > pTMinAssoTrack && aod::track::pt < pTMaxAssoTrack &&
                                                aod::track::dcaXY > dcaxyMinAssoTrack && aod::track::dcaXY < dcaxyMaxAssoTrack &&
                                                aod::track::dcaZ > dcazMinAssoTrack && aod::track::dcaZ < dcazMaxAssoTrack;
    
    void init (InitContext&){
        binningScheme = {{binsZVtx, binsMultiplicity},true};
    }
    
    void processData(/*soa::Filtered<CollisionsWMult>*/ CollisionsWMult const & collisions, // only collisions who have altleast one D*
                    aod::TracksWDca const & tracks,
                    DstarCandidates const & candidates,
                    aod::BCsWithTimestamps const & ){

        // LOG(info)<<"process data function called. Collision Table size = "<<collisions.size();

        for(const auto & collision: collisions){
            bool isDstarFound =false;
            auto bc = collision.bc_as<aod::BCsWithTimestamps>();
            auto timestamp = bc.timestamp();
            // LOG(info)<<"timestamp = "<<timestamp;

            auto selectedCandidatesGrouped = selectedCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
            // LOG(info)<< "selected candidates grouped according to current collision";
            auto selectedTracksGrouped = selectedTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache); 

            if((selectedCandidatesGrouped.size() == 0) || (selectedTracksGrouped.size() == 0)){
                continue;
            }
            // LOG(info)<<"pair creation starts";
            // pair creation
            for(auto &[triggerParticle, assocParticle] : soa::combinations(soa::CombinationsFullIndexPolicy(selectedCandidatesGrouped,selectedTracksGrouped))){
                auto gItriggerParticle = triggerParticle.globalIndex();
                auto gIassocParticle = assocParticle.globalIndex();

                //Track rejection based on daughter index
                if((triggerParticle.prong0Id() == gIassocParticle) || (triggerParticle.prong1Id() == gIassocParticle) || (triggerParticle.prongPiId() == gIassocParticle)){
                    continue; // rejected pair if associated particle is same as any of daughter particle
                }
                // Track rejection based on eta, pt cut. if partition works on expression column: remove this from here
                if(assocParticle.eta() > etaAbsMaxAssoTrack || assocParticle.pt() < pTMinAssoTrack || assocParticle.pt() > pTMaxAssoTrack){
                    continue; // reject pair of associated particle is not within kinematic range
                }
                // Trigger Particle Rejection
                if(triggerParticle.pt() > pTMaxDstar || triggerParticle.pt() < pTMinDstar){
                    continue;
                }
                auto yDstar = triggerParticle.y(constants::physics::MassDStar);
                if(std::abs(yDstar) > yMaxDstar){
                    continue;
                }

                // if(!isDstarFound){
                //     isDstarFound =true;
                // }
                isDstarFound =true;

                auto binNumber = binningScheme.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
                if(triggerParticle.signSoftPi() > 0){
                    rowsDstarHadronPair(collision.globalIndex(),
                                    gItriggerParticle,
                                    triggerParticle.phi(),
                                    triggerParticle.eta(),
                                    triggerParticle.pt(),
                                    triggerParticle.invMassDstar(),
                                    gIassocParticle,
                                    assocParticle.phi(),
                                    assocParticle.eta(),
                                    assocParticle.pt(),
                                    timestamp,
                                    binNumber
                                    );
                }else {
                    rowsDstarHadronPair(collision.globalIndex(),
                                    gItriggerParticle,
                                    triggerParticle.phi(),
                                    triggerParticle.eta(),
                                    triggerParticle.pt(),
                                    triggerParticle.invMassAntiDstar(),
                                    gIassocParticle,
                                    assocParticle.phi(),
                                    assocParticle.eta(),
                                    assocParticle.pt(),
                                    timestamp,
                                    binNumber
                                    );
                }
                LOG(info)<<"pair table filled";
            }  // D-H pair loop 
            collisionWDstar(isDstarFound); // filled collision by collision
            LOG(info)<<"collision table of Dstar only filled";
        } // collision loop
    } //processData

    PROCESS_SWITCH(HfCorrelatorDstarHadrons,processData,"process data only", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
    return WorkflowSpec{/*adaptAnalysisTask<HfCollisionSelector>(cfgc),*/
                      adaptAnalysisTask<HfCorrelatorDstarHadrons>(cfgc)};
}