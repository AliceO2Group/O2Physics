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
/// \author Deependra Sharma

// O2
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// PWGHF
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// flaging a collision if D* meson is found.
struct HfCorrelatorDstarHadronCollisionSelector{
    Produces<aod::DmesonSelection> collisionWDstar;

    Configurable<bool> selectionFlagDstar{"selectionFlagDstar",true,"selection flag for Dstar"};
    Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
    Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

    using DstarCandidates = soa::Join<aod::HfCandDstar,aod::HfSelDstarToD0Pi>;

    SliceCache cache;

    // candidates who passed the slection criteria defined in "CandidateSelectionTables.h"
    Partition<DstarCandidates> selectedDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar; 

    

    void processCollisionSelWDstar(aod::Collision const& collision,
                                    DstarCandidates const& candidates){
        bool isDstarFound = false;
        if(selectedDstarCandidates.size() > 0){
            auto selectedDstarCandidatesGrouped = selectedDstarCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
            for(const auto & selectedCandidate: selectedDstarCandidatesGrouped){
                auto yDstar = selectedCandidate.y();
                auto pTDstar = selectedCandidate.pt();
                if(yCandMax >= 0 && yDstar > yCandMax){
                    continue;
                }
                if(ptCandMin >= 0 && pTDstar < ptCandMin){
                    continue;
                }
                isDstarFound = true;
                break;
            }
        }
        collisionWDstar(isDstarFound); // compatible with collision table (filled collision by collision)
    }
    PROCESS_SWITCH(HfCorrelatorDstarHadronCollisionSelector,processCollisionSelWDstar,"process only data for dstar hadron correlation", false);
};

struct HfCorrelatorDstarHadrons{
    Produces<aod::DstarHadronPair> rowsDstarHadronPair;

    using CollisionsWMult = soa::Join<aod::Collisions,aod::Mults,aod::DmesonSelection>;
    using DstarCandidates = soa::Join<aod::HfCandDstar,aod::HfSelDstarToD0Pi>;

    Filter CollisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;

    
    void processData(aod::BCsWithTimestamps::iterator const & bc,
                    soa::Filtered<CollisionsWMult> const & collisions, // only collisions who have altleast one D*
                    aod::TracksWDca const & tracks,
                    DstarCandidates const & candidates){
        auto timestamp = bc.timestamp();
        for(const auto & collision: collisions){

        }
        

    }

};

