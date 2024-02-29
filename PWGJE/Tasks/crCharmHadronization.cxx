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

// charm lambda/d0 ratio hadronization task
//
/// \author Christian Reckziegel <christian.reckziegel@cern.ch>
//

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// creating table for storing distance data
namespace o2::aod {
  namespace DistanceSpace {
    DECLARE_SOA_COLUMN(JetPt, jetpt, float);
    DECLARE_SOA_COLUMN(JetEta, jeteta, float);
    DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
    DECLARE_SOA_COLUMN(HfMass, hfmass, float);
    DECLARE_SOA_COLUMN(HfEta, hfeta, float);
    DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
  } // end namespace DistanceTable
  DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE", 
                    DistanceSpace::JetPt,
                    DistanceSpace::JetEta,
                    DistanceSpace::JetPhi,
                    DistanceSpace::HfMass,
                    DistanceSpace::HfEta,
                    DistanceSpace::HfPhi);
} // end o2::aod namespace


// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

struct crCharmHadronizationTask {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;
  
  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // histogram configurables
  Configurable<int> distBins{"distBins", 1000, "number of bins in distance histogram"};

  void init(InitContext const&)
  {
    // define axes you want to use

    // create histograms
    //registry.add("h_track_jet_distance", ";#DeltaR_{track,jet};dN/d(#DeltaR)", {HistType::kTH1F, {{distBins, 0., 10.}}});
    registry.add("h_hf_jet_distance", ";#DeltaR_{HF,jet};dN/d(#DeltaR)", {HistType::kTH1F, {{distBins, 0., 10.}}});
  }

  void processDummy(aod::TracksIU const& tracks)
  {}
  PROCESS_SWITCH(crCharmHadronizationTask, processDummy, "Dummy process function turned on by default", false);

  void processDataChargedSubstructure(JetCollision const& collision, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets, JetTracks const& tracks, CandidatesD0Data const& hfcands) {

    double axisDistance = 0;

    for (auto& jet : jets) {
      std::cout << "This jet has pT = " << jet.pt() << " GeV/c\n";
      for (auto& jetHFCandidate : jet.hfcandidates_as<CandidatesD0Data>()) { //for jet constituents use -> auto& jetConstituent : jet.tracks_as<JetTracks>()
        /* code */
        axisDistance = sqrt(pow(jet.eta() - jetHFCandidate.eta(),2) + pow(jet.phi() - jetHFCandidate.phi(),2));
        // fill histogram
        registry.fill(HIST("h_hf_jet_distance"), axisDistance);
        // filling table
        distJetTable(jet.pt(),jet.eta(),jet.phi(),jetHFCandidate.m(),jetHFCandidate.eta(),jetHFCandidate.phi());
      }
      
    }
    
    
    
    
    
  }
  PROCESS_SWITCH(crCharmHadronizationTask, processDataChargedSubstructure, "jet substructure charged jets", true);

};

// for templates only
//using JetSubstructureD0 = crCharmHadronizationTask<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>, CandidatesD0Data, aod::JTrackD0Subs>

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<crCharmHadronizationTask>(cfgc, TaskName{"jet-charm-hadronization"})};
}
