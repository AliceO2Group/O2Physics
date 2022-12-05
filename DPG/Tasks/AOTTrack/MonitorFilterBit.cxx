// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//


#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;

struct CheckFilterBit {

  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};

  

  HistogramRegistry histos;
 

  void init(InitContext const&){
    
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/c)"};
    histos.add("Tracks/Reco/histpt", "pt", kTH1D, {axisPt});
  }

  //  COLUMN(TrackCutFlagLoose, trackCutFlagLoose, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged
  //DECLARE_SOA_COLUMN(TrackCutFlagPrimaryNoHF, trackCutFlagPrimaryNoHF, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged
  //DECLARE_SOA_COLUMN(TrackCutFlagHighPt, trackCutFlagHighPt, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged

  
  void process(soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension> const& tracks){// test whether the name can be changed using the PROCESS_SWITCH method (see qaEventTrack class)
    
    for(auto& track : tracks){
      if(track.isGlobalTrack())Printf("track is global track: %d",track.isGlobalTrack());
      else Printf("track is NO global track: %d",track.isGlobalTrack());
      if(track.trackCutFlagFb1()){
	Printf("FB1 is passed");
      }
      else Printf("FB1 not passed");
      
      histos.fill(HIST("Tracks/Reco/histpt"), track.pt());
      if(track.passedTPCRefit()){//o2::aod::track2<TrackSelectionExtension::passedTPCCrossedRows>()){
       	Printf("second table passed tpc refit ok");
      }
      if(track.passedDCAzFB1()){
	Printf("passed dcazFB1 ok");
      }
    }
    
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CheckFilterBit>(cfgc)
    
  };
}
