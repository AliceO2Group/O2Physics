// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lumiTask{
        HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
	int first_time = 1530314294062; // to be updated
	void init(o2::framework::InitContext&) {
                histos.add("vertexx", "", HistType::kTH1F,{{1000,-1,1,"x"}});
		histos.add("vertexy", "", HistType::kTH1F,{{1000,-1,1,"y"}});
		histos.add("timestamp", "", HistType::kTH1F,{{1000,0,5e7,"t"}});
		histos.add("vertexx_timestamp", "", HistType::kTH2F,{{1000,0,5e7,"t"},{1000,-1,1,"x"} });
		histos.add("vertexy_timestamp", "", HistType::kTH2F,{{1000,0,5e7,"t"},{1000,-1,1,"y"}});
	} // init

        void process(aod::Collision const& collision, aod::BCsWithTimestamps const&){
		auto bc = collision.bc_as<aod::BCsWithTimestamps>();
		histos.fill(HIST("vertexx"), collision.posX() ); 
	 	histos.fill(HIST("vertexy"), collision.posY() );	
		histos.fill(HIST("timestamp"), bc.timestamp()-1530314294062);
//		LOGF(info, "Got timestamp %llu", bc.timestamp());
		histos.fill(HIST("vertexx_timestamp"), bc.timestamp()-1530314294062, collision.posX() );
		histos.fill(HIST("vertexy_timestamp"), bc.timestamp()-1530314294062, collision.posY() );
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
        return WorkflowSpec{adaptAnalysisTask<lumiTask>(cfgc)};
}
