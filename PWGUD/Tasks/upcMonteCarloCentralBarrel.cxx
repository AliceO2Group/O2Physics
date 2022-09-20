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
///
/// \brief
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  12.07.2022

//#include <algorithm>
//#include <iterator>

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
//#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;


struct UPCFilterCentralBarrel {

	// Global varialbes
	bool isMonteCarloData = false;


	HistogramRegistry registry{
		"registry",
		{
			{"hEffectOfSelections", "Effect of cuts;Selection (-);Number of events (-)", { HistType::kTH1D, { {10,-0.5,9.5} } } },
		}
	};

	// declare shortcuts
	using MeasuredTCs = soa::Join<aod::Tracks, aod::TracksExtra>;

	// init
	void init(InitContext&){} // end init

	// process
	void analyseBarrelTracks(aod::Collision const& collision,
													 aod::McCollisions const& mcCollisions,
	                         MeasuredTCs& barrelTracks){

		registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(3);

		LOGF(info, "vtx-z (data) = %f", collision.posZ());
		if (mcCollisions != nullptr) LOGF(info, "vtx-z (data) = %f | size (MC) = %d", collision.posZ(), mcCollisions.size());

	} // end analyseBarrelTracks

	void processMeasuredData(aod::Collision const& collision,
	                         MeasuredTCs& barrelTracks){
		analyseBarrelTracks(collision, (aod::McCollisions)nullptr, barrelTracks);
	} // end processMeasuredData

	void processMonteCarloData(aod::Collision const& collision,
	                           aod::McCollisions const& mcCollisions,
	                           MeasuredTCs& barrelTracks){
		isMonteCarloData = true;
		analyseBarrelTracks(collision, mcCollisions, barrelTracks);
	} // end processMonteCarloData

	PROCESS_SWITCH(UPCFilterCentralBarrel, processMeasuredData, "Process tables with measured data only", false);
	PROCESS_SWITCH(UPCFilterCentralBarrel, processMonteCarloData, "Process tables with Monte Carlo data", false);

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<UPCFilterCentralBarrel>(cfgc, TaskName{"upc-filter-central-barrel"})
	};
}