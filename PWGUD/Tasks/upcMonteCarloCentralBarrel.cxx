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
	void processSimulatorLevel(aod::Collision const& collision,
														 aod::McCollisions const& mcCollisions,
	                           MeasuredTCs& barrelTracks){

	} // end processSimulatorLevel

	void processGeneratorLevel(aod::McCollision const& mcCollision,
														 aod::McParticles const& mcParticles){

		


	} // end processGeneratorLevel

	PROCESS_SWITCH(UPCFilterCentralBarrel, processSimulatorLevel, "Iterate MC tables with reconstructed data", false);
	PROCESS_SWITCH(UPCFilterCentralBarrel, processGeneratorLevel, "Iterate MC tables with generated data", false);

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<UPCFilterCentralBarrel>(cfgc, TaskName{"upc-mc-central-barrel"})
	};
}