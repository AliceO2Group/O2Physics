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

// \brief   Task for the calculation of SPC with filtered data.
// \author  Maxim Virta (maxim.virta@cern.ch), Cindy Mordasini (cindy.mordasini@cern.ch)

// Standard headers.
#include <chrono>
#include <string>
#include <vector>
#include <TRandom3.h>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers. //
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"


#include "PWGCF/Flow/Core/FlowJSPCAnalysis.h"
#include "PWGCF/Flow/Core/FlowJSPCObservables.h" 
#include "PWGCF/Flow/Core/FlowJHistManager.h"

// Namespaces and definitions.
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
								 aod::FT0sCorrected, aod::CentFT0Ms,
								 aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
								 aod::CentFDDMs, aod::CentNTPVs>;
								 
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;



struct flowJSPCAnalysisTask
{
	HistogramRegistry SPCHistograms{"SPCResults", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
	FlowJSPCAnalysis spcAnalysis;

	HistogramRegistry qaHistRegistry{"qaHistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
	FlowJHistManager histManager;

	FlowJSPCObservables SPCobservables;

	// Set Configurables here
	Configurable<bool> cfgUseNUA{"cfgUseNUA", false, "Use NUA correction"}; 
	Configurable<bool> cfgUseNUE{"cfgUseNUE", false, "Use NUE correction"}; 
	Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill QA plots"};

	Configurable<Int_t> cfgWhichSPC{"cfgWhichSPC", 0, "Which SPC observables to compute."};

	struct : ConfigurableGroup {
		Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
		Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
		Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta used for track selection."};
	} cfgTrackCuts;


	// The centrality estimators are the ones available for Run 3.
	enum centEstimators {FT0M, FT0A, FT0C, FDDM, NTPV};
	struct : ConfigurableGroup {
		Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator."};
		Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.0f, "Maximum primary vertex cut applied for the events."};
		Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum number of particles required for the event to have."};
	} cfgEventCuts;

	// Set the access to the CCDB for the NUA/NUE weights.
	struct : ConfigurableGroup {
	Configurable<bool> cfgUseCCDB{"cfgUseCCDB", true, "Use CCDB for NUA/NUE corrections."};
	Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch",
									 "Address of the CCDB to get the NUA/NUE."};
	Configurable<long> cfgTime{"ccdb-no-later-than",
								std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
								"Latest acceptable timestamp of creation for the object."};
	} cfgCCDB;
	Service<o2::ccdb::BasicCCDBManager> ccdb;


	// // Filters to be applied to the received data.
	// // The analysis assumes the data has been subjected to a QA of its selection,
	// // and thus only the final distributions of the data for analysis are saved.
	Filter collFilter = (nabs(aod::collision::posZ) < cfgEventCuts.cfgZvtxMax);
	Filter trackFilter =  (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);


	void init(InitContext const&) {
		// Add histomanager here
		spcAnalysis.SetHistRegistry(&SPCHistograms);
		spcAnalysis.CreateHistos();

		SPCobservables.SetSPCObservables(cfgWhichSPC);
		spcAnalysis.SetFullCorrSet(SPCobservables.harmonicArray);

		histManager.SetHistRegistryQA(&qaHistRegistry);
		histManager.SetDebugLog(false);
		histManager.CreateHistQA();

		/////////////////////

		// Analysis class initialization
		// blaa blaa spcAnalysis.CreateOutputObjects()
			
		ccdb->setURL(cfgCCDB.cfgURL);
		ccdb->setCaching(true);
		ccdb->setLocalObjectValidityChecking();
		ccdb->setCreatedNotAfter(cfgCCDB.cfgTime.value);
		
		// Add CCDB access here

	}

	void process(soa::Filtered<MyCollisions>::iterator const& coll, soa::Filtered<MyTracks> const& tracks) {
		if (tracks.size() < cfgEventCuts.cfgMultMin) return;
		// Int_t nTracks = tracks.size();

		float cent = -1.;
		switch (cfgEventCuts.cfgCentEst) {
			case FT0M : cent = coll.centFT0M(); break;
			case FT0A : cent = coll.centFT0A(); break;
			case FT0C : cent = coll.centFT0C(); break;
			case FDDM : cent = coll.centFDDM(); break;
			case NTPV : cent = coll.centNTPV(); break;
		}
		if (cent < 0. || cent > 70.) {return;}
		Int_t cBin = histManager.GetCentBin(cent);
		// printf("Centrality: %.2f,\tBin: %d\n", cent, cBin);
		SPCHistograms.fill(HIST("FullCentrality"), cent);
		vector<float> trackPhi;
		vector<float> trackWeight;
		int nTracks = 0;

		for (auto& track : tracks) {
			if (cfgFillQA) histManager.FillTrackQA<1>(track, cBin, coll.posZ());
			trackPhi.push_back(track.phi());
			// printf("Track phi: %.3f, %.2f\n", track.phi(), cent);
			// We get the NUE and NUA weight values from the objects fetched earlier
			// in the CCDB. These weights still need to be inverted to be used in
			// the computations of the Q-vectors.
			float wNUE = 1.;  // Unit NUE weight for pT.
			float wNUA = 1.;  // Unit NUA weight for phi.

			if (cfgUseNUE) {;}
			if (cfgUseNUA) {;}
			nTracks++;
			trackWeight.push_back(1./(wNUE*wNUA));
		}

		if (cfgFillQA) histManager.FillEventQA<1>(coll, cBin, cent, nTracks);

		spcAnalysis.CalculateQvectors(trackPhi, trackWeight);
		spcAnalysis.CalculateCorrelators(cBin);

		/* Reset the variables for the next collision. */
		// This ensures no mixing between collision can happen accidentally.
		trackPhi.clear();
		trackWeight.clear();

		LOGF(info, "Collision analysed. Next...");
	}
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{
		adaptAnalysisTask<flowJSPCAnalysisTask>(cfgc)
	};
}
