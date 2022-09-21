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
#include "PWGUD/DataModel/UPCFilterCentralBarrel.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;

template <typename TCs>
int TestTPConlyPIDhypothesis(TCs barrelTrack)
// Choose, which particle it is according to PID of TPC only
{

	float pidTPC[5];
	pidTPC[P_ELECTRON] = TMath::Abs(barrelTrack.tpcNSigmaEl());
	pidTPC[P_MUON] = TMath::Abs(barrelTrack.tpcNSigmaMu());
	pidTPC[P_PION] = TMath::Abs(barrelTrack.tpcNSigmaPi());
	pidTPC[P_KAON] = TMath::Abs(barrelTrack.tpcNSigmaKa());
	pidTPC[P_PROTON] = TMath::Abs(barrelTrack.tpcNSigmaPr());

	return std::distance(std::begin(pidTPC), std::min_element(std::begin(pidTPC), std::end(pidTPC)));

}


template <typename TCs>
int TestPIDhypothesis(TCs barrelTrack)
// Choose, which particle it is according to PID of TPC and TOF
{

	float pidTPC[5], pidTOF[5];
	pidTPC[P_ELECTRON] = TMath::Abs(barrelTrack.tpcNSigmaEl());
	pidTPC[P_MUON] = TMath::Abs(barrelTrack.tpcNSigmaMu());
	pidTPC[P_PION] = TMath::Abs(barrelTrack.tpcNSigmaPi());
	pidTPC[P_KAON] = TMath::Abs(barrelTrack.tpcNSigmaKa());
	pidTPC[P_PROTON] = TMath::Abs(barrelTrack.tpcNSigmaPr());
	pidTOF[P_ELECTRON] = TMath::Abs(barrelTrack.tofNSigmaEl());
	pidTOF[P_MUON] = TMath::Abs(barrelTrack.tofNSigmaMu());
	pidTOF[P_PION] = TMath::Abs(barrelTrack.tofNSigmaPi());
	pidTOF[P_KAON] = TMath::Abs(barrelTrack.tofNSigmaKa());
	pidTOF[P_PROTON] = TMath::Abs(barrelTrack.tofNSigmaPr());

	int _id_TPC = std::distance(std::begin(pidTPC), std::min_element(std::begin(pidTPC), std::end(pidTPC)));
	float pidTPCpick[3] = {pidTPC[0], pidTPC[1], pidTPC[2]};
	int _id_TPC_pick = std::distance(std::begin(pidTPCpick), std::min_element(std::begin(pidTPCpick), std::end(pidTPCpick)));

	int _id_TOF = std::distance(std::begin(pidTOF), std::min_element(std::begin(pidTOF), std::end(pidTOF)));
//
	if (pidTOF[_id_TOF] < 999.){
		if      (_id_TPC == P_ELECTRON || _id_TPC == P_MUON || _id_TPC == P_PION){
			if      (_id_TOF == P_KAON  ) return P_KAON; // probably kaon
			else if (_id_TOF == P_PROTON) return P_PROTON; // probably proton
			else {
				if      (_id_TPC == P_ELECTRON) return P_ELECTRON; // probably electron
				else if (_id_TPC == P_MUON)     return P_MUON; // probably muon
				else                            return P_PION; // probably pion
			}
		}
		else if (_id_TPC == P_KAON){
			if      (_id_TOF == P_KAON  ) return P_KAON; // probably kaon
			else if (_id_TOF == P_PROTON) return P_PROTON; // probably proton
			else {
				if      (_id_TPC_pick == P_ELECTRON) return P_ELECTRON; // probably misidentified electron
				else if (_id_TPC_pick == P_MUON)     return P_MUON; // probably misidentified muon
				else if (_id_TPC_pick == P_PION)     return P_PION; // probably misidentified pion
			}
		}
		else {
			if      (_id_TOF == P_KAON  ) return P_KAON; // probably kaon
			else if (_id_TOF == P_PROTON) return P_PROTON; // probably proton
			else {
				if      (_id_TPC_pick == P_ELECTRON) return P_ELECTRON; // probably misidentified electron
				else if (_id_TPC_pick == P_MUON)     return P_MUON; // probably misidentified muon
				else if (_id_TPC_pick == P_PION)     return P_PION; // probably misidentified pion
			}
		}
	}
	else return _id_TPC;

	return -1;
}


struct UPCFilterCentralBarrel {

	// Global varialbes
	bool isMonteCarloData = false;


	HistogramRegistry registry{
		"registry",
		{
			{"hEffectOfSelections", "Effect of cuts;Selection (-);Number of events (-)", { HistType::kTH1D, { {10,-0.5,9.5} } } },
			{"hTPCsignalVsMom", "All tracks;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCelectronIdentified", "Tracks identified as electrons;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCmuonIdentified", "Tracks identified as muons;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCpionIdentified", "Tracks identified as pions;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCkaonIdentified", "Tracks identified as kaons;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCprotonIdentified", "Tracks identified as protons;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } }
		}
	};

	// declare shortcuts
	using MeasuredTCs = soa::Join<aod::Tracks, aod::TracksExtra,
		aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
		aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

	// init
	void init(InitContext&){} // end init

	// declare production of table
	Produces<aod::UPCTrackCandidates> selectedUPCcandidateTracks;

	// declate configurables
	Configurable<int> maxNSigmaTPC{"maxNSigmaTPC", 4, {"Maximum allowed TPC PID sigma for the most probable particle; default it 4."}};
	Configurable<int> maxNSigmaTOF{"maxNSigmaTOF", 4, {"Maximum allowed TOF PID sigma for the most probable particle; default it 4."}};

	// process
	void analyseBarrelTracks(aod::Collision const& collision,
													 aod::McCollisions const& mcCollisions,
	                         MeasuredTCs& barrelTracks)
		{

//		auto & mcCollision: mcCollisions::iterator;

		registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(3);

		LOGF(info, "vtx-z (data) = %f", collision.posZ());
		//if (mcCollisions != nullptr) LOGF(info, "vtx-z (data) = %f | size (MC) = %d", collision.posZ(), mcCollisions.size());

		for (auto & barrelTrack: barrelTracks) {

			registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(0);

			bool isTwoTracks = false;
			bool isFourTracks = false;
			int whatTrack = 0;
			bool reachedTOF = barrelTrack.hasTOF();

			// Selection criteria
			if (barrelTrack.p() < 0.1 ) continue;
			registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(1);

			if (collision.numContrib() == 2) isTwoTracks = true;
			if (collision.numContrib() == 4) isFourTracks = true;

			if (isTwoTracks && reachedTOF){
				registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(2);
				whatTrack = TestTPConlyPIDhypothesis(barrelTrack);
//				whatTrack = TestPIDhypothesis(barrelTrack);
				LOGP(info,"momentum={}",barrelTrack.p());
				LOGP(info,"TPC");
				LOGP(info,"El={}, Mu={}, Pi={}, Ka={}, Pr={}, ",barrelTrack.tpcNSigmaEl(), barrelTrack.tpcNSigmaMu(), barrelTrack.tpcNSigmaPi(), barrelTrack.tpcNSigmaKa(), barrelTrack.tpcNSigmaPr());
				LOGP(info,"tpcExpSigmaEl={}  tpcExpSignalDiffEl={}, tpcSignal={}",barrelTrack.tpcExpSigmaEl(),barrelTrack.tpcExpSignalDiffEl(),barrelTrack.tpcSignal());
				LOGP(info,"TOF");
				LOGP(info,"El={}, Mu={}, Pi={}, Ka={}, Pr={}, ",barrelTrack.tofNSigmaEl(), barrelTrack.tofNSigmaMu(), barrelTrack.tofNSigmaPi(), barrelTrack.tofNSigmaKa(), barrelTrack.tofNSigmaPr());
				LOGP(info,"tofExpSigmaEl={}  tofExpSignalEl={}  tofExpSignalDiffEl={}",barrelTrack.tofExpSigmaEl(),barrelTrack.tofExpSignalDiffEl(),barrelTrack.tofExpSignalDiffEl());
				LOGP(info,"************* this tracks is ={}", whatTrack);

				registry.get<TH2>(HIST("hTPCsignalVsMom"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());
				if (whatTrack == 0) registry.get<TH2>(HIST("hTPCelectronIdentified"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());
				if (whatTrack == 1) registry.get<TH2>(HIST("hTPCmuonIdentified"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());
				if (whatTrack == 2) registry.get<TH2>(HIST("hTPCpionIdentified"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());
				if (whatTrack == 3) registry.get<TH2>(HIST("hTPCkaonIdentified"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());
				if (whatTrack == 4) registry.get<TH2>(HIST("hTPCprotonIdentified"))->Fill(barrelTrack.p(),barrelTrack.tpcSignal());

			}



			selectedUPCcandidateTracks(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(),
																 barrelTrack.pt(), barrelTrack.p(),
																 isTwoTracks, isFourTracks, whatTrack, reachedTOF);

		} // end loop over barrel tracks

	} // end data analysis

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