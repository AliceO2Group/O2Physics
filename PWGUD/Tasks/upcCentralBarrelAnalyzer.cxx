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

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
#include "PWGUD/DataModel/UPCFilterCentralBarrel.h"

using namespace o2;
using namespace o2::framework;

struct UPCCentralBarrelAnalyzer {

	HistogramRegistry registry{
		"registry",
		{
			{"hRegNtracks", ";n tracks; Entries", { HistType::kTH1F, { { 201,-0.5,200.5 } } } },
			{"hRegChannels", ";channel; Entries", { HistType::kTH1I, { { 4,-0.5,3.5 } } } },
			{"hRegPvsPt", ";Track #it{p} (GeV/#it{c}); Track #it{p}_{T} (GeV/#it{c})", { HistType::kTH2F, { { 2000,0.,2. }, { 2000,0.,2. } } } },
			{"hRegPvsPt2trks", ";Track #it{p} (GeV/#it{c}); Track #it{p}_{T} (GeV/#it{c})", { HistType::kTH2F, { { 2000,0.,2. }, { 2000,0.,2. } } } },
			{"hRegPvsPt4trks", ";Track #it{p} (GeV/#it{c}); Track #it{p}_{T} (GeV/#it{c})", { HistType::kTH2F, { { 2000,0.,2. }, { 2000,0.,2. } } } }
		}
	};

	void init(InitContext&){}

	void process(aod::UPCTrackCandidates const& upcTracks) {

		for (auto & upcTrack: upcTracks) {

			registry.get<TH1>(HIST("hRegNtracks"))->Fill(upcTrack.numContrib());
			registry.get<TH2>(HIST("hRegPvsPt"))->Fill(upcTrack.p(),upcTrack.pt());
			if (upcTrack.isTwoTracks()) {
				registry.get<TH1>(HIST("hRegChannels"))->Fill(1);
				registry.get<TH2>(HIST("hRegPvsPt2trks"))->Fill(upcTrack.p(),upcTrack.pt());
			}
			else if (upcTrack.isFourTracks()) {
				registry.get<TH1>(HIST("hRegChannels"))->Fill(2);
				registry.get<TH2>(HIST("hRegPvsPt4trks"))->Fill(upcTrack.p(),upcTrack.pt());
			}
			else if (upcTrack.isTwoTracks() && upcTrack.isFourTracks()) registry.get<TH1>(HIST("hRegChannels"))->Fill(3);
			else registry.get<TH1>(HIST("hRegChannels"))->Fill(0);

		} // end loop over barrel tracks


	} // end process

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<UPCCentralBarrelAnalyzer>(cfgc, TaskName{"upc-central-barrel-analyzer"})
	};
}