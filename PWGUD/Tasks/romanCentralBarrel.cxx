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
/// \since  10.11.2022

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCMonteCarloCentralBarrelHelper.h"
#include "PWGUD/Core/RomanPIDhelper.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct RomanCentralBarrel {

	// Global varialbes
	int countCollisions;
	Service<TDatabasePDG> pdg;

	HistogramRegistry registry{
		"registry",
		{
			{"hEffectOfSelections", "Effect of cuts;Selection (-);Number of events (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEffectOfTrackSelections", "Effect of track cuts;Selection (-);Number of tracks (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEventSelectionTask", ";Selection (-);Number of events (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEventSelectionTaskParameters", ";Selection (-);Number of events (-)", { HistType::kTH1D, { {kNsel,-0.5,kNsel-0.5} } } },
			{"hNcontributionsToCollision", ";Number of contributions in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNtracks", ";Number of tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNtracksWide", ";Number of tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {1000,-0.5,999.5} } } },
			{"hNcontributionsToVertex", ";Number of contributors to vertex (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNprimaryTracks", ";Number of primary tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", { HistType::kTH1D, { {1200,2.8,3.4} } } },
			{"hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", { HistType::kTH1D, { {1200,1.0,5.} } } },
			{"hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,2.} } } },
			{"hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,10.} } } },
			{"hMotherPt", ";Mother #it{p_{#rm T}} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,2.} } } },
			{"hMotherPhi", ";Mother #phi (rad);Number of events (-)", { HistType::kTH1D, { {50,-TMath::Pi(),TMath::Pi()} } } },
			{"hMotherRapidity", ";Mother #y (-);Number of events (-)", { HistType::kTH1D, { {500,-2.,2.} } } },
			{"hVtxZ", ";Vertex z-position (cm);Number of events (-)", { HistType::kTH1D, { {1000,-30.,30.} } } },
			{"hVtxTransversal", ";Vertex x-position (cm);Vertex y-position (cm)", { HistType::kTH2D, { {1000,-0.1,0.1}, {1000,-0.1,0.1} } } },
			{"hTPCsignalVsMom", "All tracks;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCsignalVtxContributors", ";Vtx contributor 1 - TPC d#it{E}/d#it{x} (arb. units);Vtx contributor 2 - TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 950,10.,200. }, { 950,10.,200. } } } },
			{"hTPCsignalElectrons", ";Identified electron 1 - TPC d#it{E}/d#it{x} (arb. units);Identified electron 2 - TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 950,10.,200. }, { 950,10.,200. } } } },
			{"hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,10.} } } },
			{"hTrackPt", ";Track #it{p_{#rm T}} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,5.} } } },
			{"hTrackPhi", ";Track #phi (rad);Number of events (-)", { HistType::kTH1D, { {50,-2*TMath::Pi(),2*TMath::Pi()} } } },
			{"hTrackEta", ";Track #eta (-);Number of events (-)", { HistType::kTH1D, { {500,-7.,7.} } } },
			{"hTrackEnergy", ";Track Energy (GeV);Number of events (-)", { HistType::kTH1D, { {500,0.,100.} } } },
			{"hTrackRapidity", ";Track rapidity (-);Number of events (-)", { HistType::kTH1D, { {500,-7.,7.} } } }
		}
	};

	// declare configurables
	Configurable<bool> verboseInfo{"verboseInfo", true, {"Print general info to terminal; default it true."}};
	Configurable<bool> printMetaInfo{"printMetaInfo", false, {"Print general info to terminal about collision, tracks, particles...; default it false."}};
	Configurable<bool> verboseDebug{"verboseDebug", false, {"Print debug info to terminal; default it false."}};
	Configurable<int> applyTrackCuts{"applyTrackCuts", 0, {"Apply n selections on track; default it no cut."}};
	Configurable<int> applySingleTrackCut{"applySingleTrackCut", 0, {"Apply selection n on track, applyTrackCuts must be at maximum for full usage; default it no cut."}};

	// declare filters
//	Filter nCollisionContributorsFilter = aod::collision::numContrib > 2;

	// declare shortcuts
	using FullTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID>;
	using FullCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;

	// init
	void init(InitContext&){
		if (verboseInfo) printLargeMessage("INIT METHOD");
		countCollisions = 0;

	} // end init

	// run
	void run(ProcessingContext& context){

		if (verboseInfo) printLargeMessage("RUN METHOD");
		if (verboseDebug) LOGF(info,"countCollisions = %d",countCollisions);

	} // end run


	// process
	void processCentralBarrel(FullCollision const& collision,
														 FullTracks const& tracks,
														 soa::Join<aod::Collisions, aod::EvSels> const& infoEventSelection){

		const float P_MASS[5] = {constants::physics::MassElectron,
													 constants::physics::MassMuon,
													 constants::physics::MassPionCharged,
													 constants::physics::MassKaonCharged,
													 constants::physics::MassProton};

		if (verboseInfo) printLargeMessage("NEW COLLISION");
		countCollisions++;
		registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(0);

		auto thisCollisionEventSelectionInfo = infoEventSelection.iteratorAt(getEvSelsIndexOfThisCollisionBC(infoEventSelection,collision.globalBC()));
		fillEventSelectionHistogram(registry,thisCollisionEventSelectionInfo);
		fillTrackSelectionHistogram(registry,tracks,collision);

		registry.get<TH1>(HIST("hNcontributionsToCollision"))->Fill(collision.numContrib());

		int nTracksInThisCollision = tracks.size();
		if (verboseDebug) LOGF(info,"collision number = %d has %d contributors",countCollisions,collision.numContrib());

		registry.get<TH1>(HIST("hNtracks"))->Fill(nTracksInThisCollision);
		registry.get<TH1>(HIST("hNtracksWide"))->Fill(nTracksInThisCollision);
		registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
		registry.get<TH2>(HIST("hVtxTransversal"))->Fill(collision.posX(),collision.posY());

		int countPrimaryTracks = 0;
		int countContributorsToVertex = 0;
		int countGlobalContributorsToVertex = 0;
		for (auto& track : tracks){
			if (!selectTrack(track,applyTrackCuts)) continue;
			if (isUDprimaryTrack(track)) countContributorsToVertex++;
			if (selectTrack(track,2)) countGlobalContributorsToVertex++;
			if (isUDprimaryTrack(track)) countPrimaryTracks++;
			if (track.hasTPC()) {
				registry.get<TH2>(HIST("hTPCsignalVsMom"))->Fill(momentum(track.px(),track.py(),track.pz()),track.tpcSignal());
			}
		}
		registry.get<TH1>(HIST("hNcontributionsToVertex"))->Fill(countContributorsToVertex);
		registry.get<TH1>(HIST("hNprimaryTracks"))->Fill(countPrimaryTracks);

		//
		// Playing around
		//
		if (verboseInfo) printLargeMessage("Primary + Global + FIT");
		for (auto& track : tracks){
			if (!trackSelection(track,2)) continue;
			if (!trackSelection(track,5)) continue;
			if (isFITempty(collision)) continue;
			if (printMetaInfo) printTrackData(track);
			float trkPx = track.px();
			float trkPy = track.py();
			float trkPz = track.pz();
			registry.get<TH1>(HIST("hTrackP"))->Fill(momentum(trkPx,trkPy,trkPz));
			registry.get<TH1>(HIST("hTrackPt"))->Fill(track.pt());
			registry.get<TH1>(HIST("hTrackPhi"))->Fill(phi(trkPx,trkPy));
			registry.get<TH1>(HIST("hTrackEta"))->Fill(eta(trkPx,trkPy,trkPz));
			int typeParticle = testPIDhypothesis(track);
			if (typeParticle >= 0 && typeParticle < (sizeof(P_MASS)/sizeof(float))){
				float mass = P_MASS[typeParticle];
				registry.get<TH1>(HIST("hTrackEnergy"))->Fill(energy(mass,trkPx,trkPy,trkPz));
				registry.get<TH1>(HIST("hTrackRapidity"))->Fill(rapidity(mass,trkPx,trkPy,trkPz));
			}
		}

		//
		// fetching FT0, FDD, FV0 information
		//
		bool FITisEmpty = isFITempty(collision);
		if (FITisEmpty) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(1);
		FITisEmpty = isEvSelFITempty(thisCollisionEventSelectionInfo);
		if (FITisEmpty) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(2);
		if (countContributorsToVertex == 2) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(3);
		if (countGlobalContributorsToVertex == 2) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(4);


		//
		// TPC signal info starts here
		//
		TLorentzVector sumOfVertexTracks, trackFromVertex;
		float dedxVC[2] = {-999.,-999.};
		const int NPARTICLES = 5;//TODO think of more elegant way how to pass this info
		float dedx[2][NPARTICLES] = {{-999.,-999.,-999.,-999.,-999.},{-999.,-999.,-999.,-999.,-999.}};
		int countTrks = 0;
		if (countContributorsToVertex == 2) {
			for (auto& track : tracks){
				if (!isUDprimaryTrack(track)) continue;
				countTrks++;
				trackFromVertex.SetPxPyPzE(track.px(),track.py(),track.pz(),
																	 energy(constants::physics::MassElectron,track.px(),track.py(),track.pz()));
				sumOfVertexTracks+=trackFromVertex;
				if (track.hasTPC()){
					if (countTrks == 1) dedxVC[0] = track.tpcSignal();
					else if (countTrks == 2) dedxVC[1] = track.tpcSignal();
					else LOGP(warning,"!!!Strange behaviour in the loop, check!!!");
					int typeParticle = testPIDhypothesis(track);
					if (typeParticle >= 0 && typeParticle <= NPARTICLES){
							if (countTrks == 1) dedx[0][typeParticle] = track.tpcSignal();
							else if (countTrks == 2) dedx[1][typeParticle] = track.tpcSignal();
							else LOGP(warning,"!!!Strange behaviour in the loop, check!!!");
					}
				}
			}
			registry.get<TH2>(HIST("hTPCsignalVtxContributors"))->Fill(dedxVC[0],dedxVC[1]);
			// pairs of different particles are filled in the following histograms, but outside of the range
			registry.get<TH2>(HIST("hTPCsignalElectrons"))->Fill(dedx[0][P_ELECTRON],dedx[1][P_ELECTRON]);
			registry.get<TH2>(HIST("hTPCsignalMuons"))->Fill(dedx[0][P_MUON],dedx[1][P_MUON]);
			registry.get<TH2>(HIST("hTPCsignalPions"))->Fill(dedx[0][P_PION],dedx[1][P_PION]);
			registry.get<TH2>(HIST("hTPCsignalKaons"))->Fill(dedx[0][P_KAON],dedx[1][P_KAON]);
			registry.get<TH2>(HIST("hTPCsignalProtons"))->Fill(dedx[0][P_PROTON],dedx[1][P_PROTON]);
			registry.get<TH1>(HIST("hInvariantMass"))->Fill(sumOfVertexTracks.M());
			registry.get<TH1>(HIST("hInvariantMassWide"))->Fill(sumOfVertexTracks.M());
			registry.get<TH1>(HIST("hMotherP"))->Fill(sumOfVertexTracks.P());
			registry.get<TH1>(HIST("hMotherPwide"))->Fill(sumOfVertexTracks.P());
			registry.get<TH1>(HIST("hMotherPt"))->Fill(sumOfVertexTracks.Pt());
			registry.get<TH1>(HIST("hMotherPhi"))->Fill(sumOfVertexTracks.Phi());
			registry.get<TH1>(HIST("hMotherRapidity"))->Fill(sumOfVertexTracks.Rapidity());
		}

	} // end processCentralBarrel


	PROCESS_SWITCH(RomanCentralBarrel, processCentralBarrel, "Iterate central barrel data using UD tables", true);

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<RomanCentralBarrel>(cfgc, TaskName{"roman-central-barrel"})
	};
}