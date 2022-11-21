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
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGUD/Core/UPCMonteCarloCentralBarrelHelper.h"
#include "PWGUD/DataModel/UDTables.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcMCCentralBarrel {

	// Global varialbes
	bool isFirstReconstructedCollisions, isFirstGeneratedCollisions;
	int countCollisions;
	int32_t indexCollision;
	Service<TDatabasePDG> pdg;

	HistogramRegistry registry{
		"registry",
		{
			{"hEffectOfSelections", "Effect of cuts;Selection (-);Number of events (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEffectOfTrackSelections", "Effect of track cuts;Selection (-);Number of tracks (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEventSelectionTask", ";Selection (-);Number of events (-)", { HistType::kTH1D, { {20,-0.5,19.5} } } },
			{"hEventSelectionTaskParameters", ";Selection (-);Number of events (-)", { HistType::kTH1D, { {kNsel,-0.5,kNsel-0.5} } } },
			{"hNcontributionsToReconstructedCollision", ";Number of contributions in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedParticles", ";Number of particles in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedTracks", ";Number of tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedTracksWide", ";Number of tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {1000,-0.5,999.5} } } },
			{"hNcontributionsToVertex", ";Number of contributors to vertex (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedTracksWithGeneratedParticle", ";Number of tracks in a collision with asociated generated particle (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedTracksWithoutGeneratedParticle", ";Number of tracks in a collision without asociated generated particle (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedPrimaryTracks", ";Number of primary tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedPrimaryParticlesOfReconstructedTracks", ";Number of primary generated particles in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedPrimaryElectronsOfReconstructedTracks", ";Number of primary generated electron in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedElectrons", ";Number of electrons in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNrecoVSgene", ";Number of reconstructed tracks in a collision (-);Number of generated particles in a collisions", { HistType::kTH2D, { {30,-0.5,29.5}, {30,-0.5,29.5} } } },
			{"hNrecoVSgeneMany", ";Number of reconstructed tracks in a collision (-);Number of generated particles in a collisions", { HistType::kTH2D, { {3000,-0.5,2999.5}, {3000,-0.5,2999.5} } } },
			{"hNgeneratedPrimaryElectrons", ";Number of primary electrons from mother decay in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedCollisions", ";Number of asociated reconstructed collisions in a generated collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedParticlesWithRecoColl", ";Number of generated particles in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedPrimaryParticlesWithRecoColl", ";Number of generated primary particles from mother decay in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedElectronsWithRecoColl", ";Number of generated electrons in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedPrimaryElectronsWithRecoColl", ";Number of generated primary electrons from mother decay in a reconstructed collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", { HistType::kTH1D, { {1200,2.8,3.4} } } },
			{"hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", { HistType::kTH1D, { {1200,1.0,5.} } } },
			{"hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,2.} } } },
			{"hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,10.} } } },
			{"hMotherPt", ";Mother #it{p_{#rm T}} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,2.} } } },
			{"hMotherPhi", ";Mother #phi (rad);Number of events (-)", { HistType::kTH1D, { {50,-TMath::Pi(),TMath::Pi()} } } },
			{"hMotherRapidity", ";Mother #y (-);Number of events (-)", { HistType::kTH1D, { {500,-2.,2.} } } },
			{"hVtxZ", ";Vertex z-position (cm);Number of events (-)", { HistType::kTH1D, { {1000,-30.,30.} } } },
			{"hVtxZrecToGen", ";Vertex z-position: (reconstructed collision) - (generated collisions) (cm);Number of events (-)", { HistType::kTH1D, { {1000,-.5,.5} } } },
			{"hVtxTransversal", ";Vertex x-position (cm);Vertex y-position (cm)", { HistType::kTH2D, { {1000,-0.1,0.1}, {1000,-0.1,0.1} } } },
			{"hTPCsignalVsMom", "All tracks;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCelectronIdentified", ";Track from generated electron #it{p} (GeV/#it{c});Track from generated electron TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 200,0.,2. }, { 950,10.,200. } } } },
			{"hTPCsignalVtxContributors", ";Vtx contributor 1 - TPC d#it{E}/d#it{x} (arb. units);Vtx contributor 2 - TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 950,10.,200. }, { 950,10.,200. } } } },
			{"hTPCsignalGeneratedElectrons", ";Track from generated electron 1 - TPC d#it{E}/d#it{x} (arb. units);Track from generated electron 2 - TPC d#it{E}/d#it{x} (arb. units)", { HistType::kTH2F, { { 950,10.,200. }, { 950,10.,200. } } } },
			{"hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,10.} } } },
			{"hTrackPt", ";Track #it{p_{#rm T}} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,5.} } } },
			{"hTrackPhi", ";Track #phi (rad);Number of events (-)", { HistType::kTH1D, { {50,-2*TMath::Pi(),2*TMath::Pi()} } } },
			{"hTrackEta", ";Track #eta (-);Number of events (-)", { HistType::kTH1D, { {500,-7.,7.} } } },
			{"hTrackEnergy", ";Track Energy (GeV);Number of events (-)", { HistType::kTH1D, { {500,0.,100.} } } },
			{"hTrackRapidity", ";Track rapidity (-);Number of events (-)", { HistType::kTH1D, { {500,-7.,7.} } } },
			{"hPDGcodes", ";PDG codes (-);Number of events (-)", { HistType::kTH1D, { {6001,-3000.,3000.} } } }
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
	using ReconstructedTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDMcTrackLabels, aod::UDTracksFlags>;
	using ReconstructedCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;

	// init
	void init(InitContext&){
		if (verboseInfo) printLargeMessage("INIT METHOD");
		countCollisions = 0;
		indexCollision = -1;
		isFirstReconstructedCollisions = true;
		isFirstGeneratedCollisions = true;

	} // end init

	// run (always called before process :( )
	void run(ProcessingContext& context){

		if (verboseInfo) printLargeMessage("RUN METHOD");
		if (verboseDebug) LOGF(info,"countCollisions = %d",countCollisions);

	} // end run

	// declare preslices = table cross-connections
	Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

	// process
	void processSimulatorLevel(ReconstructedCollision const& reconstructedCollision,
														 aod::UDMcCollisions const& generatedCollisions,
														 ReconstructedTracks const& reconstructedBarrelTracks,
	                           aod::UDMcParticles const& generatedParticles,
														 soa::Join<aod::Collisions, aod::EvSels> const& infoEventSelection){


		if(isFirstReconstructedCollisions){
			isFirstReconstructedCollisions = false;
			if (verboseInfo) printLargeMessage("START LOOPING OVER RECONSTRUCTED COLLISIONS");
		}
		if (verboseInfo) printLargeMessage("NEW COLLISION");
		countCollisions++;
		registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(0);

		auto thisCollisionEventSelectionInfo = infoEventSelection.iteratorAt(getEvSelsIndexOfThisCollisionBC(infoEventSelection,reconstructedCollision.globalBC()));
		fillEventSelectionHistogram(registry,thisCollisionEventSelectionInfo);
		fillTrackSelectionHistogram(registry,reconstructedBarrelTracks,reconstructedCollision);

		registry.get<TH1>(HIST("hNcontributionsToReconstructedCollision"))->Fill(reconstructedCollision.numContrib());

		// TODO fix when you can touch mcccollision // int nGeneratedParticlesInThisReconstructedCollision = ????????;
		int nReconstructedTracksInThisReconstructedCollision = reconstructedBarrelTracks.size();
		if (verboseDebug) LOGF(info,"collision number = %d has %d contributors",countCollisions,reconstructedCollision.numContrib());
		// TODO fix when you can touch mcccollision //if (verboseDebug) LOGF(info,"This slice: number of reconstructed tracks = %d, number of generated particles in this slice %d", reconstructedBarrelTracks.size(),nGeneratedParticlesInThisReconstructedCollision);

		// TODO fix when you can touch mcccollision //registry.get<TH1>(HIST("hNgeneratedParticles"))->Fill(nGeneratedParticlesInThisReconstructedCollision);
		registry.get<TH1>(HIST("hNreconstructedTracks"))->Fill(nReconstructedTracksInThisReconstructedCollision);
		registry.get<TH1>(HIST("hNreconstructedTracksWide"))->Fill(nReconstructedTracksInThisReconstructedCollision);
		// TODO fix when you can touch mcccollision //registry.get<TH2>(HIST("hNrecoVSgene"))->Fill(nReconstructedTracksInThisReconstructedCollision,nGeneratedParticlesInThisReconstructedCollision);
		// TODO fix when you can touch mcccollision //registry.get<TH2>(HIST("hNrecoVSgeneMany"))->Fill(nReconstructedTracksInThisReconstructedCollision,nGeneratedParticlesInThisReconstructedCollision);
		registry.get<TH1>(HIST("hVtxZ"))->Fill(reconstructedCollision.posZ());
		// TODO find a way how to fill in this histogram //registry.get<TH1>(HIST("hVtxZrecToGen"))->Fill(reconstructedCollision.posZ()-reconstructedCollision.mcCollision().posZ());
		registry.get<TH2>(HIST("hVtxTransversal"))->Fill(reconstructedCollision.posX(),reconstructedCollision.posY());

		int countReconstructedPrimaryTracks = 0;
		int countReconstructedTracksWithGeneratedParticle = 0;
		int countReconstructedTracksWithoutGeneratedParticle = 0;
		int countPrimaryGeneratedParticlesWithReconstructedTrack = 0;
		int countPrimaryGeneratedElectronsWithReconstructedTrack = 0;
		int countContributorsToVertex = 0;
		int countGlobalContributorsToVertex = 0;
		for (auto& reconstructedBarrelTrack : reconstructedBarrelTracks){
			float thisTrackMomentum = momentum(reconstructedBarrelTrack.px(),reconstructedBarrelTrack.py(),reconstructedBarrelTrack.pz());
			float thisTrackTPCsignal = reconstructedBarrelTrack.tpcSignal();
			if (!selectTrack(reconstructedBarrelTrack,applyTrackCuts)) continue;
			// TODO get tag for PV contributor //if (reconstructedBarrelTrack.isPVContributor()) countContributorsToVertex++;
			if (selectTrack(reconstructedBarrelTrack,2)) countGlobalContributorsToVertex++;
			if (reconstructedBarrelTrack.has_udMcParticle()) {
				if (reconstructedBarrelTrack.hasTPC()) registry.get<TH2>(HIST("hTPCsignalVsMom"))->Fill(thisTrackMomentum,thisTrackTPCsignal);
				countReconstructedTracksWithGeneratedParticle++;
				auto generatedParticle = reconstructedBarrelTrack.udMcParticle();
				if (generatedParticle.isPhysicalPrimary()) {
					countPrimaryGeneratedParticlesWithReconstructedTrack++;
					if (TMath::Abs(generatedParticle.pdgCode()) == 11) {
						countPrimaryGeneratedElectronsWithReconstructedTrack++;
						if (reconstructedBarrelTrack.hasTPC()) registry.get<TH2>(HIST("hTPCelectronIdentified"))->Fill(thisTrackMomentum,thisTrackTPCsignal);
					}
				}
			}
			else countReconstructedTracksWithoutGeneratedParticle++;
			// TODO get tag for PV contributor //if (reconstructedBarrelTrack.isPVContributor()) countReconstructedPrimaryTracks++;
		}

		if (printMetaInfo) printCollisionTracksData(reconstructedBarrelTracks,applyTrackCuts);

		registry.get<TH1>(HIST("hNcontributionsToVertex"))->Fill(countContributorsToVertex);
		registry.get<TH1>(HIST("hNreconstructedTracksWithGeneratedParticle"))->Fill(countReconstructedTracksWithGeneratedParticle);
		registry.get<TH1>(HIST("hNreconstructedTracksWithoutGeneratedParticle"))->Fill(countReconstructedTracksWithoutGeneratedParticle);
		registry.get<TH1>(HIST("hNreconstructedPrimaryTracks"))->Fill(countReconstructedPrimaryTracks);
		registry.get<TH1>(HIST("hNgeneratedPrimaryParticlesOfReconstructedTracks"))->Fill(countPrimaryGeneratedParticlesWithReconstructedTrack);
		registry.get<TH1>(HIST("hNgeneratedPrimaryElectronsOfReconstructedTracks"))->Fill(countPrimaryGeneratedElectronsWithReconstructedTrack);

		//
		// Playing around
		//
		printLargeMessage("Primary + Global + FIT");
		for (auto& reconstructedBarrelTrack : reconstructedBarrelTracks){
			if (!trackSelection(reconstructedBarrelTrack,2)) continue;
			if (!trackSelection(reconstructedBarrelTrack,5)) continue;
			if (isFITempty(reconstructedCollision)) continue;
			printTrackData(reconstructedBarrelTrack);
			float trkPx = reconstructedBarrelTrack.px();
			float trkPy = reconstructedBarrelTrack.py();
			float trkPz = reconstructedBarrelTrack.pz();
			registry.get<TH1>(HIST("hTrackP"))->Fill(momentum(trkPx,trkPy,trkPz));
			registry.get<TH1>(HIST("hTrackPt"))->Fill(reconstructedBarrelTrack.pt());
			registry.get<TH1>(HIST("hTrackPhi"))->Fill(phi(trkPx,trkPy));
			registry.get<TH1>(HIST("hTrackEta"))->Fill(eta(trkPx,trkPy,trkPz));
			if (reconstructedBarrelTrack.has_udMcParticle()) {
				auto mcparticle = reconstructedBarrelTrack.udMcParticle();
				registry.get<TH1>(HIST("hPDGcodes"))->Fill(mcparticle.pdgCode());
				auto pdgInfo = pdg->GetParticle(mcparticle.pdgCode());
				if (pdgInfo != nullptr) {
					float mass = pdgInfo->Mass();
					registry.get<TH1>(HIST("hTrackEnergy"))->Fill(energy(mass,trkPx,trkPy,trkPz));
					registry.get<TH1>(HIST("hTrackRapidity"))->Fill(rapidity(mass,trkPx,trkPy,trkPz));
				}
			}
		}

		//
		// fetching FT0, FDD, FV0 information
		//
		bool FITisEmpty = isFITempty(reconstructedCollision);
		if (FITisEmpty) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(1);
		FITisEmpty = isEvSelFITempty(thisCollisionEventSelectionInfo);
		if (FITisEmpty) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(2);
		if (countContributorsToVertex == 2) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(3);
		if (countGlobalContributorsToVertex == 2) registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(4);


		//
		// TPC signal info starts here
		//
		TLorentzVector sumOfVertexTracks, trackFromVertex;
		float dedxVC1 = -999.;
		float dedxVC2 = -999.;
		float dedxEl1 = -999.;
		float dedxEl2 = -999.;
		int countTrks = 0;
		if (countContributorsToVertex == 2) {
			for (auto& reconstructedBarrelTrack : reconstructedBarrelTracks){
				// TODO get tag for PV contributor //if (!reconstructedBarrelTrack.isPVContributor()) continue;
				countTrks++;
				trackFromVertex.SetPxPyPzE(reconstructedBarrelTrack.px(),reconstructedBarrelTrack.py(),reconstructedBarrelTrack.pz(),
																	 energy(constants::physics::MassElectron,reconstructedBarrelTrack.px(),reconstructedBarrelTrack.py(),reconstructedBarrelTrack.pz()));
				sumOfVertexTracks+=trackFromVertex;
				if (reconstructedBarrelTrack.hasTPC()){
					if (countTrks == 1) dedxVC1 = reconstructedBarrelTrack.tpcSignal();
					else if (countTrks == 2) dedxVC2 = reconstructedBarrelTrack.tpcSignal();
					else LOGP(warning,"!!!Strange behaviour in the loop, check!!!");
					if (reconstructedBarrelTrack.has_udMcParticle()) {
						auto generatedParticle = reconstructedBarrelTrack.udMcParticle();
						if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11) {
							if (countTrks == 1) dedxEl1 = reconstructedBarrelTrack.tpcSignal();
							else if (countTrks == 2) dedxEl2 = reconstructedBarrelTrack.tpcSignal();
							else LOGP(warning,"!!!Strange behaviour in the loop, check!!!");
						}
					}
				}
			}
			registry.get<TH2>(HIST("hTPCsignalVtxContributors"))->Fill(dedxVC1,dedxVC2);
			registry.get<TH2>(HIST("hTPCsignalGeneratedElectrons"))->Fill(dedxEl1,dedxEl2);
			registry.get<TH1>(HIST("hInvariantMass"))->Fill(sumOfVertexTracks.M());
			registry.get<TH1>(HIST("hInvariantMassWide"))->Fill(sumOfVertexTracks.M());
			registry.get<TH1>(HIST("hMotherP"))->Fill(sumOfVertexTracks.P());
			registry.get<TH1>(HIST("hMotherPwide"))->Fill(sumOfVertexTracks.P());
			registry.get<TH1>(HIST("hMotherPt"))->Fill(sumOfVertexTracks.Pt());
			registry.get<TH1>(HIST("hMotherPhi"))->Fill(sumOfVertexTracks.Phi());
			registry.get<TH1>(HIST("hMotherRapidity"))->Fill(sumOfVertexTracks.Rapidity());
		}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//   MC TRUTH STARTS HERE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO fix when you can touch mcccollision
//  if (printMetaInfo) printCollisionGeneratedParticles(generatedParticles,perMcCollision,reconstructedCollision);
//		TLorentzVector candidateJpsi, electron;
//		int countGeneratedElectrons = 0;
//		int countGeneratedPrimaryElectrons = 0;
//		int countGeneratedPrimaryParticles = 0;
//		for (auto & generatedParticle: ???){
//			if (verboseDebug) LOGF(info,"Particle type = %d",generatedParticle.pdgCode());
//			if (verboseDebug) LOGF(info,"Particle isPhysicalPrimary = %d",generatedParticle.isPhysicalPrimary());
//			if (generatedParticle.isPhysicalPrimary()) countGeneratedPrimaryParticles++;
//			if (TMath::Abs(generatedParticle.pdgCode()) == 11) {
//				countGeneratedElectrons++;
//			}
//			if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11){
//				countGeneratedPrimaryElectrons++;
//				electron.SetPxPyPzE(generatedParticle.px(),generatedParticle.py(),generatedParticle.pz(),generatedParticle.e());
//			}
//		}
//		if(countGeneratedPrimaryElectrons==2){
//			int electronCharge = 0;
//			for (auto & generatedParticle: ???){
//				if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11){
//					electron.SetPxPyPzE(generatedParticle.px(),generatedParticle.py(),generatedParticle.pz(),generatedParticle.e());
//					candidateJpsi+=electron;
//					if (electronCharge == 0) electronCharge = getElectronCharge(generatedParticle);
//					else if (electronCharge == getElectronCharge(generatedParticle))continue;
//				}
//			}
//		}
//		registry.get<TH1>(HIST("hNgeneratedElectronsWithRecoColl"))->Fill(countGeneratedElectrons);
//		registry.get<TH1>(HIST("hNgeneratedPrimaryElectronsWithRecoColl"))->Fill(countGeneratedPrimaryElectrons);
//		registry.get<TH1>(HIST("hNgeneratedPrimaryParticlesWithRecoColl"))->Fill(countGeneratedPrimaryParticles);


	} // end processSimulatorLevel

	void processGeneratorLevel(aod::McCollision const& generatedCollision,
														 soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& reconstructedCollisions,
														 aod::McParticles const& generatedParticles){

		if(isFirstGeneratedCollisions){
			isFirstGeneratedCollisions = false;
			if (verboseInfo) LOGF(info,"####################################### START LOOPING OVER GENERATED COLLISIONS #######################################");
		}

		if (verboseInfo) LOGF(info,"####################################### NEW COLLISION #######################################");
		countCollisions++;
		int countElectrons = 0;
		int countElectronsWithRecoColl = 0;
		int countPrimaryElectrons = 0;
		int countPrimaryElectronsWithRecoColl = 0;
		if (verboseDebug) LOGF(info,"countCollisions = %d, number of particles in this collision %d",countCollisions,generatedParticles.size());
		if (verboseDebug) LOGF(info, "(truth) Generators id = %d, MC weight = %f", generatedCollision.generatorsID(), generatedCollision.weight());
		if (verboseDebug) LOGF(info, "(reco) collisions.size = %d", reconstructedCollisions.size());
		if (verboseDebug) LOGF(info, "(truth) BC id = %d, collision time = %f,  impact parameter = %f", generatedCollision.bcId(), generatedCollision.t(), generatedCollision.impactParameter());
		if (verboseDebug) LOGF(info, "(truth) vtx-x = %f,  vtx-y = %f,  vtx-z = %f", generatedCollision.posX(), generatedCollision.posY(), generatedCollision.posZ());

		registry.get<TH1>(HIST("hNgeneratedParticles"))->Fill(generatedParticles.size());
		registry.get<TH1>(HIST("hVtxZ"))->Fill(generatedCollision.posZ());
		registry.get<TH2>(HIST("hVtxTransversal"))->Fill(generatedCollision.posX(),generatedCollision.posY());

		registry.get<TH1>(HIST("hNreconstructedCollisions"))->Fill(reconstructedCollisions.size());
		if (reconstructedCollisions.size()>0) {
			for (auto& reconstructedCollision : reconstructedCollisions) {
				registry.get<TH1>(HIST("hNgeneratedParticlesWithRecoColl"))->Fill(generatedParticles.size());
				registry.get<TH1>(HIST("hVtxZrecToGen"))->Fill(reconstructedCollision.posZ()-generatedCollision.posZ());
			}
		}

		TLorentzVector candidateJpsi, electron;
		for (auto & generatedParticle: generatedParticles){
			if (verboseDebug) LOGF(info,"Particle type = %d",generatedParticle.pdgCode());
			if (verboseDebug) LOGF(info,"Particle isPhysicalPrimary = %d",generatedParticle.isPhysicalPrimary());
			if (TMath::Abs(generatedParticle.pdgCode()) == 11) {
				countElectrons++;
				if (reconstructedCollisions.size()>0) countElectronsWithRecoColl++;
			}
			if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11){
				countPrimaryElectrons++;
				if (reconstructedCollisions.size()>0) countPrimaryElectronsWithRecoColl++;
				electron.SetPxPyPzE(generatedParticle.px(),generatedParticle.py(),generatedParticle.pz(),generatedParticle.e());
				candidateJpsi+=electron;
			}
		}
		registry.get<TH1>(HIST("hNgeneratedElectrons"))->Fill(countElectrons);
		registry.get<TH1>(HIST("hNgeneratedElectronsWithRecoColl"))->Fill(countElectronsWithRecoColl);
		registry.get<TH1>(HIST("hNgeneratedPrimaryElectrons"))->Fill(countPrimaryElectrons);
		registry.get<TH1>(HIST("hNgeneratedPrimaryElectronsWithRecoColl"))->Fill(countPrimaryElectronsWithRecoColl);
		registry.get<TH1>(HIST("hInvariantMass"))->Fill(candidateJpsi.M());
		registry.get<TH1>(HIST("hMotherP"))->Fill(candidateJpsi.P());
		registry.get<TH1>(HIST("hMotherPwide"))->Fill(candidateJpsi.P());
		registry.get<TH1>(HIST("hMotherPt"))->Fill(candidateJpsi.Pt());
		registry.get<TH1>(HIST("hMotherPhi"))->Fill(candidateJpsi.Phi());
		registry.get<TH1>(HIST("hMotherRapidity"))->Fill(candidateJpsi.Rapidity());


	} // end processGeneratorLevel

	void processTest1(aod::McCollision const& generatedCollision, aod::McParticles const& generatedParticles){
		if (verboseDebug) LOGF(info,"n particles %d",generatedParticles.size());
	}// end processTest1

	void processTest2(aod::McCollision const& generatedCollision){

		LOGF(info,"MC generators ID %d",generatedCollision.generatorsID());

	}// end processTest2

	void processTest3(aod::McParticles const& generatedParticles){
		if (verboseDebug) LOGF(info,"n particles %d",generatedParticles.size());
	}// end processTest3

	void processTest4(aod::McParticle const& generatedParticle){

		if (generatedParticle.mcCollisionId() < indexCollision) LOGF(warning,"WARNING: previous index %d, this index %d",indexCollision,generatedParticle.mcCollisionId());
		indexCollision = generatedParticle.mcCollisionId();
	}// end processTest4

	int largestIndex = 0;

	void processTest5(soa::Join<aod::Collisions,aod::McCollisionLabels>::iterator const& collision,
	                                aod::McCollisions const& generatedCollisions){

		if (collision.has_mcCollision()) {
			auto generatedCollision = collision.mcCollision();
			if (largestIndex <= generatedCollision.globalIndex()) largestIndex = generatedCollision.globalIndex();
		}

	}// end processTest5

	void processTest6(aod::Collision const& collision){

		LOGF(info,"Collision BC id %d",collision.bcId());

	}// end processTest6


	void processTestUDtable1(soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDMcTrackLabels, aod::UDTracksFlags>::iterator const& track,
													 aod::UDMcParticles const& particles){

		LOGF(info,"IDX %d, mask %d",track.udMcParticleId(),track.mcMask());
		printTrackParticleData(track);

	}// end processTestUDtable1

	void processTestUDtable2(soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA>::iterator const& track){

//		printTrackData(track);

		if (isUDglobalTrack(track)) {
			LOGF(info,"++++++");
			LOGF(info,"isUDprimaryTrack(track) %b",isUDprimaryTrack(track));
			LOGF(info,"isUDinAcceptanceTrack(track) %b",isUDinAcceptanceTrack(track));
			LOGF(info,"isUDqualityTrack(track) %b",isUDqualityTrack(track));
		}


	}// end processTestUDtable2

	void processTestUDtable3(aod::UDCollision const& collision){

		LOGF(info,"IDX %d, globalBC %d",collision.globalIndex(),collision.globalBC());

	}// end processTestUDtable3


	PROCESS_SWITCH(UpcMCCentralBarrel, processSimulatorLevel, "Iterate MC tables with reconstructed data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processGeneratorLevel, "Iterate MC tables with generated data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest1, "Test MC collision + particles", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest2, "Test MC collisions", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest3, "Test MC particles", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest4, "Test single MC particle", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest5, "Test index sorting of MC collisions in reco collision", false)
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest6, "Test single collision", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTestUDtable1, "Test tracks and its particles from UD table", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTestUDtable2, "Test tracks from UD table", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTestUDtable3, "Test collisions from UD table", false);

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<UpcMCCentralBarrel>(cfgc, TaskName{"upc-mc-central-barrel"})
	};
}