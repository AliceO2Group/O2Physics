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
//#include "CommonConstants/PhysicsConstants.h"

// O2Physics headers
//#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

// ROOT headers
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename T>
int getElectronCharge(T generatedElectron)
// Check if particle is electron or positron and return charge accordingly. Return zero if particle is not electron/positron
{
	if (generatedElectron.pdgCode() == 11) return -1;
	else if (generatedElectron.pdgCode() == -11) return 1;
	else return 0;
}

struct UpcMCCentralBarrel {

	// helper function
	float invariantMass(float E, float px, float py, float pz){
		float im = E*E - px*px - py*py - pz*py;
		return TMath::Sqrt(im);
	}
	// Global varialbes
	bool isFirstReconstructedCollisions, isFirstGeneratedCollisions;
	int countCollisions;
	int32_t indexCollision;

		HistogramRegistry registry{
		"registry",
		{
			{"hEffectOfSelections", "Effect of cuts;Selection (-);Number of events (-)", { HistType::kTH1D, { {10,-0.5,9.5} } } },
			{"hNcontributionsToReconstructedCollision", ";Number of contributions in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNgeneratedParticles", ";Number of particles in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
			{"hNreconstructedTracks", ";Number of tracks in a collision (-);Number of events (-)", { HistType::kTH1D, { {30,-0.5,29.5} } } },
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
			{"hMotherPt", ";Mother #it{p_{#rm T}} (GeV/c);Number of events (-)", { HistType::kTH1D, { {500,0.,2.} } } },
			{"hMotherPhi", ";Mother #phi (rad);Number of events (-)", { HistType::kTH1D, { {50,-TMath::Pi(),TMath::Pi()} } } },
			{"hMotherRapidity", ";Mother #y (-);Number of events (-)", { HistType::kTH1D, { {500,-2.,2.} } } },
			{"hVtxZ", ";Vertex z-position (cm);Number of events (-)", { HistType::kTH1D, { {1000,-30.,30.} } } },
			{"hVtxZrecToGen", ";Vertex z-position: (reconstructed collision) - (generated collisions) (cm);Number of events (-)", { HistType::kTH1D, { {1000,-.5,.5} } } },
			{"hVtxTransversal", ";Vertex x-position (cm);Vertex y-position (cm)", { HistType::kTH2D, { {1000,-0.1,0.1}, {1000,-0.1,0.1} } } },
		}
	};

	// declare configurables
	Configurable<bool> verboseInfo{"verboseInfo", true, {"Print general info to terminal; default it true."}};
	Configurable<bool> verboseDebug{"verboseDebug", false, {"Print debug info to terminal; default it false."}};
	Configurable<bool> verboseCurrent{"verboseCurrent", false, {"Print current info to terminal; default it false."}};

	// declare filters
	Filter nCollisionContributorsFilter = aod::collision::numContrib > 2;

	// declare shortcuts
	using ReconstructedTCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels>;
	using ReconstructedCollision = soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator;
	//	using ReconstructedCollision = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>>::iterator;


	// init
	void init(InitContext&){
		if (verboseInfo) LOGF(info,"####################################### INIT METHOD #######################################");
		countCollisions = 0;
		indexCollision = -1;
		isFirstReconstructedCollisions = true;
		isFirstGeneratedCollisions = true;
	} // end init

	// run (always called before process :( )
	void run(ProcessingContext& context){

		if (verboseInfo) LOGF(info,"####################################### RUN METHOD #######################################");
		if (verboseDebug) LOGF(info,"countCollisions = %d",countCollisions);

	} // end run

	// declare preslices = table cross-connections
	Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
	Preslice<aod::McCollisions> perBC = aod::mccollision::bcId;
	Preslice<aod::Tracks> perCollision = aod::track::collisionId;

	// process
	void processSimulatorLevel(ReconstructedCollision const& reconstructedCollision,
														 aod::McCollisions const& generatedCollisions,
														 ReconstructedTCs const& reconstructedBarrelTracks,
	                           aod::McParticles const& generatedParticles){

		if(isFirstReconstructedCollisions){
			isFirstReconstructedCollisions = false;
			if (verboseInfo) LOGF(info,"####################################### START LOOPING OVER RECONSTRUCTED COLLISIONS #######################################");
		}

		registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(0);
		registry.get<TH1>(HIST("hNcontributionsToReconstructedCollision"))->Fill(reconstructedCollision.numContrib());

		auto thisReconstructedCollisionGeneratedParticles = generatedParticles.sliceBy(perMcCollision, reconstructedCollision.mcCollision().globalIndex());
		auto thisBunchCrossingGeneratedCollisions = generatedCollisions.sliceBy(perBC,reconstructedCollision.bcId());
//		auto thisReconstructedCollisionReconstructedTracks = reconstructedBarrelTracks.sliceBy(perCollision, reconstructedCollision.globalIndex());

		if (verboseInfo) LOGF(info,"####################################### NEW COLLISION #######################################");
		countCollisions++;
		if (verboseCurrent) LOGF(info,"Global idx of reconstructedCollision %d; global idx of associated reconstructedCollision.mcCollision() %d",
														 reconstructedCollision.globalIndex(),reconstructedCollision.mcCollision().globalIndex());
		if (verboseCurrent) LOGF(info,"%d generated collisions in bunch crossing %d related to this reconstructed collisions",
														 thisBunchCrossingGeneratedCollisions.size(),reconstructedCollision.bcId());

		int nGeneratedParticlesInThisReconstructedCollision = thisReconstructedCollisionGeneratedParticles.size();
		int nReconstructedTracksInThisReconstructedCollision = reconstructedBarrelTracks.size();
//		int nReconstructedTracksInThisReconstructedCollisionAlternative = thisReconstructedCollisionReconstructedTracks.size();
		if (verboseDebug) LOGF(info,"collision number = %d has %d contributors",countCollisions,reconstructedCollision.numContrib());
		if (verboseDebug) LOGF(info,"This slice: number of reconstructed tracks = %d, number of generated particles in this slice %d",
													reconstructedBarrelTracks.size(),nGeneratedParticlesInThisReconstructedCollision);
//		if (verboseCurrent) LOGF(info,"%d reconstructed tracks according to preslice",nReconstructedTracksInThisReconstructedCollisionAlternative);
//
		registry.get<TH1>(HIST("hNgeneratedParticles"))->Fill(nGeneratedParticlesInThisReconstructedCollision);
		registry.get<TH1>(HIST("hNreconstructedTracks"))->Fill(nReconstructedTracksInThisReconstructedCollision);
		registry.get<TH2>(HIST("hNrecoVSgene"))->Fill(nReconstructedTracksInThisReconstructedCollision,nGeneratedParticlesInThisReconstructedCollision);
		registry.get<TH2>(HIST("hNrecoVSgeneMany"))->Fill(nReconstructedTracksInThisReconstructedCollision,nGeneratedParticlesInThisReconstructedCollision);
		registry.get<TH1>(HIST("hVtxZ"))->Fill(reconstructedCollision.posZ());
		registry.get<TH1>(HIST("hVtxZrecToGen"))->Fill(reconstructedCollision.posZ()-reconstructedCollision.mcCollision().posZ());
		registry.get<TH2>(HIST("hVtxTransversal"))->Fill(reconstructedCollision.posX(),reconstructedCollision.posY());

		int countReconstructedPrimaryTracks = 0;
		int countReconstructedTracksWithGeneratedParticle = 0;
		int countReconstructedTracksWithoutGeneratedParticle = 0;
		int countPrimaryGeneratedParticlesWithReconstructedTrack = 0;
		int countPrimaryGeneratedElectronsWithReconstructedTrack = 0;
		int countContributorsToVertex = 0;
		for (auto& reconstructedBarrelTrack : reconstructedBarrelTracks){
//			if (!reconstructedBarrelTrack.has_mcParticle())  {
//				LOGF(warning, "No MC particle for track, skip...");
//				continue;
//			}
			if (reconstructedBarrelTrack.isPVContributor()) countContributorsToVertex++;
			if (reconstructedBarrelTrack.has_mcParticle()) {
				countReconstructedTracksWithGeneratedParticle++;
				auto generatedParticle = reconstructedBarrelTrack.mcParticle();
				if (verboseCurrent) LOGF(info,"Global idx of track in reconstructedBarrelTrack %d, is Vtx contributor %d, global idx of particle in reconstructedBarrelTrack.mcParticle() %d, is Primary %d",
																 reconstructedBarrelTrack.globalIndex(),reconstructedBarrelTrack.isPVContributor(),generatedParticle.globalIndex(),generatedParticle.isPhysicalPrimary());
				if (generatedParticle.isPhysicalPrimary()) {
					countPrimaryGeneratedParticlesWithReconstructedTrack++;
					if (TMath::Abs(generatedParticle.pdgCode()) == 11) countPrimaryGeneratedElectronsWithReconstructedTrack++;
				}
			}
			else countReconstructedTracksWithoutGeneratedParticle++;
			if (reconstructedBarrelTrack.isPVContributor()) countReconstructedPrimaryTracks++;
		}

//		for (auto& reconstructedBarrelTrack : thisReconstructedCollisionReconstructedTracks){
//			if (reconstructedBarrelTrack.has_mcParticle()) {
//				auto generatedParticle = reconstructedBarrelTrack.mcParticle();
//				if (verboseCurrent) LOGF(info,"Global idx of track in thisReconstructedCollisionReconstructedTracks %d",reconstructedBarrelTrack.globalIndex());
//				if (verboseCurrent) LOGF(info,"Global idx of particle in thisReconstructedCollisionReconstructedTracks.mcParticle() %d, is Primary %d",generatedParticle.globalIndex(),generatedParticle.isPhysicalPrimary());
//			}
//		}

		registry.get<TH1>(HIST("hNcontributionsToVertex"))->Fill(countContributorsToVertex);
		registry.get<TH1>(HIST("hNreconstructedTracksWithGeneratedParticle"))->Fill(countReconstructedTracksWithGeneratedParticle);
		registry.get<TH1>(HIST("hNreconstructedTracksWithoutGeneratedParticle"))->Fill(countReconstructedTracksWithoutGeneratedParticle);
		registry.get<TH1>(HIST("hNreconstructedPrimaryTracks"))->Fill(countReconstructedPrimaryTracks);
		registry.get<TH1>(HIST("hNgeneratedPrimaryParticlesOfReconstructedTracks"))->Fill(countPrimaryGeneratedParticlesWithReconstructedTrack);
		registry.get<TH1>(HIST("hNgeneratedPrimaryElectronsOfReconstructedTracks"))->Fill(countPrimaryGeneratedElectronsWithReconstructedTrack);

		TLorentzVector sumOfVertexTracks, trackFromVertex;
		if (countContributorsToVertex == 2) {
			for (auto& reconstructedBarrelTrack : reconstructedBarrelTracks){
				if (!reconstructedBarrelTrack.isPVContributor()) continue;
				trackFromVertex.SetPxPyPzE(reconstructedBarrelTrack.px(),reconstructedBarrelTrack.py(),reconstructedBarrelTrack.pz(),
																	 reconstructedBarrelTrack.energy(constants::physics::MassElectron));
				sumOfVertexTracks+=trackFromVertex;
			}
		}

		registry.get<TH1>(HIST("hInvariantMass"))->Fill(sumOfVertexTracks.M());
		registry.get<TH1>(HIST("hInvariantMassWide"))->Fill(sumOfVertexTracks.M());
		registry.get<TH1>(HIST("hMotherP"))->Fill(sumOfVertexTracks.P());
		registry.get<TH1>(HIST("hMotherPt"))->Fill(sumOfVertexTracks.Pt());
		registry.get<TH1>(HIST("hMotherPhi"))->Fill(sumOfVertexTracks.Phi());
		registry.get<TH1>(HIST("hMotherRapidity"))->Fill(sumOfVertexTracks.Rapidity());


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//   MC TRUTH STARTS HERE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (auto& generatedParticle : thisReconstructedCollisionGeneratedParticles){
			if (verboseCurrent) LOGF(info,"Global idx of particle in thisReconstructedCollisionGeneratedParticles %d, is Primary %d",generatedParticle.globalIndex(),generatedParticle.isPhysicalPrimary());
		}

		TLorentzVector candidateJpsi, electron;
		int countGeneratedElectrons = 0;
		int countGeneratedPrimaryElectrons = 0;
		int countGeneratedPrimaryParticles = 0;
		for (auto & generatedParticle: thisReconstructedCollisionGeneratedParticles){
			if (verboseDebug) LOGF(info,"Particle type = %d",generatedParticle.pdgCode());
			if (verboseDebug) LOGF(info,"Particle isPhysicalPrimary = %d",generatedParticle.isPhysicalPrimary());
			if (generatedParticle.isPhysicalPrimary()) countGeneratedPrimaryParticles++;
			if (TMath::Abs(generatedParticle.pdgCode()) == 11) {
				registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(1);
				countGeneratedElectrons++;
			}
			if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11){
				registry.get<TH1>(HIST("hEffectOfSelections"))->Fill(2);
				countGeneratedPrimaryElectrons++;
				electron.SetPxPyPzE(generatedParticle.px(),generatedParticle.py(),generatedParticle.pz(),generatedParticle.e());
			}
		}
		if(countGeneratedPrimaryElectrons==2){
			int electronCharge = 0;
			for (auto & generatedParticle: thisReconstructedCollisionGeneratedParticles){
				if (generatedParticle.isPhysicalPrimary() && TMath::Abs(generatedParticle.pdgCode()) == 11){
					electron.SetPxPyPzE(generatedParticle.px(),generatedParticle.py(),generatedParticle.pz(),generatedParticle.e());
					candidateJpsi+=electron;
					if (electronCharge == 0) electronCharge = getElectronCharge(generatedParticle);
					else if (electronCharge == getElectronCharge(generatedParticle))continue;
				}
			}
		}
		registry.get<TH1>(HIST("hNgeneratedElectronsWithRecoColl"))->Fill(countGeneratedElectrons);
		registry.get<TH1>(HIST("hNgeneratedPrimaryElectronsWithRecoColl"))->Fill(countGeneratedPrimaryElectrons);
		registry.get<TH1>(HIST("hNgeneratedPrimaryParticlesWithRecoColl"))->Fill(countGeneratedPrimaryParticles);

//		registry.get<TH1>(HIST("hInvariantMass"))->Fill(candidateJpsi.M());
//		registry.get<TH1>(HIST("hMotherP"))->Fill(candidateJpsi.P());
//		registry.get<TH1>(HIST("hMotherPt"))->Fill(candidateJpsi.Pt());
//		registry.get<TH1>(HIST("hMotherPhi"))->Fill(candidateJpsi.Phi());
//		registry.get<TH1>(HIST("hMotherRapidity"))->Fill(candidateJpsi.Rapidity());

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
		registry.get<TH1>(HIST("hMotherPt"))->Fill(candidateJpsi.Pt());
		registry.get<TH1>(HIST("hMotherPhi"))->Fill(candidateJpsi.Phi());
		registry.get<TH1>(HIST("hMotherRapidity"))->Fill(candidateJpsi.Rapidity());



	} // end processGeneratorLevel

	void processTest1(aod::McCollision const& generatedCollision,
	                           aod::McParticles const& generatedParticles){
		if (verboseDebug) LOGF(info,"n particles %d",generatedParticles.size());

	}// end processGeneratorLevelTest

	void processTest2(aod::McCollision const& generatedCollision){

	}// end processGeneratorLevelTest

	void processTest3(aod::McParticles const& generatedParticles){
		if (verboseDebug) LOGF(info,"n particles %d",generatedParticles.size());
	}// end processGeneratorLevelTest

	void processTest4(aod::McParticle const& generatedParticle){

		if (generatedParticle.mcCollisionId() < indexCollision) LOGF(warning,"WARNING: previous index %d, this index %d",indexCollision,generatedParticle.mcCollisionId());
		indexCollision = generatedParticle.mcCollisionId();
	}// end processGeneratorLevelTest

	int largestIndex = 0;

	void processTest5(soa::Join<aod::Collisions,aod::McCollisionLabels>::iterator const& collision,
	                                aod::McCollisions const& generatedCollisions){

		if (collision.has_mcCollision()) {
			auto generatedCollision = collision.mcCollision();
			if (largestIndex <= generatedCollision.globalIndex()) largestIndex = generatedCollision.globalIndex();
		}

	}// end processGeneratorLevelTest5

	void processAnalysisFinished(aod::Collisions const& collisions){

		if (verboseInfo) LOGF(info,"####################################### END #######################################");
		if (verboseDebug) LOGF(info,"countCollisions = %d, largestIndex = %d",countCollisions,largestIndex);
		isFirstReconstructedCollisions = true;
		isFirstGeneratedCollisions = true;

	} // end processAnalysisFinished

	PROCESS_SWITCH(UpcMCCentralBarrel, processSimulatorLevel, "Iterate MC tables with reconstructed data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processGeneratorLevel, "Iterate MC tables with generated data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest1, "Test with data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest2, "Test with data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest3, "Test with data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest4, "Test with data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processTest5, "Test with data", false);
	PROCESS_SWITCH(UpcMCCentralBarrel, processAnalysisFinished, "Simply runs in the end", true);

};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
	return WorkflowSpec{
		adaptAnalysisTask<UpcMCCentralBarrel>(cfgc, TaskName{"upc-mc-central-barrel"})
	};
}