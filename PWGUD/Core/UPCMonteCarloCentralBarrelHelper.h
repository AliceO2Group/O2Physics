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
/// \since  27.10.2022

#ifndef ALISW_UPCMONTECARLOCENTRALBARRELHELPER_H
#define ALISW_UPCMONTECARLOCENTRALBARRELHELPER_H


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename C>
bool isEvSelFITempty(C collision)
// Get FIT information from EventSelection task for each collision
{
	if (collision.has_foundFT0() || collision.has_foundFV0() || collision.has_foundFDD()) return false;
	else return true;
}

template <typename T>
bool trackSelection(T track, int selection)
// Do selection of reconstructed track
{
	if (selection==0) return true;
	// Is central barrel propagated track
	if (selection==1 && track.trackType()!=1) return false;
	// Is central barrel vertex contributor
	if (selection==2 && track.isPVContributor()!=1) return false;
	// Is central barrel track selection global track
	if (selection==3 && track.isQualityTrack()!=1) return false;
	// Is central barrel track selection global track
	if (selection==4 && track.isPrimaryTrack()!=1) return false;
	// Is central barrel track selection global track
	if (selection==5 && track.isGlobalTrack()!=1) return false;

	return true;
}

template <typename T>
bool selectTrack(T track, int setOfCuts)
// Do selection of reconstructed track
{
	if (setOfCuts<1) return true;
	// Is central barrel propagated track)
	if (setOfCuts<2 && trackSelection(track,1)!=1) return false;
	// Is central barrel vertex contributor
	if (setOfCuts<3 && trackSelection(track,2)!=1) return false;

	return true;
}

template <typename E>
int getElectronCharge(E generatedElectron)
// Check if particle is electron or positron and return charge accordingly. Return zero if particle is not electron/positron
{
	if (generatedElectron.pdgCode() == 11) return -1;
	else if (generatedElectron.pdgCode() == -11) return 1;
	else return 0;
}

template < typename C>
void fillEventSelectionHistogram(HistogramRegistry &registry, C collision)
// Fill into histogram information from EventSelection task for each collision
{

	for (int i = 0; i < kNsel; i++) {
		if(collision.selection()[i]) registry.get<TH1>(HIST("hEventSelectionTaskParameters"))->Fill(i);
	}

	const int nXbins = registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->GetNbins();
	registry.get<TH1>(HIST("hEventSelectionTask"))->SetNdivisions(nXbins, "X");
	if (collision.sel7()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(0);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(1,70,0.02,30,-1,-1,"sel7");
	if (collision.sel8()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(1);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(2,70,0.02,30,-1,-1,"sel8");
	if (collision.bbV0A()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(2);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(3,70,0.02,30,-1,-1,"bbV0A");
	if (collision.bbV0C()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(3);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(4,70,0.02,30,-1,-1,"bbV0C");
	if (collision.bgV0A()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(4);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(5,70,0.02,30,-1,-1,"bgV0A");
	if (collision.bgV0C()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(5);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(6,70,0.02,30,-1,-1,"bgV0C");
	if (collision.bbFDA()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(6);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(7,70,0.02,30,-1,-1,"bbFDA");
	if (collision.bbFDC()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(7);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(8,70,0.02,30,-1,-1,"bbFDC");
	if (collision.bgFDA()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(8);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(9,70,0.02,30,-1,-1,"bgFDA");
	if (collision.bgFDC()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(9);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(10,70,0.02,30,-1,-1,"bgFDC");
	if (collision.has_foundBC()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(10);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(11,70,0.02,30,-1,-1,"has_foundBC");
	if (collision.has_foundFT0()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(11);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(12,70,0.02,30,-1,-1,"has_foundFT0");
	if (collision.has_foundFV0()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(12);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(13,70,0.02,30,-1,-1,"has_foundFV0");
	if (collision.has_foundFDD()) registry.get<TH1>(HIST("hEventSelectionTask"))->Fill(13);
	registry.get<TH1>(HIST("hEventSelectionTask"))->GetXaxis()->ChangeLabel(14,70,0.02,30,-1,-1,"has_foundFDD");

}

template <typename Ts, typename C>
void fillTrackSelectionHistogram(HistogramRegistry &registry, Ts tracks, C collision)
// Fill into histogram effect of track selection for all tracks
{

	const int nXbins = registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->GetNbins();
	registry.get<TH1>(HIST("hEffectOfTrackSelections"))->SetNdivisions(nXbins, "X");

	for (auto& track : tracks){
		if (trackSelection(track,0)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(0);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(1,70,0.02,30,-1,-1,"no cut");
		if (trackSelection(track,1)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(1);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(2,70,0.02,30,-1,-1,"propagated");
		if (trackSelection(track,3)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(2);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(3,70,0.02,30,-1,-1,"quality");
		if (trackSelection(track,4)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(3);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(4,70,0.02,30,-1,-1,"primary");
		if (trackSelection(track,5)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(4);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(5,70,0.02,30,-1,-1,"global");
		if (trackSelection(track,2)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(5);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(6,70,0.02,30,-1,-1,"PV contributor");


		if (!trackSelection(track,2)) continue;
		if (!trackSelection(track,5)) continue;
		bool hasITS = track.hasITS();
		bool hasTPC = track.hasTPC();
		bool hasTOF = track.hasTOF();
		if (hasITS) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(10);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(11,70,0.02,30,-1,-1,"ITS");
		if (hasTPC) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(11);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(12,70,0.02,30,-1,-1,"TPC");
		if (hasTOF) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(12);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(13,70,0.02,30,-1,-1,"TOF");
		if (hasITS && hasTPC) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(13);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(14,70,0.02,30,-1,-1,"ITS+TPC");
		if (hasITS && hasTOF) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(14);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(15,70,0.02,30,-1,-1,"ITS+TOF");
		if (hasTPC && hasTOF) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(15);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(16,70,0.02,30,-1,-1,"TPC+TOF");
		if (hasITS && hasTPC && hasTOF) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(16);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(17,70,0.02,30,-1,-1,"ITS+TPC+TOF");
		if (isEvSelFITempty(collision)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(17);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(18,70,0.02,30,-1,-1,"FIT empty");

	}
}

template <typename FTs, typename FDs, typename FVs>
bool isFITempty(uint64_t bc,
                std::map<uint64_t, int32_t>& bcsWithFT0,
                std::map<uint64_t, int32_t>& bcsWithFDD,
                std::map<uint64_t, int32_t>& bcsWithFV0A,
                FTs ft0s, FDs fdds, FVs fv0as)
// Return true if FIT had no signal
{
	// FIT beam-gas flags
	bool isBGFT0A, isBGFT0C, isBGFV0A, isBGFDDA, isBGFDDC;
	// FIT beam-beam flags
	bool isBBFT0A, isBBFT0C, isBBFV0A, isBBFDDA, isBBFDDC;

	// use "default" parameters
	float fV0ABBlower = -3.0;  // ns
	float fV0ABBupper = +2.0;  // ns
	float fV0ABGlower = 2.0;   // ns
	float fV0ABGupper = 5.0;   // ns

	float fFDABBlower = -3.0;  // ns
	float fFDABBupper = +3.0;  // ns
	float fFDABGlower = 10.0;  // ns
	float fFDABGupper = 13.0;  // ns
	float fFDCBBlower = -3.0;  // ns
	float fFDCBBupper = +3.0;  // ns
	float fFDCBGlower = -10.0; // ns
	float fFDCBGupper = -3.0;  // ns

	float fT0ABBlower = -1.0;  // ns
	float fT0ABBupper = +1.0;  // ns
	float fT0CBBlower = -1.0;  // ns
	float fT0CBBupper = +1.0;  // ns
	float fT0ABGlower = +1.0;  // ns
	float fT0ABGupper = +4.0;  // ns
	float fT0CBGlower = -4.0;  // ns
	float fT0CBGupper = -1.0;  // ns

	float timeFV0A = -999.f;
	float timeFT0A = -999.f;
	float timeFT0C = -999.f;
	float timeFDDA = -999.f;
	float timeFDDC = -999.f;
	float timeV0ABG = -999.f;
	float timeT0ABG = -999.f;
	float timeT0CBG = -999.f;
	float timeFDABG = -999.f;
	float timeFDCBG = -999.f;

	// check FIT info in the same BC
	auto it = bcsWithFT0.find(bc);
	if (it != bcsWithFT0.end()) {
		const auto& ft0 = ft0s.iteratorAt(it->second);
		timeFT0A = ft0.timeA();
		timeFT0C = ft0.timeC();
	}

	it = bcsWithFDD.find(bc);
	if (it != bcsWithFDD.end()) {
		const auto& fdd = fdds.iteratorAt(it->second);
		timeFDDA = fdd.timeA();
		timeFDDC = fdd.timeC();
	}

	it = bcsWithFV0A.find(bc);
	if (it != bcsWithFV0A.end()) {
		const auto& fv0a = fv0as.iteratorAt(it->second);
		timeFV0A = fv0a.time();
	}

	// check beam-gas
	it = bcsWithFT0.find(bc - 1);
	if (it != bcsWithFT0.end()) {
		const auto& ft0 = ft0s.iteratorAt(it->second);
		timeT0ABG = ft0.timeA();
		timeT0CBG = ft0.timeC();
	}

	it = bcsWithFDD.find(bc - 5);
	if (it != bcsWithFDD.end()) {
		const auto& ft0 = fdds.iteratorAt(it->second);
		timeFDABG = ft0.timeA();
		timeFDCBG = ft0.timeC();
	}

	it = bcsWithFV0A.find(bc - 1);
	if (it != bcsWithFV0A.end()) {
		const auto& fv0a = fv0as.iteratorAt(it->second);
		timeV0ABG = fv0a.time();
	}

	// beam-gas flags
	isBGFV0A = timeV0ABG > fV0ABGlower && timeV0ABG < fV0ABGupper;
	isBGFDDA = timeFDABG > fFDABGlower && timeFDABG < fFDABGupper;
	isBGFDDC = timeFDCBG > fFDCBGlower && timeFDCBG < fFDCBGupper;
	isBGFT0A = timeT0ABG > fT0ABGlower && timeT0ABG < fT0ABGupper;
	isBGFT0C = timeT0CBG > fT0CBGlower && timeT0CBG < fT0CBGupper;

	// beam-beam flags
	isBBFT0A = timeFT0A > fT0ABBlower && timeFT0A < fT0ABBupper;
	isBBFT0C = timeFT0C > fT0CBBlower && timeFT0C < fT0CBBupper;
	isBBFV0A = timeFV0A > fV0ABBlower && timeFV0A < fV0ABBupper;
	isBBFDDA = timeFDDA > fFDABBlower && timeFDDA < fFDABBupper;
	isBBFDDC = timeFDDC > fFDCBBlower && timeFDDC < fFDCBBupper;

	// check FT0 signal
	bool hasNoFT0 = true;
	bool isBB = isBBFT0A || isBBFT0C;
	bool isBG = isBGFT0A || isBGFT0C;
	hasNoFT0 = !isBB && !isBG;
	// check FV0 signal
	bool hasNoFV0A = true;
	isBB = isBBFV0A;
	isBG = isBGFV0A;
	hasNoFV0A = !isBB && !isBG;
	// check FDD signal
	bool hasNoFDD = true;
	isBB = isBBFDDA || isBBFDDC;
	isBG = isBGFDDA || isBGFDDC;
	hasNoFDD = !isBB && !isBG;

	if (hasNoFT0 && hasNoFV0A && hasNoFDD) return true;
	else return false;
}

float invariantMass(float E, float px, float py, float pz)
// Just a simple function to return invariant mass
{
	float im = E*E - px*px - py*py - pz*py;
	return TMath::Sqrt(im);
}

template <typename C, typename MCs, typename P>
void printCollisionData(C collision, MCs generatedCollisions, P slice)
// Function to print collision info
{

	// sliced TF generated collisions according to BC id of the current reconstructed collision
	auto slicedGeneratedCollisions = generatedCollisions.sliceBy(slice,collision.bcId());

	LOGF(info,"Reconstructed collision idx: %d; Associated generated collision idx: %d",
			 collision.globalIndex(),collision.mcCollision().globalIndex());
	LOGF(info,"%d generated collisions in bunch crossing %d related to this reconstructed collisions",
	     slicedGeneratedCollisions.size(),collision.bcId());
}

template <typename T>
void printTrackData(T track)
// Function to print basic info on track and its associated mc particle
{
	if (track.has_mcParticle()) {
		auto mcparticle = track.mcParticle();
		LOGF(info,"Track idx %d, vtx contributor %d, hasITS %d, hasTPC %d, hasTOF %d;"
							" Associated MC particle idx %d, primary %d, PDG code %d",
		     track.globalIndex(),track.isPVContributor(),track.hasITS(),track.hasTPC(),track.hasTOF(),
		     mcparticle.globalIndex(),mcparticle.isPhysicalPrimary(), mcparticle.pdgCode());
	}
}

template <typename Ts>
void printCollisionTracksData(Ts tracks, int setOfCuts)
// Function to loop over tracks associated to a collision and print basic info
{
	int countNoMCparticle = 0;
	for (auto& track : tracks){
		if (!selectTrack(track,setOfCuts)) continue;
		if (track.has_mcParticle()) {
			printTrackData(track);
		}
		else countNoMCparticle++;
	}
	if (countNoMCparticle > 0) LOGF(warning,"This collision has %d tracks without associated mc particle",countNoMCparticle);
}

template <typename MPs, typename P, typename C>
void printCollisionGeneratedParticles(MPs particles, P slice, C collision){
	auto slicedParticles = particles.sliceBy(slice, collision.mcCollision().globalIndex());
	for (auto& slicedParticle : slicedParticles){
		LOGF(info,"Particle idx %d, primary %d",slicedParticle.globalIndex(),slicedParticle.isPhysicalPrimary());
	}
}

void printLargeMessage(std::string info)
// Helper to printf info message to terminal
{
	LOGF(info,"################################### %s ###################################",info);
}


#endif //ALISW_UPCMONTECARLOCENTRALBARRELHELPER_H
