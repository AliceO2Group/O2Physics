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

#ifndef ALISW_RLHELPER_H
#define ALISW_RLHELPER_H

#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum MyParticle { P_ELECTRON = 0, P_MUON = 1, P_PION = 2, P_KAON = 3, P_PROTON = 4};

template <typename T>
int testPIDhypothesis(T trackPIDinfo)
// Choose, which particle it is according to PID
{
	float nSigmaTPC[5];
	nSigmaTPC[P_ELECTRON] = trackPIDinfo.tpcNSigmaEl();
	nSigmaTPC[P_MUON] = trackPIDinfo.tpcNSigmaMu();
	nSigmaTPC[P_PION] = trackPIDinfo.tpcNSigmaPi();
	nSigmaTPC[P_KAON] = trackPIDinfo.tpcNSigmaKa();
	nSigmaTPC[P_PROTON] = trackPIDinfo.tpcNSigmaPr();
	int enumChoiceTPC = std::distance(std::begin(nSigmaTPC), std::min_element(std::begin(nSigmaTPC), std::end(nSigmaTPC)));
	float nSigmaTPCpick[3] = {nSigmaTPC[P_ELECTRON], nSigmaTPC[P_MUON], nSigmaTPC[P_PION]};
	int enumChoiceTPCpick = std::distance(std::begin(nSigmaTPCpick), std::min_element(std::begin(nSigmaTPCpick), std::end(nSigmaTPCpick)));

	float nSigmaTOF[5];
	nSigmaTOF[P_ELECTRON] = trackPIDinfo.tofNSigmaEl();
	nSigmaTOF[P_MUON] = trackPIDinfo.tofNSigmaMu();
	nSigmaTOF[P_PION] = trackPIDinfo.tofNSigmaPi();
	nSigmaTOF[P_KAON] = trackPIDinfo.tofNSigmaKa();
	nSigmaTOF[P_PROTON] = trackPIDinfo.tofNSigmaPr();
	int enumChoiceTOF = std::distance(std::begin(nSigmaTOF), std::min_element(std::begin(nSigmaTOF), std::end(nSigmaTOF)));

	// do PID using TPC+TOF
	if (trackPIDinfo.hasTPC() && trackPIDinfo.hasTOF()){
		if      (enumChoiceTPC == P_ELECTRON || enumChoiceTPC == P_MUON || enumChoiceTPC == P_PION){
			if      (enumChoiceTOF == P_KAON  ) return P_KAON; // probably kaon
			else if (enumChoiceTOF == P_PROTON) return P_PROTON; // probably proton
			else {
				if      (enumChoiceTPC == P_ELECTRON) return P_ELECTRON; // probably electron
				else if (enumChoiceTPC == P_MUON)     return P_MUON; // probably muon
				else                            return P_PION; // probably pion
			}
		}
		else {
			if      (enumChoiceTOF == P_KAON  ) return P_KAON; // probably kaon
			else if (enumChoiceTOF == P_PROTON) return P_PROTON; // probably proton
			else {
				if      (enumChoiceTPCpick == P_ELECTRON) return P_ELECTRON; // probably misidentified electron
				else if (enumChoiceTPCpick == P_MUON)     return P_MUON; // probably misidentified muon
				else if (enumChoiceTPCpick == P_PION)     return P_PION; // probably misidentified pion
			}
		}
	}
	// do PID using TPC only
	else if (trackPIDinfo.hasTPC()) return enumChoiceTPC;
	// do PID using TOF only
	else if (trackPIDinfo.hasTOF()) return enumChoiceTOF;
	// give warning and return non-sense
	else {
		LOGF(warning,"testPIDhypothesis failed - track did not leave information in TPC or TOF");
		return -1;
	}
	return -1;
}


float momentum(float px, float py, float pz)
// Just a simple function to return momentum
{
	return TMath::Sqrt(px*px+py*py+pz*py);
}

float invariantMass(float E, float px, float py, float pz)
// Just a simple function to return invariant mass
{
	return TMath::Sqrt(E*E-px*px-py*py-pz*py);
}

float phi(float px, float py)
// Just a simple function to return azimuthal angle
{
	if (px!=0) return TMath::ATan(py/px);
	return -999.;
}

float eta(float px, float py, float pz)
// Just a simple function to return pseudorapidity
{
	float eta = -999.;
	float mom = momentum(px,py,pz);
	if (mom!=0) eta = TMath::ATanH(pz/mom);
	if (-1.<eta && eta<1.) return eta;
	return -999.;
}

float energy(float mass, float px, float py, float pz)
// Just a simple function to return track energy
{
	return TMath::Sqrt(mass*mass+px*px+py*py+pz*py);
}

float rapidity(float mass, float px, float py, float pz)
// Just a simple function to return track rapidity
{
	return 0.5*TMath::Log((energy(mass,px,py,pz)+pz)/(energy(mass,px,py,pz)-pz));
}

template <typename E>
int getElectronCharge(E const& generatedElectron)
// Check if particle is electron or positron and return charge accordingly. Return zero if particle is not electron/positron
{
	if (generatedElectron.pdgCode() == 11) return -1;
	else if (generatedElectron.pdgCode() == -11) return 1;
	else return 0;
}

template <typename Es>
int64_t getEvSelsIndexOfThisCollisionBC(Es const& infoEvSels, uint64_t globalBC)
// reads full event selection table end return global index corresponding to given global BC. Return -1 when fails.
{
	for (auto& infoEvSel : infoEvSels){
		if (infoEvSel.bcId() == globalBC) return infoEvSel.globalIndex();
	}
	return -1;
}

template <typename C>
bool isEvSelFITempty(C const& collision)
// Get FIT information from EventSelection task for each collision
{
	if (collision.has_foundFT0() || collision.has_foundFV0() || collision.has_foundFDD()) return false;
	else return true;
}

template <typename T>
bool isFITempty(T const& FITinfo)
// Return true if FIT had no signal
{
	// check FT0 signal
	bool hasNoFT0 = true;
	bool isBB = FITinfo.bbFT0A() || FITinfo.bbFT0C();
	bool isBG = FITinfo.bgFT0A() || FITinfo.bgFT0C();
	hasNoFT0 = !isBB && !isBG;
	// check FV0 signal
	bool hasNoFV0A = true;
	isBB = FITinfo.bbFV0A();
	isBG = FITinfo.bgFV0A();
	hasNoFV0A = !isBB && !isBG;
	// check FDD signal
	bool hasNoFDD = true;
	isBB = FITinfo.bbFDDA() || FITinfo.bbFDDC();
	isBG = FITinfo.bgFDDA() || FITinfo.bgFDDC();
	hasNoFDD = !isBB && !isBG;

	if (hasNoFT0 && hasNoFV0A && hasNoFDD) return true;
	else return false;
}

template <typename T>
bool isUDprimaryTrack(T const& udtrack)
// TrackSelection::kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz;
{
	// temporary hardcoded input
	float maxDcaXY = 0.5;// in mm
	float maxDcaZ = 2.;// in mm
	// selection GoldenChi2
	// seems to be always true for Run 3... at least it is not coded
	// selection DCAxy
	if (abs(udtrack.dcaXY()) > maxDcaXY) return false;
	// selection DCAz
	if (abs(udtrack.dcaZ()) > maxDcaZ) return false;
	// passed all selections
	return true;
}

template <typename T>
bool isUDinAcceptanceTrack(T const& udtrack)
// TrackSelection::kInAcceptanceTracks = kPtRange | kEtaRange;
{
	// temporary hardcoded input
	float minPt = 0.;// in GeV
	float minEta = -10.;
	float maxPt = 10.;// in GeV
	float maxEta = 10.;
	float currentEta = eta(udtrack.px(),udtrack.py(),udtrack.pz());
	// selection pt
	if (udtrack.pt() < minPt || udtrack.pt() > maxPt) return false;
	// selection eta
	if (currentEta < minEta || currentEta > maxEta) return false;
	// passed all selections
	return true;
}

template <typename T>
bool isUDqualityTrack(T const& udtrack)
// TrackSelection::kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits;
{
	// temporary hardcoded input
	float cutITSchi2ndf = 6.;
	uint8_t cutNitsClusters = 4;
	float cutTPCchi2ndf = 4.;
	int16_t cutNtpcClusters = 0;
	int16_t cutNtpcCrossedRows = 70;
	float cutCrossedRowsOverNclusters = 0.8;
	// track type
	// TODO //ignoring for the moment (its either innermost update track (0) or propagated track (1), the rest is Run 2).
	// ITS hits
	// TODO //ignoring for the moment (some crazy function which I haven't found anywhere to be used)
	// ITS refit
	if (!udtrack.hasITS()) return false;
	// ITS chi2/ndf
	if (udtrack.itsChi2NCl() > cutITSchi2ndf) return false;
	// ITS n clusters
	if (udtrack.itsNCls() < cutNitsClusters) return false;
	// TPC refit
	if (!udtrack.hasTPC()) return false;
	// TPC chi2/ndf
	if (udtrack.tpcChi2NCl() > cutTPCchi2ndf) return false;
	// TPC n clusters
	int16_t nFoundClusters = udtrack.tpcNClsFindable() - udtrack.tpcNClsFindableMinusFound();
	if (nFoundClusters < cutNtpcClusters) return false;
	// TPC crossed rows
	if (udtrack.tpcNClsCrossedRows() < cutNtpcCrossedRows) return false;
	// TPC crossed rows over n clusters
	float crossedRowsOverNclusters = (float)udtrack.tpcNClsCrossedRows()/nFoundClusters;
	if (crossedRowsOverNclusters < cutCrossedRowsOverNclusters) return false;
	// passed all selections
	return true;
}

template <typename T>
bool isUDglobalTrack(T const& udtrack)
// combine quality+primary+acceptance
{
	if (!isUDinAcceptanceTrack(udtrack)) return false;
	if (!isUDprimaryTrack(udtrack)) return false;
	if (!isUDqualityTrack(udtrack)) return false;
	return true;
}

template <typename T>
bool trackSelection(T const& udtrack, int selection)
// Do selection of reconstructed track
{

	if (selection==0) return true;
	// Is central barrel propagated track
	// TODO //if (selection==1 && udtrack.trackType()!=1) return false;
	// Is central barrel vertex contributor
	if (selection==2 && udtrack.isPVContributor()!=1) return false;
	// Is central barrel track selection global track
	if (selection==3 && isUDqualityTrack(udtrack)!=1) return false;
	// Is central barrel track selection global track
	if (selection==4 && isUDprimaryTrack(udtrack)!=1) return false;
	// Is central barrel track selection global track
	if (selection==5 && isUDglobalTrack(udtrack)!=1) return false;

	return true;
}

template <typename T>
bool selectTrack(T const& track, int setOfCuts)
// Do selection of reconstructed track
{
	if (setOfCuts<1) return true;
	// Is central barrel propagated track)
	if (setOfCuts<2 && trackSelection(track,1)!=1) return false;
	// Is central barrel vertex contributor
	if (setOfCuts<3 && trackSelection(track,2)!=1) return false;

	return true;
}

template < typename C>
void fillEventSelectionHistogram(HistogramRegistry &registry, C const& collision)
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
void fillTrackSelectionHistogram(HistogramRegistry &registry, Ts const& tracks, C const& collision)
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
		if (isFITempty(collision)) registry.get<TH1>(HIST("hEffectOfTrackSelections"))->Fill(17);
		registry.get<TH1>(HIST("hEffectOfTrackSelections"))->GetXaxis()->ChangeLabel(18,70,0.02,30,-1,-1,"FIT empty");

	}
}

template <typename T>
void printTrackData(T const& track)
// Function to print basic info on track
{
	LOGF(info,"Track idx %d, vtx contributor %d, hasITS %d, hasTPC %d, hasTOF %d",
			 track.globalIndex(),track.isPVContributor(),track.hasITS(),track.hasTPC(),track.hasTOF());
}

template <typename T>
void printTrackParticleData(T const& track)
// Function to print basic info on track and its associated mc particle
{
	printTrackData(track);
	if (track.has_udMcParticle()) {
		auto mcparticle = track.udMcParticle();
		LOGF(info," Associated MC particle idx %d, primary %d, PDG code %d",
				 mcparticle.globalIndex(),mcparticle.isPhysicalPrimary(), mcparticle.pdgCode());
	}
}

template <typename Ts>
void printCollisionTracksData(Ts const& tracks, int setOfCuts)
// Function to loop over tracks associated to a collision and print basic info
{
	int countNoMCparticle = 0;
	for (auto& track : tracks){
		if (!selectTrack(track,setOfCuts)) continue;
		if (track.has_udMcParticle()) {
			printTrackParticleData(track);
		}
		else countNoMCparticle++;
	}
	if (countNoMCparticle > 0) LOGF(warning,"This collision has %d tracks without associated mc particle",countNoMCparticle);
}

template <typename MPs, typename P, typename C>
void printCollisionGeneratedParticles(MPs const& particles, P const& slice, C const& collision){
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


#endif //ALISW_RLHELPER_H
