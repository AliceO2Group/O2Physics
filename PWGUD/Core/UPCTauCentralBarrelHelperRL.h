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

#ifndef ALISW_UPCTAUCENTRALBARRELHELPERRL_H
#define ALISW_UPCTAUCENTRALBARRELHELPERRL_H

#include <TRandom.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum MyParticle {
	P_ELECTRON = 0,
	P_MUON = 1,
	P_PION = 2,
	P_KAON = 3,
	P_PROTON = 4
};

void printLargeMessage(std::string info)
// Helper to printf info message to terminal
{
	LOGF(info,"################################### %s ###################################",info);
}

void printMediumMessage(std::string info)
// Helper to printf info message to terminal
{
	LOGF(info,"+++++++++++++ %s +++++++++++++",info);
}

template <typename T>
int testPIDhypothesis(T trackPIDinfo, float nSigmaShift = 0., bool isMC = false)
// Choose, which particle it is according to PID
{
	float nSigmaTPC[5];
	nSigmaTPC[P_ELECTRON] = abs(trackPIDinfo.tpcNSigmaEl());
	nSigmaTPC[P_MUON] = abs(trackPIDinfo.tpcNSigmaMu());
	nSigmaTPC[P_PION] = abs(trackPIDinfo.tpcNSigmaPi());
	nSigmaTPC[P_KAON] = abs(trackPIDinfo.tpcNSigmaKa());
	nSigmaTPC[P_PROTON] = abs(trackPIDinfo.tpcNSigmaPr());
    // Correction if TPC tuneOnData is wrong
    if (isMC) {
        for (int i=0; i<5; i++) nSigmaTPC[i] -= nSigmaShift;
    }
	int enumChoiceTPC = std::distance(std::begin(nSigmaTPC),
	                                  std::min_element(std::begin(nSigmaTPC), std::end(nSigmaTPC)));

//	if (trackPIDinfo.hasTPC()) {
//		return enumChoiceTPC;
//	} else {
//		LOGF(warning, "testPIDhypothesis failed - track did not leave information in TPC");
//		return -1;
//	}

//	if (enumChoiceTPC == P_MUON) {
//	if (enumChoiceTPC == P_ELECTRON || enumChoiceTPC == P_MUON) {
//		return P_PION;
//	if (enumChoiceTPC == P_PION || enumChoiceTPC == P_MUON) {
//		return P_ELECTRON;
//	} else {
//		return enumChoiceTPC;
//	}

  float nSigmaTOF[5];
  nSigmaTOF[P_ELECTRON] = trackPIDinfo.tofNSigmaEl();
  nSigmaTOF[P_MUON] = trackPIDinfo.tofNSigmaMu();
  nSigmaTOF[P_PION] = trackPIDinfo.tofNSigmaPi();
  nSigmaTOF[P_KAON] = trackPIDinfo.tofNSigmaKa();
  nSigmaTOF[P_PROTON] = trackPIDinfo.tofNSigmaPr();
  int enumChoiceTOF = std::distance(std::begin(nSigmaTOF),
                                    std::min_element(std::begin(nSigmaTOF), std::end(nSigmaTOF)));

	if (trackPIDinfo.hasTPC() || trackPIDinfo.hasTOF()) {
		if (trackPIDinfo.hasTOF()) {
			return enumChoiceTOF;
		} else {
			return enumChoiceTPC;
		}
	} else {
		LOGF(debug, "testPIDhypothesis failed - track did not leave information in TPC or TOF");
		return -1;
	}

//	float nSigmaTPCpick[3] = {nSigmaTPC[P_ELECTRON], nSigmaTPC[P_MUON], nSigmaTPC[P_PION]};
//	int enumChoiceTPCpick = std::distance(std::begin(nSigmaTPCpick),
//	                                      std::min_element(std::begin(nSigmaTPCpick), std::end(nSigmaTPCpick)));
//
////   do PID using TPC+TOF
//  if (trackPIDinfo.hasTPC() && trackPIDinfo.hasTOF()) {
//    if (enumChoiceTPC == P_ELECTRON || enumChoiceTPC == P_MUON || enumChoiceTPC == P_PION) {
//      if (enumChoiceTOF == P_KAON) {
//        return P_KAON; // probably kaon
//      } else if (enumChoiceTOF == P_PROTON) {
//        return P_PROTON; // probably proton
//      } else {
//        if (enumChoiceTPC == P_ELECTRON) {
//          return P_ELECTRON; // probably electron
//        } else if (enumChoiceTPC == P_MUON) {
//          return P_MUON; // probably muon
//        } else {
//          return P_PION; // probably pion
//        }
//      }
//    } else {
//      if (enumChoiceTOF == P_KAON) {
//        return P_KAON; // probably kaon
//      } else if (enumChoiceTOF == P_PROTON) {
//        return P_PROTON; // probably proton
//      } else {
//        if (enumChoiceTPCpick == P_ELECTRON) {
//          return P_ELECTRON; // probably misidentified electron
//        } else if (enumChoiceTPCpick == P_MUON) {
//          return P_MUON; // probably misidentified muon
//        } else {
//          return P_PION; // probably misidentified pion
//        }
//      }
//    }
//    // do PID using TPC only
//  } else if (trackPIDinfo.hasTPC()) {
//    return enumChoiceTPC;
//    // do PID using TOF only
//  } else if (trackPIDinfo.hasTOF()) {
//    return enumChoiceTOF;
//    // give warning and return non-sense
//  } else {
//    LOGF(warning, "testPIDhypothesis failed - track did not leave information in TPC or TOF");
//    return -1;
//  }

}

template <typename T>
int trackPDG(T trackPIDinfo)
// using testPIDhypothesis, reads enumMyParticle and return pdg value
{
	if (testPIDhypothesis(trackPIDinfo) == P_ELECTRON) return 11;
	else if (testPIDhypothesis(trackPIDinfo) == P_MUON) return 13;
	else if (testPIDhypothesis(trackPIDinfo) == P_PION) return 211;
	else if (testPIDhypothesis(trackPIDinfo) == P_KAON) return 321;
	else if (testPIDhypothesis(trackPIDinfo) == P_PROTON) return 2212;
	else {
		printMediumMessage("Something is wrong with track PDG selector");
		return -1.;
	}
}

int enumMyParticle(int valuePDG)
// reads pdg value and returns particle number as in enumMyParticle
{
	if (abs(valuePDG) == 11) return P_ELECTRON;
	else if (abs(valuePDG) == 13) return P_MUON;
	else if (abs(valuePDG) == 211) return P_PION;
	else if (abs(valuePDG) == 321) return P_KAON;
	else if (abs(valuePDG) == 2212) return P_PROTON;
	else {
		printMediumMessage("PDG value not found in enumMyParticle. Returning -1.");
		return -1.;
	}
}

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
//	if (selection==1 && track.trackType()!=1) return false;
	// Is central barrel vertex contributor
	if (selection==2 && track.isPVContributor()!=1) return false;
	// Is central barrel track selection global track
//	if (selection==3 && track.isQualityTrack()!=1) return false;
	// Is central barrel track selection global track
//	if (selection==4 && track.isPrimaryTrack()!=1) return false;
	// Is central barrel track selection global track
//	if (selection==5 && track.isGlobalTrack()!=1) return false;
	// Is central barrel track selection has TOF
	if (selection==6 && track.hasTOF()!=1) return false;
	// Is central barrel track selection nSigmaElectron cut
	if (selection==7 && abs(track.tpcNSigmaEl())>3.) return false;

	return true;
}

template <typename T>
bool selectTrack(T track, int setOfCuts)
// Do selection of reconstructed track
{
	if (setOfCuts<1) return true;
	// Is central barrel propagated track
	if (setOfCuts<2 && trackSelection(track,1)!=1) return false;
	// Is central barrel vertex contributor
	if (setOfCuts<3 && trackSelection(track,2)!=1) return false;
	// Is central barrel track selection global track
	if (setOfCuts<4 && trackSelection(track,3)!=1) return false;
	// Is central barrel track selection global track
	if (setOfCuts<5 && trackSelection(track,4)!=1) return false;
	// Is central barrel track selection global track
	if (setOfCuts<6 && trackSelection(track,5)!=1) return false;
	// Is central barrel track selection global track has not TOF
	if (setOfCuts<7 && trackSelection(track,6)==1) return false;

	return true;
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


float momentum(float px, float py, float pz)
// Just a simple function to return momentum
{
	return TMath::Sqrt(px * px + py * py + pz * pz);
}

float invariantMass(float E, float px, float py, float pz)
// Just a simple function to return invariant mass
{
	return TMath::Sqrt(E * E - px * px - py * py - pz * pz);
}

float phi(float px, float py)
// Just a simple function to return azimuthal angle
{
	if (px != 0) return TMath::ATan(py / px);
	return -999.;
}

float eta(float px, float py, float pz)
// Just a simple function to return pseudorapidity
{
	float eta = -999.;
	float mom = momentum(px, py, pz);
	if (mom != 0) eta = TMath::ATanH(pz / mom);
	if (-1. < eta && eta < 1.) return eta;
	return -999.;
}

float energy(float mass, float px, float py, float pz)
// Just a simple function to return track energy
{
	return TMath::Sqrt(mass * mass + px * px + py * py + pz * pz);
}

float rapidity(float mass, float px, float py, float pz)
// Just a simple function to return track rapidity
{
	return 0.5 * TMath::Log((energy(mass, px, py, pz) + pz) / (energy(mass, px, py, pz) - pz));
}

double calculateAcoplanarity(double phi_trk1, double phi_trk2)
// Function to calculate acoplanarity of two tracks based on phi of both tracks, which is in interval (0,2*pi)
{
    double aco = TMath::Abs(phi_trk1-phi_trk2);
    if (aco <= o2::constants::math::PI) return aco;
    else return (o2::constants::math::TwoPI-aco);
}

template <typename Ts>
int countPVcontributors(Ts tracks)
// Function to loop over tracks associated to a collision and return total of PV contributors
{
	int nTotal = 0;
	for (auto& track : tracks){
		if (!trackSelection(track,2)) continue;
		nTotal++;
	}
	return nTotal;
}

template <typename Ts>
int countGlobalTracks(Ts tracks)
// Function to loop over tracks associated to a collision and return total of PV contributors
{
	int nTotal = 0;
	for (auto& track : tracks){
		if (!trackSelection(track,5)) continue;
		nTotal++;
	}
	return nTotal;
}

template <typename TTracks>
void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
{
    for (const auto& tr : tracks) {
        int32_t candId = tr.udCollisionId();
        if (candId < 0) {
            continue;
        }
        tracksPerCand[candId].push_back(tr.globalIndex());
    }
}

#endif //ALISW_UPCTAUCENTRALBARRELHELPERRL_H
