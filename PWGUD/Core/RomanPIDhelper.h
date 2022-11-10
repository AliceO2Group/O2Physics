//
// Created by Roman Lavička on 10.11.2022.
//

#ifndef ALISW_ROMANPIDHELPER_H
#define ALISW_ROMANPIDHELPER_H

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
}


#endif //ALISW_ROMANPIDHELPER_H
