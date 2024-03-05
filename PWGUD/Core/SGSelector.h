#ifndef PWGUD_CORE_SGSELECTOR_H_
#define PWGUD_CORE_SGSELECTOR_H_

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGCutParHolder.h"

class SGSelector {
public:
    SGSelector() : fPDG(TDatabasePDG::Instance()) {}

    template <typename CC, typename BCs, typename TCs, typename FWs>
    int Print(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks) {
        LOGF(info, "Size of array %i", collision.size());
        return 1;
    }

    template <typename CC, typename BCs>
    //int IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks) {
    int IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange) {
        LOGF(debug, "Collision %f", collision.collisionTime());
        LOGF(debug, "Number of close BCs: %i", bcRange.size());

        bool gA = true, gC = true;
	int bA = 0, bC = 0;
	int64_t bcA[5], bcC[5];
        for (auto const& bc : bcRange) {
	    bcA[bA] = 0;
	    bcC[bC] = 0;
            if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {gA = false; bcA[bA] = bc.globalBC(); bA++;}
            if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {gC = false; bcC[bC] = bc.globalBC(); bC++;}
        }
        if (!gA && !gC) return 3;
	if (bA > 1 || bC > 1){
        LOGF(info, "Number of BCs with FT0 signal: %i, %i", bA, bC);
	for (int i = 0; i < bA; i++){
		LOGF(info, "Busy BC A-side: %i", bcA[i]);
	}
	for (int i = 0; i < bC; i++){
		LOGF(info, "Busy BC C-side: %i", bcC[i]);
	}
	}
        if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
            return 4;
        }
        return gA && gC ? 2 : (gA ? 0 : 1);
}
  template <typename TFwdTrack>
    int FwdTrkSelector(TFwdTrack const& fwdtrack) {
    if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) return 1;
    else return 0;        
    }
  template <typename TTrack>
    int TrkSelector(SGCutParHolder diffCuts, TTrack const& track) {
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }
    auto lvtmp = TLorentzVector();
    if (!track.isPVContributor()) return 1;
    lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
    if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) return 2;
    if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) return 3;
    if (diffCuts.globalTracksOnly() && !track.isGlobalTrack()) return 4;
    if (!diffCuts.ITSOnlyTracks() && !track.hasTPC()) return 5;
    return 0;
    }

private:
    TDatabasePDG* fPDG;
};

#endif // PWGUD_CORE_SGSELECTOR_H_

