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

    template <typename CC, typename BCs, typename TCs, typename FWs>
    int IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks) {
        LOGF(debug, "Collision %f", collision.collisionTime());
        LOGF(debug, "Number of close BCs: %i", bcRange.size());

        bool gA = true, gC = true;
        for (auto const& bc : bcRange) {
            if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) gA = false;
            if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) gC = false;
        }
        if (!gA && !gC) return 3;

        LOGF(debug, "FwdTracks %i", fwdtracks.size());
        for (auto& fwdtrack : fwdtracks) {
            if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) {
                return 4;
            }
        }

        double rgtrwTOF = 0.0;
        for (auto& track : tracks) {
            if (track.isGlobalTrack() && !track.isPVContributor()) {
                return 5;
            }
            if (diffCuts.globalTracksOnly() && track.isPVContributor() && !track.isGlobalTrack()) {
                return 6;
            }
            if (!diffCuts.ITSOnlyTracks() && track.isPVContributor() && !track.hasTPC()) {
                return 7;
            }
            if (track.isPVContributor() && track.hasTOF()) {
                rgtrwTOF += 1.0;
            }
        }
        if (collision.numContrib() > 0) {
            rgtrwTOF /= static_cast<double>(collision.numContrib());
        }
        if (rgtrwTOF < diffCuts.minRgtrwTOF()) {
            return 8;
        }

        if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
            return 9;
        }

// PID, pt, and eta of tracks, invariant mass, and net charge
// consider only vertex tracks

    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      if (track.isPVContributor()) {
        // pt
        lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
        if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
          return 10;
        }
        // eta
        if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
          return 11;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }


        return gA && gC ? 2 : (gA ? 0 : 1);
    }

private:
    TDatabasePDG* fPDG;
};

#endif // PWGUD_CORE_SGSELECTOR_H_

