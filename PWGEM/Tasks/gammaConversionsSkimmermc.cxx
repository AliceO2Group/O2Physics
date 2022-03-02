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

/// \brief write relevant information for photon conversion analysis to an AO2D.root file. This file is then the only necessary input to perform
/// pcm analysis.
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "gammaTables.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct SkimmerMc {

  Produces<aod::GammaConversionTracks> fFuncTableGammaTracks;
  Produces<aod::GammaConversionsInfoTrue> fFuncTableV0InfoTrue;

  // ============================ FUNCTION DEFINITIONS ====================================================
  void process(aod::Collisions::iterator const& theCollision,
               aod::V0s,
               aod::V0Datas const& theV0s,
               tracksAndTPCInfoMC const& theTracks,
               aod::McParticles const& theMcParticles)
  {
    auto fillTrackTable = [&](auto& theV0, auto& theTrack, bool theIsPositive, bool theIsFromConversionPhoton) {
      fFuncTableGammaTracks(
        theV0.v0Id(),
        theIsPositive,
        theIsFromConversionPhoton,
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.eta(),
        theTrack.p(),
        theTrack.phi(),
        theTrack.pt(),
        theTrack.dcaXY(),
        theTrack.tpcNClsCrossedRows(),
        theTrack.tpcSignal(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi());
    };

    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter

      bool lIsConversionPhoton = isConversionPhoton(lV0,
                                                    lTrackPos,
                                                    lTrackNeg,
                                                    theMcParticles);

      fillTrackTable(lV0, lTrackPos, true, lIsConversionPhoton);
      fillTrackTable(lV0, lTrackNeg, false, lIsConversionPhoton);
    }
  }

  // SFS todo: make pretty and short
  template <typename TV0, typename TTRACK, typename TMC>
  bool isConversionPhoton(const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg, const TMC& theMcParticles)
  {
    bool result = false;
    // todo: verify it is enough to check only mother0 being equal

    /* example from https://aliceo2group.github.io/analysis-framework/docs/tutorials/indexTables.html?highlight=_as:
     * track0.collision_as<myCol>().mult() : access multiplicity of collission associated with track0
     */

    auto lMcPos = theTrackPos.template mcParticle_as<aod::McParticles_001>();
    auto lMcNeg = theTrackNeg.template mcParticle_as<aod::McParticles_001>();

    //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
    std::vector<int> lMothers;
    // SFS todo: remove all those mothers_as<aod::McParticles_001>, are the loops even necesarry?
    for (auto& mP : lMcPos.template mothers_as<aod::McParticles_001>()) {
      LOGF(info, "   mother index mP: %d", mP.globalIndex());
      lMothers.push_back(mP.globalIndex());
    }

    if (lMothers.size() > 0) {
      for (auto& mN : lMcNeg.template mothers_as<aod::McParticles_001>()) {
        LOGF(info, "   mother index mN: %d", mN.globalIndex());
        lMothers.push_back(mN.globalIndex());
      }
    }

    // SFS verify theyre all the same category and remove
    int lSame = 0;
    {
      if (lMothers.size() == 2) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size2: 01");
          lSame = 1;
        }
      }

      if (lMothers.size() == 3) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size3: 01");
          lSame = 2;
        }
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size2: 02");
          lSame = 3;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size2: 12");
          lSame = 4;
        }
      }

      if (lMothers.size() == 4) {
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size4 02");
          lSame = 4;
        }
        if (lMothers[1] == lMothers[3]) {
          LOGF(info, "size4 13");
          lSame = 5;
        }
        if (lMothers[0] == lMothers[3]) {
          LOGF(info, "size4 03");
          lSame = 6;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size4 12");
          lSame = 7;
        }
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size4 01");
          lSame = 8;
        }
        if (lMothers[2] == lMothers[3]) {
          LOGF(info, "size4 23");
          lSame = 9;
        }
      }
    }

    if (lSame) {
      // SFS todo: actually no loop required here, for this
      for (auto& lMother : lMcNeg.template mothers_as<aod::McParticles_001>()) {

        if ((result = lMother.pdgCode() == 22)) {
          fFuncTableV0InfoTrue(
            theV0.v0Id(),
            lMcPos.vx(), lMcPos.vy(), lMcPos.vz(),
            lMcPos.px(), lMcPos.py(), lMcPos.pz(),
            lMcNeg.px(), lMcNeg.py(), lMcNeg.pz(),
            lMother.px(), lMother.py(), lMother.pz(),
            lMother.eta(),
            lMother.p(),
            lMother.phi(),
            lMother.pt());
        }
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SkimmerMc>(cfgc)};
}
