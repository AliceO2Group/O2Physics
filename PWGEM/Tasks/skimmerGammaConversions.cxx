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
using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi>;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct skimmerGammaConversions {

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionZ_MCRec", "hCollisionZ_MCRec", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_all_MCTrue", "hCollisionZ_all_MCTrue", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_MCTrue", "hCollisionZ_MCTrue", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hMcParticlesSize", "hMcParticlesSize", {HistType::kTH1F, {{100, 0.f, 1000000.f}}}},
      {"hMotherSameNess", "hMotherSameNess", {HistType::kTH1F, {{13, 0.f, 14.f}}}},
    },
  };

  Produces<aod::V0DaughterTracks> fFuncTableV0DaughterTracks;
  Produces<aod::McGammasTrue> fFuncTableMcGammasFromConfirmedV0s;

  // ============================ FUNCTION DEFINITIONS ====================================================
  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::V0s,
                  aod::V0Datas const& theV0s,
                  tracksAndTPCInfo const& theTracks)
  {
    auto fillTrackTable = [&](auto& theV0, auto& theTrack, bool theIsPositive, bool theIsFromConversionPhoton) {
      fFuncTableV0DaughterTracks(
        theV0.v0Id(),
        theIsFromConversionPhoton,
        theTrack.dcaXY(),
        theTrack.eta(),
        theTrack.p(),
        theTrack.phi(),
        theTrack.pt(),
        theIsPositive,
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcNClsCrossedRows(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi(),
        theTrack.tpcSignal());
    };

    registry.fill(HIST("hCollisionZ"), theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfo>(); // negative daughter

      fillTrackTable(lV0, lTrackPos, true, false);
      fillTrackTable(lV0, lTrackNeg, false, false);
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processRec, "process reconstructed info only", true);

  void processMc(aod::McCollision const& theMcCollision,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                            aod::Collisions>> const& theCollisions,
                 aod::V0s,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& theMcParticles)
  {
    auto fillTrackTable = [&](auto& theV0, auto& theTrack, bool theIsPositive, bool theIsFromConversionPhoton) {
      fFuncTableV0DaughterTracks(
        theV0.v0Id(),
        theIsFromConversionPhoton,
        theTrack.dcaXY(),
        theTrack.eta(),
        theTrack.p(),
        theTrack.phi(),
        theTrack.pt(),
        theIsPositive,
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcNClsCrossedRows(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi(),
        theTrack.tpcSignal());
    };

    registry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());
    if (theCollisions.size() == 0) {
      return;
    }

    registry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());
    registry.fill(HIST("hMcParticlesSize"), theMcParticles.size());

    for (auto& lCollision : theCollisions) {
      registry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());

      // todo: replace by sliceByCached
      auto lGroupedV0s = theV0s.sliceBy(aod::v0data::collisionId, lCollision.globalIndex());
      for (auto& lV0 : lGroupedV0s) {

        auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
        auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter

        bool lIsConversionPhoton = isConversionPhoton(lV0,
                                                      lTrackPos,
                                                      lTrackNeg);

        fillTrackTable(lV0, lTrackPos, true, lIsConversionPhoton);
        fillTrackTable(lV0, lTrackNeg, false, lIsConversionPhoton);
      }
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processMc, "process reconstructed and mc info ", false);

  // SFS todo: make pretty and short
  template <typename TV0, typename TTRACK>
  bool isConversionPhoton(TV0 const& theV0,
                          TTRACK const& theTrackPos,
                          TTRACK const& theTrackNeg)
  {
    bool result = false;
    // todo: verify it is enough to check only mother0 being equal

    /* example from https://aliceo2group.github.io/analysis-framework/docs/tutorials/indexTables.html?highlight=_as:
     * track0.collision_as<myCol>().mult() : access multiplicity of collission associated with track0
     */

    auto lMcPos = theTrackPos.template mcParticle_as<aod::McParticles>();
    auto lMcNeg = theTrackNeg.template mcParticle_as<aod::McParticles>();

    // get all mc mothers from positive and negative tracks
    //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
    std::vector<int> lMothers;
    // SFS todo: remove all those mothers_as<aod::McParticles_001>, are the loops even necesarry?
    for (auto& mP : lMcPos.template mothers_as<aod::McParticles>()) {
      LOGF(info, "   mother index mP: %d", mP.globalIndex());
      lMothers.push_back(mP.globalIndex());
    }

    if (lMothers.size() > 0) {
      for (auto& mN : lMcNeg.template mothers_as<aod::McParticles>()) {
        LOGF(info, "   mother index mN: %d", mN.globalIndex());
        lMothers.push_back(mN.globalIndex());
      }
    }

    // SFS verify theyre all the same category and remove
    int lMotherSameNess = 0;
    {
      if (lMothers.size() == 2) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size2: 01");
          lMotherSameNess = 1;
        }
      }

      if (lMothers.size() == 3) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size3: 01");
          lMotherSameNess = 2;
        }
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size2: 02");
          lMotherSameNess = 3;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size2: 12");
          lMotherSameNess = 4;
        }
      }

      if (lMothers.size() == 4) {
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size4 02");
          lMotherSameNess = 4;
        }
        if (lMothers[1] == lMothers[3]) {
          LOGF(info, "size4 13");
          lMotherSameNess = 5;
        }
        if (lMothers[0] == lMothers[3]) {
          LOGF(info, "size4 03");
          lMotherSameNess = 6;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size4 12");
          lMotherSameNess = 7;
        }
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size4 01");
          lMotherSameNess = 8;
        }
        if (lMothers[2] == lMothers[3]) {
          LOGF(info, "size4 23");
          lMotherSameNess = 9;
        }
      }
    }

    registry.fill(HIST("hMotherSameNess"), 0.5 + (float)lMotherSameNess);

    // if both tracks have exactly one and the same mother
    if (lMotherSameNess == 1) {
      // SFS todo: actually no loop required here, for this
      for (auto& lMcMother : lMcNeg.template mothers_as<aod::McParticles>()) {
        if ((result = lMcMother.pdgCode() == 22)) {

          auto lDaughters = lMcMother.template daughters_as<aod::McParticles>();
          float lDaughter0Vx = -1.;
          float lDaughter0Vy = -1.;
          float lDaughter0Vz = -1.;
          float lV0Radius = -1.;
          if (lDaughters.size()) {
            auto lDaughter0 = lDaughters.begin();
            lDaughter0Vx = lDaughter0.vx();
            lDaughter0Vy = lDaughter0.vy();
            lDaughter0Vz = lDaughter0.vz();
            lV0Radius = sqrt(pow(lDaughter0Vx, 2) + pow(lDaughter0Vy, 2));
          }

          fFuncTableMcGammasFromConfirmedV0s(
            lMcMother.mcCollisionId(),
            lMcMother.globalIndex(),
            theV0.v0Id(),
            lMcMother.statusCode(),
            lMcMother.flags(),
            lMcMother.px(), lMcMother.py(), lMcMother.pz(),
            lMcMother.vx(), lMcMother.vy(), lMcMother.vz(), lMcMother.vt(),
            lDaughters.size(),
            lMcMother.eta(), lMcMother.phi(), lMcMother.p(), lMcMother.pt(), lMcMother.y(),
            lDaughter0Vx, lDaughter0Vy, lDaughter0Vz,
            lV0Radius);
        }
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversions>(cfgc)};
}
