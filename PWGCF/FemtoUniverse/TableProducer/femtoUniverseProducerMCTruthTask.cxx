// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseProducerMCTruthTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;

} // namespace o2::aod

/// \todo fix how to pass array to setSelection, getRow() passing a different
/// type!
// static constexpr float arrayV0Sel[3][3] = {{100.f, 100.f, 100.f}, {0.2f,
// 0.2f, 0.2f}, {100.f, 100.f, 100.f}}; unsigned int rows = sizeof(arrayV0Sel) /
// sizeof(arrayV0Sel[0]); unsigned int columns = sizeof(arrayV0Sel[0]) /
// sizeof(arrayV0Sel[0][0]);

template <typename T>
int getRowDaughters(int daughID, T const& vecID)
{
  int rowInPrimaryTrackTableDaugh = -1;
  for (size_t i = 0; i < vecID.size(); i++) {
    if (vecID.at(i) == daughID) {
      rowInPrimaryTrackTableDaugh = i;
      break;
    }
  }
  return rowInPrimaryTrackTableDaugh;
}

struct FemtoUniverseProducerMCTruthTask {
  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  // Tables being produced
  Produces<aod::FdCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  // Produces<aod::FDMCLabels> outputPartsMCLabels;
  // Produces<aod::FdMCParticles> outputPartsMC;

  // Analysis configs
  Configurable<bool> confIsRun3{"confIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> confIsMC{"confIsMC", false, "Running on MC; implemented only for Run3"};
  Configurable<bool> confIsForceGRP{"confIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};
  Configurable<std::vector<int>> confPDGCodes{"confPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
  Configurable<bool> confAnalysisWithPID{"confAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles"};

  /// Event cuts
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<float> confCentFT0Min{"confCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> confCentFT0Max{"confCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};
  Configurable<bool> confDoSpher{"confDoSpher", false, "Calculate sphericity. If false sphericity will take value of 2."};

  // Track cuts
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPtLowFilterCut{"confPtLowFilterCut", 0.14, "Lower limit for Pt for the filtering tracks"};   // pT low
    Configurable<float> confPtHighFilterCut{"confPtHighFilterCut", 5.0, "Higher limit for Pt for the filtering tracks"}; // pT high
    Configurable<float> confEtaFilterCut{"confEtaFilterCut", 0.8, "Eta cut for the filtering tracks"};
  } ConfFilteringTracks;

  // D0/D0bar cuts
  Configurable<float> yD0CandGenMax{"yD0CandGenMax", 0.5, "Rapidity cut for the D0/D0bar mesons"};

  FemtoUniverseCollisionSelection colCuts;
  FemtoUniverseTrackSelection trackCuts;
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  void init(InitContext&)
  {
    if ((doprocessTrackMC) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }

    colCuts.setCuts(confEvtZvtx, confEvtTriggerCheck, confEvtTriggerSel, confEvtOfflineCheck, confIsRun3, confCentFT0Min, confCentFT0Max);

    colCuts.init(&qaRegistry);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::CutContainerType>(&qaRegistry);

    mRunNumber = 0;
    mMagField = 0.0;

    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  template <typename CollisionType, typename TrackType>
  void fillCollisions(CollisionType const& col, TrackType const& tracks)
  {
    for (const auto& c : col) {
      const auto vtxZ = c.posZ();
      float mult = confIsRun3 ? c.multFV0M() : 0.5 * (c.multFV0M());
      int multNtr = confIsRun3 ? c.multNTracksPV() : c.multTracklets();

      // Removing collisions with Zvtx > 10 cm
      if (std::abs(vtxZ) > confEvtZvtx) {
        continue;
      }
      // colCuts.fillQA(c); //for now, TODO: create a configurable so in the FemroUniverseCollisionSelection.h there is an option to plot QA just for the posZ
      if (confDoSpher) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, mult, multNtr, 2, mMagField);
      }
    }
  }

  template <typename TrackType>
  void fillParticles(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children

    for (const auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (particle.pt() < ConfFilteringTracks.confPtLowFilterCut || particle.pt() > ConfFilteringTracks.confPtHighFilterCut)
        continue;

      int pdgCode = particle.pdgCode();

      if (confAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = confPDGCodes; // necessary due to some features of the Configurable
        for (auto const& pdg : tmpPDGCodes) {
          if (pdgCode == Pdg::kPhi) { // phi meson
            pass = true;
          } else if (std::abs(pdgCode) == Pdg::kD0) { // D0(bar) meson
            pass = true;
          } else if (pdgCode == Pdg::kDPlus) { // D+ meson
            pass = true;
          } else if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (particle.isPhysicalPrimary())
              pass = true;
          }
        }
        if (!pass)
          continue;
      }

      // check if D0/D0bar mesons pass the rapidity cut
      // if pass then saving the orgin of D0/D0bar
      // check if tracks (besides D0/D0bar) pass pseudorapidity cut
      int8_t origin = -99;
      if (std::abs(particle.pdgCode()) == Pdg::kD0) {
        if (std::abs(particle.y()) > yD0CandGenMax) {
          continue;
        } else {
          origin = RecoDecay::getCharmHadronOrigin(tracks, particle);
        }
      } else {
        if (std::abs(particle.eta()) > ConfFilteringTracks.confEtaFilterCut) {
          continue;
        } else {
          origin = -99;
        }
      }

      /// check if we end-up with the correct final state using MC info
      int8_t sign = 0;
      if (std::abs(pdgCode) == Pdg::kD0 && !RecoDecay::isMatchedMCGen(tracks, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
        /// check if we have D0(bar) → π± K∓
        continue;
      }

      // we cannot use isSelectedMinimal since it takes Ncls
      // if (!trackCuts.isSelectedMinimal(track)) {
      //   continue;
      // }

      // trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
      //                  aod::femtouniverseparticle::TrackType::kNoChild>(track);
      //  the bit-wise container of the systematic variations is obtained
      // auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(track);
      // instead of the bitmask, the PDG of the particle is stored as uint32_t

      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  particle.pt(),
                  particle.eta(),
                  particle.phi(),
                  aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                  0,
                  pdgCode,
                  pdgCode,
                  childIDs,
                  origin,
                  -999.);
    }
  }

  void
    processTrackMC(aod::McCollision const&,
                   soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> const& collisions,
                   aod::McParticles const& mcParticles,
                   aod::BCsWithTimestamps const&)
  {
    // magnetic field for run not needed for mc truth
    // fill the tables
    fillCollisions(collisions, mcParticles);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(FemtoUniverseProducerMCTruthTask, processTrackMC, "Provide MC data for track analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUniverseProducerMCTruthTask>(cfgc)};
  return workflow;
}
