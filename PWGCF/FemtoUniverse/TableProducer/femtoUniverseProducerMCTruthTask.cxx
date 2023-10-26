// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
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
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePhiSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Math/Vector4D.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "TMath.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;

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

struct femtoUniverseProducerMCTruthTask {
  int evCount = 0;
  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  // Tables being produced
  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;

  // Analysis configs
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"}; // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};
  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", true, "Activate filling of V0 into femtouniverse tables"};
  Configurable<bool> ConfIsActivatePhi{"ConfIsActivatePhi", true, "Activate filling of Phi into femtouniverse tables"};
  Configurable<std::vector<int>> ConfPDGCodes{"ConfPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
  Configurable<bool> ConfAnalysisWithPID{"ConfAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles"};

  /// Event cuts
  Configurable<bool> ConfEvtUseTPCmult{"ConfEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  // Track cuts
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfPtLowFilterCut{"ConfPtLowFilterCut", 0.14, "Lower limit for Pt for the filtering tracks"};   // pT low
    Configurable<float> ConfPtHighFilterCut{"ConfPtHighFilterCut", 5.0, "Higher limit for Pt for the filtering tracks"}; // pT high
    Configurable<float> ConfEtaFilterCut{"ConfEtaFilterCut", 0.8, "Eta cut for the filtering tracks"};
  } ConfFilteringTracks;

  FemtoUniverseCollisionSelection colCuts;
  FemtoUniverseTrackSelection trackCuts;
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  void init(InitContext&)
  {
    if ((doprocessTrackMC) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);

    colCuts.init(&qaRegistry);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);

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
    const auto vtxZ = col.posZ();
    const auto spher = 0; // colCuts.computeSphericity(col, tracks);
    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in
                                     /// FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfEvtUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions

    // CHECK WHAT CUTS SHOULD BE USED FOR MC TRUTH
    //  if (!colCuts.isSelected(col)) {
    //    if (ConfIsTrigger) {
    //      outputCollision(vtxZ, mult, multNtr, spher, mMagField);
    //    }
    //    return;
    //  }

    colCuts.fillQA(col);
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
  }

  template <typename TrackType>
  void fillParticles(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children

    for (auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (particle.eta() < -ConfFilteringTracks.ConfEtaFilterCut || particle.eta() > ConfFilteringTracks.ConfEtaFilterCut)
        continue;
      if (particle.pt() < ConfFilteringTracks.ConfPtLowFilterCut || particle.pt() > ConfFilteringTracks.ConfPtHighFilterCut)
        continue;

      uint32_t pdgCode = (uint32_t)particle.pdgCode();

      if (ConfAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfPDGCodes; // necessary due to some features of the Configurable
        for (uint32_t pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (pdgCode == 333) { // ATTENTION: workaround for now, because all Phi mesons are NOT primary particles for now.
              pass = true;
            } else {
              if (particle.isPhysicalPrimary())
                pass = true;
            }
          }
        }
        if (!pass)
          continue;
      }

      // we cannot use isSelectedMinimal since it takes Ncls
      // if (!trackCuts.isSelectedMinimal(track)) {
      //   continue;
      // }

      // trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
      //                  aod::femtouniverseparticle::TrackType::kNoChild>(track);
      //  the bit-wise container of the systematic variations is obtained
      // auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(track);
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
                  0,
                  0);
    }
  }

  void
    processTrackMC(aod::FemtoFullCollisionMC const& col,
                   aod::BCsWithTimestamps const&,
                   aod::McCollisions const& mcCollisions,
                   aod::McParticles const& mcParticles)
  {
    if (evCount > 30)
      return;
    evCount++;
    // magnetic field for run not needed for mc truth
    // fill the tables
    fillCollisions(col, mcParticles);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(femtoUniverseProducerMCTruthTask, processTrackMC, "Provide MC data for track analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerMCTruthTask>(cfgc)};
  return workflow;
}
