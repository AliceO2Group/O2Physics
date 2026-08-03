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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <set>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::rctsel;

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

  float mMagField = 0.;
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; /// Accessing the CCDB

  // Tables being produced
  Produces<aod::FdCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  // Produces<aod::FDMCLabels> outputPartsMCLabels;
  // Produces<aod::FdMCParticles> outputPartsMC;

  RCTFlagsChecker rctChecker;

  static constexpr int kHadElasticScatt = 20;   // kPHElastic - hadronic elastic scattering
  static constexpr int kHadInelasticScatt = 23; // kPHInhelastic - hadronic inelastic scattering

  // Analysis configs
  Configurable<bool> confIsRun3{"confIsRun3", false, "Running on Run3 or pilot"};
  Configurable<std::vector<int>> confPDGCodes{"confPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
  Configurable<bool> confAnalysisWithPID{"confAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles"};
  Configurable<bool> confStoreMotherPDG{"confStoreMotherPDG", false, "Store mother's PDG in tempFitVar."};

  /// Event cuts
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<float> confCentFT0Min{"confCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> confCentFT0Max{"confCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};
  Configurable<bool> confDoSpher{"confDoSpher", false, "Calculate sphericity. If false sphericity will take value of 2."};
  Configurable<bool> confIsCheckRCTFlags{"confIsCheckRCTFlags", true, "Use RCTFlags"};

  // Track cuts
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPtLowFilterCut{"confPtLowFilterCut", 0.14, "Lower limit for Pt for the filtering tracks"};   // pT low
    Configurable<float> confPtHighFilterCut{"confPtHighFilterCut", 5.0, "Higher limit for Pt for the filtering tracks"}; // pT high
    Configurable<float> confEtaFilterCut{"confEtaFilterCut", 0.8, "Eta cut for the filtering tracks"};
  } ConfFilteringTracks;

  // D0/D0bar cuts
  Configurable<float> yD0CandGenMax{"yD0CandGenMax", 0.5, "Rapidity cut for the D0/D0bar mesons"};

  FemtoUniverseCollisionSelection colCuts;
  FemtoUniverseCollisionSelection recoCollCuts;
  FemtoUniverseTrackSelection trackCuts;
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  HistogramRegistry qaRecoRegistry{"QARecoHistos", {}, OutputObjHandlingPolicy::QAObject};

  void init(InitContext&)
  {
    if (!doprocessTrackMC && !doprocessTrackMcOnlyRecoColl) {
      LOGF(fatal, "Neither processTrackMC nor processTrackMcOnlyRecoColl enabled. Please choose one.");
    }
    rctChecker.init("CBT_hadronPID", false, true);
    // MC Truth collisions
    colCuts.setCuts(confEvtZvtx, confEvtTriggerCheck, confEvtTriggerSel, confEvtOfflineCheck, confIsRun3, confCentFT0Min, confCentFT0Max);
    colCuts.init(&qaRegistry);
    // MC Reco collisions
    recoCollCuts.setCuts(confEvtZvtx, confEvtTriggerCheck, confEvtTriggerSel, confEvtOfflineCheck, confIsRun3, confCentFT0Min, confCentFT0Max);
    recoCollCuts.init(&qaRecoRegistry);

    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::CutContainerType>(&qaRegistry);

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

  template <typename CollisionType>
  bool checkMcRecoCollisions(CollisionType const& col)
  {
    if (confIsCheckRCTFlags && !rctChecker(col)) {
      return false;
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled MC Truth collision is rejected
    if (!recoCollCuts.isSelected(col)) {
      return false;
    }
    // additional checks on the reconstructed collisions
    if (col.selection_bit(aod::evsel::kNoSameBunchPileup) && col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      recoCollCuts.fillQA(col);
      return true;
    }
    return false;
  }

  template <typename ParticleType>
  int getMotherPDG(ParticleType const& particle)
  {
    if (particle.isPhysicalPrimary()) {
      return 0;
    }
    if (particle.has_mothers()) {
      if (particle.getProcess() == kHadElasticScatt || particle.getProcess() == kHadInelasticScatt) { // treat particles from hadronic scattering (20, 23) as primary
        return 0;
      }
      auto motherparticlesMC = particle.template mothers_as<aod::McParticles>();
      const auto& motherparticleMC = motherparticlesMC.front();
      return motherparticleMC.pdgCode();
    }
    return 999;
  }

  template <typename TrackType>
  void fillParticles(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children

    for (const auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (particle.pt() < ConfFilteringTracks.confPtLowFilterCut || particle.pt() > ConfFilteringTracks.confPtHighFilterCut) {
        continue;
      }

      int pdgCode = particle.pdgCode();

      if (confAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = confPDGCodes; // necessary due to some features of the Configurable
        for (auto const& pdg : tmpPDGCodes) {
          if (pdgCode == Pdg::kPhi ||          // phi meson
              std::abs(pdgCode) == Pdg::kD0 || // D0(bar) meson
              pdgCode == Pdg::kDPlus ||        // D+ meson
              (pdg == pdgCode && particle.isPhysicalPrimary())) {
            pass = true;
            break; // Exit early once a match is found
          }
        }
        if (!pass) {
          continue;
        }
      }

      // check if D0/D0bar mesons pass the rapidity cut
      // if pass then saving the orgin of D0/D0bar
      // check if tracks (besides D0/D0bar) pass pseudorapidity cut
      int8_t origin = -99;
      if (std::abs(particle.pdgCode()) == Pdg::kD0) {
        if (std::abs(particle.y()) > yD0CandGenMax) {
          continue;
        }
        origin = RecoDecay::getCharmHadronOrigin(tracks, particle);
      } else if (std::abs(particle.eta()) > ConfFilteringTracks.confEtaFilterCut) {
        continue;
      }

      /// check if we end-up with the correct final state using MC info
      int8_t sign = 0;
      if (std::abs(pdgCode) == Pdg::kD0 && !RecoDecay::isMatchedMCGen(tracks, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
        /// check if we have D0(bar) → π± K∓
        continue;
      }
      // Getting the PDG code of the mother particle
      int motherPDG = confStoreMotherPDG ? getMotherPDG(particle) : particle.pdgCode();
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
                  motherPDG);
    }
  }

  void processTrackMC(aod::McCollision const&,
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

  Preslice<aod::McParticles> perMCCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> recoCollsPerMCColl = aod::mcparticle::mcCollisionId;

  void processTrackMcOnlyRecoColl(aod::McCollisions const& mccols,
                                  soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
                                  aod::McParticles const& mcParticles,
                                  aod::BCsWithTimestamps const&)
  {
    // check on the reco collisions
    static std::set<int> mcColIds;
    mcColIds.clear();

    for (const auto& col : collisions) {
      const auto colcheck = checkMcRecoCollisions(col);
      if (colcheck) {
        mcColIds.insert(col.mcCollisionId());
      }
    }
    // filling truth collisions and particles
    for (const auto& mccol : mccols) {
      if (!mcColIds.contains(mccol.globalIndex())) {
        continue;
      }
      auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCColl, mccol.globalIndex());
      // magnetic field for run not needed for mc truth
      // fill the tables
      fillCollisions(groupedCollisions, groupedMCParticles);
      fillParticles(groupedMCParticles);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerMCTruthTask, processTrackMcOnlyRecoColl, "Provide MC data for track analysis from truth coll which were recontructed ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUniverseProducerMCTruthTask>(cfgc)};
  return workflow;
}
