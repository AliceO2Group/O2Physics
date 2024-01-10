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

/// \author Anton Riedel

/// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "TDatabasePDG.h"

/// FemtoDream includes
#include "PWGCF/FemtoDream/FemtoDreamMath.h"
#include "PWGCF/FemtoDream/FemtoUtils.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFTutorialTask5 {

  // Define additional analysis level cuts applied as filters
  Configurable<float> ConfZvtxMin{"ConfZvtxMin", -10, "Min Z vtx cut"};
  Configurable<float> ConfZvtxMax{"ConfZvtxMax", 10, "Max Z vtx cut"};
  Configurable<float> ConfEtaMin{"ConfEtaMin", -0.8, "Pseudorapidity cut"};
  Configurable<float> ConfEtaMax{"ConfEtaMax", 0.8, "Pseudorapidity cut"};
  Configurable<float> ConfPtMin{"ConfPtMin", 0.5, "Max Pt cut"};
  Configurable<float> ConfPtMax{"ConfPtMax", 4.0, "Min Pt cut"};

  // Defining filters
  Filter collisionFilter = (aod::collision::posZ > ConfZvtxMin) && (aod::collision::posZ < ConfZvtxMax);
  Filter trackFilter = (aod::femtodreamparticle::eta > ConfEtaMin) && (aod::femtodreamparticle::eta < ConfEtaMax) && (aod::femtodreamparticle::pt > ConfPtMin) && (aod::femtodreamparticle::pt < ConfPtMax);

  // Apply filters
  using FilteredFDCollisions = soa::Filtered<aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  using FilteredFDParts = soa::Filtered<aod::FDParticles>;
  using FilteredFDPart = FilteredFDParts::iterator;

  // selections for particles
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 3191978, "Particle 1 - Selection bit from cutCulator"};
  Configurable<int> ConfPIDPartOne{"ConfPIDPartOne", 0, "Particle 1 - Index in ConfTrkPIDspecies of producer task"};
  Configurable<int> ConfPIDValuePartOne{"ConfPIDValuePartOne", 3, "Particle 1 - Read from cutCulator"};

  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 3191978, "Particle 2 - Selection bit"};
  Configurable<int> ConfPIDPartTwo{"ConfPIDPartTwo", 0, "Particle 2 - Index in ConfTrkPIDspecies of producer task"};
  Configurable<int> ConfPIDValuePartTwo{"ConfPIDValuePartTwo", 3, "Particle 1 - Read from cutCulator"};

  Configurable<float> ConfPIDThreshold{"ConfPIDThreshold", 0.75, "Momentum threshold for TPC to TPCTOF PID"};
  Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{3.f, 3.5f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};

  /// Partitions for particle 1 and particle 2
  Partition<FilteredFDParts> PartsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Partition<FilteredFDParts> PartsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  HistogramRegistry HistRegistry{"FemtoTutorial", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// mixing
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  Configurable<int> ConfMixingDepth{"ConfMixingDepth", 10, "Number of events for mixing"};

  // TODO:
  // add binning policy for event mixing using bins in both multiplicity and position of z vertex

  // create analysis objects like histograms
  void init(o2::framework::InitContext&)
  {

    // Add histograms to histogram registry
    HistRegistry.add("Event/hZvtx", ";Z (cm)", kTH1F, {{240, -12, 12}});

    HistRegistry.add("Particle1/hPt", ";#it{p_{T}} (GeV/#it{c})", kTH1F, {{100, 0, 4}});
    HistRegistry.add("Particle1/hEta", ";#eta", kTH1F, {{100, -1., 1.}});
    HistRegistry.add("Particle1/hPhi", ";#phi", kTH1F, {{360, 0, 6.28}});

    HistRegistry.add("Particle2/hPt", ";#it{p_{T}} (GeV/#it{c})", kTH1F, {{100, 0, 4}});
    HistRegistry.add("Particle2/hEta", ";#eta", kTH1F, {{100, -1., 1.}});
    HistRegistry.add("Particle2/hPhi", ";#phi", kTH1F, {{360, 0, 6.28}});

    HistRegistry.add("Pair/hSE", ";k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});
    HistRegistry.add("Pair/hME", ";k^{*} (GeV/#it{c})", kTH1F, {{1000, 0., 5.}});
  }

  // TODO:
  // implement a process switch to run both same and mixed event processing

  // process same event
  void process(FilteredFDCollision const& col, FilteredFDParts const& parts)
  {

    /// event QA
    HistRegistry.fill(HIST("Event/hZvtx"), col.posZ());

    // generate partition of particels
    auto GroupPartsOne = PartsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto GroupPartsTwo = PartsTwo->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    /// QA for particle 1
    for (auto& part : GroupPartsOne) {
      /// check PID of particle 1 using function from FemtoUtils using PID bit
      if (isFullPIDSelected(part.pidcut(),
                            part.p(),
                            ConfPIDThreshold.value,
                            ConfPIDPartOne.value,
                            ConfNspecies.value,
                            ConfTrkPIDnSigmaMax.value,
                            ConfPIDValuePartOne.value,
                            ConfPIDValuePartOne.value)) {
        HistRegistry.fill(HIST("Particle1/hPt"), part.pt());
        HistRegistry.fill(HIST("Particle1/hEta"), part.eta());
        HistRegistry.fill(HIST("Particle1/hPhi"), part.phi());
      }
    }

    /// QA for particle 2
    /// skip QA if particle 1 & 2 are the same
    if (ConfIsSame.value == false) {
      for (auto& part : GroupPartsTwo) {
        /// check PID of particle 1 using function from FemtoUtils using PID bit
        if (isFullPIDSelected(part.pidcut(),
                              part.p(),
                              ConfPIDThreshold.value,
                              ConfPIDPartTwo.value,
                              ConfNspecies.value,
                              ConfTrkPIDnSigmaMax.value,
                              ConfPIDValuePartTwo.value,
                              ConfPIDValuePartTwo.value)) {
          HistRegistry.fill(HIST("Particle2/hPt"), part.pt());
          HistRegistry.fill(HIST("Particle2/hEta"), part.eta());
          HistRegistry.fill(HIST("Particle2/hPhi"), part.phi());
        }
      }
    }

    float kstar = 0.;
    float m0 = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartOne.value)->Mass();
    float m1 = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartTwo.value)->Mass();

    /// particle combinations
    /// if particles are the same or not determines the combination stratety
    if (ConfIsSame) {
      for (auto& [p0, p1] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(GroupPartsOne, GroupPartsTwo))) {
        if (isFullPIDSelected(p0.pidcut(),
                              p0.p(),
                              ConfPIDThreshold.value,
                              ConfPIDPartOne.value,
                              ConfNspecies.value,
                              ConfTrkPIDnSigmaMax.value,
                              ConfPIDValuePartOne.value,
                              ConfPIDValuePartOne.value) &&
            isFullPIDSelected(p1.pidcut(),
                              p1.p(),
                              ConfPIDThreshold.value,
                              ConfPIDPartOne.value,
                              ConfNspecies.value,
                              ConfTrkPIDnSigmaMax.value,
                              ConfPIDValuePartOne.value,
                              ConfPIDValuePartOne.value)

        ) {
          kstar = FemtoDreamMath::getkstar(p0, m0, p1, m1);
          HistRegistry.fill(HIST("Pair/hSE"), kstar);
        }
      }
    } else {
      for (auto& [p0, p1] : combinations(soa::CombinationsFullIndexPolicy(GroupPartsOne, GroupPartsTwo))) {
        if (isFullPIDSelected(p0.pidcut(),
                              p0.p(),
                              ConfPIDThreshold.value,
                              ConfPIDPartOne.value,
                              ConfNspecies.value,
                              ConfTrkPIDnSigmaMax.value,
                              ConfPIDValuePartOne.value,
                              ConfPIDValuePartOne.value) &&
            isFullPIDSelected(p1.pidcut(),
                              p1.p(),
                              ConfPIDThreshold.value,
                              ConfPIDPartOne.value,
                              ConfNspecies.value,
                              ConfTrkPIDnSigmaMax.value,
                              ConfPIDValuePartOne.value,
                              ConfPIDValuePartOne.value)

        ) {
          kstar = FemtoDreamMath::getkstar(p0, m0, p1, m1);
          HistRegistry.fill(HIST("Pair/hSE"), kstar);
        }
      }
    }
  }
  // PROCESS_SWITCH(...);

  // TODO:
  // Implement mixed event processing
  // abstract or reuse code from the same processing

  // void processMixedEvent(...)
  // {
  // }
  // PROCESS_SWITCH(...);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask5>(cfgc)};
  return workflow;
}
