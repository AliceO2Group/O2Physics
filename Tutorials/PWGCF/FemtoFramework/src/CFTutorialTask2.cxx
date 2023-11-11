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

/// FemtoDream includes
#include "PWGCF/FemtoDream/FemtoDreamMath.h"
#include "PWGCF/FemtoDream/FemtoUtils.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFTutorialTask2 {

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

  // TODO:
  // add configurables for track selection of particle 1
  // Configurable<uint32_t> cut bit particle 1

  // add configurables for track selection of particle 2
  // Configurable<uint32_t> cut bit particle 2

  // TODO:
  // create partitions for particles 1 and 2 applying the cut bit

  HistogramRegistry HistRegistry{"FemtoTutorial", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// add cache for particle mixing
  // SliceCache ...
  // Preslice<...> perCol = ...

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

  // process same event
  void process(FilteredFDCollision const& col, FilteredFDParts const& parts)
  {

    /// event QA
    HistRegistry.fill(HIST("Event/hZvtx"), col.posZ());

    // TODO
    // generate partition of particles 1&2 with sliceByCached method

    /// TODO:
    /// loop over particle group 1
    // for (part in particles 1) {
    // fill some histograms
    // ...
    // }

    /// TODO:
    /// loop over particle group 2
    // for (part in particles 2) {
    // fill some histograms
    // ...
    // }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask2>(cfgc)};
  return workflow;
}
