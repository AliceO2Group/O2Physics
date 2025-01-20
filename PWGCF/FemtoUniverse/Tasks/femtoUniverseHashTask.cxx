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

/// \file femtoUniverseHashTask.cxx
/// \brief Tasks that reads the track tables used for the pairing
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Common/Core/EventMixing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;

struct FemtoUniverseHashTask {

  Configurable<std::vector<float>> cfgVtxBins{"cfgVtxBins", std::vector<float>{-10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<std::vector<float>> cfgMultBins{"cfgMultBins", std::vector<float>{0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  // Configurable<std::vector<float>> cfgMultBins{"cfgMultBins", std::vector<float>{0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};

  std::vector<float> castCfgVtxBins, castCfgMultBins;

  Produces<aod::MixingHashes> hashes;

  void init(InitContext&)
  {
    /// here the Configurables are passed to std::vectors
    castCfgVtxBins = (std::vector<float>)cfgVtxBins;
    castCfgMultBins = (std::vector<float>)cfgMultBins;
  }

  void process(o2::aod::FdCollision const& col)
  {
    /// the hash of the collision is computed and written to table
    hashes(eventmixing::getMixingBin(castCfgVtxBins, castCfgMultBins, col.posZ(), col.multV0M()));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniverseHashTask>(cfgc)};

  return workflow;
}
