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
///
/// \brief Consume task processed intermediate tables created by Preprocess task
/// \author
/// \since

#include "IntermediateTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct Consume {
  Configurable<float> ptCutOff{"ptcutoff", 10.f, "Pt cut off"};
  Configurable<float> etaRange{"etarange", 1.2f, "eta range"};

  void process(aod::Collision const&, soa::Join<aod::Tracks, aod::Decisions> const& tracks)
  {
    for (auto& track : tracks) {
      if (track.globalIndex() % 500 == 0 && track.inpt() && track.ineta()) {
        LOG(info) << "Pt = " << track.pt() << " [0, " << static_cast<float>(ptCutOff) << "); eta = " << track.eta() << "(-" << static_cast<float>(etaRange) << ", " << static_cast<float>(etaRange) << ");";
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return {
    adaptAnalysisTask<Consume>(context),
  };
}
