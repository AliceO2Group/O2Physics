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
/// \brief Filters are used to select specific rows of a table.
/// \author Anton Alkin (anton.alkin@cern.ch)
/// \file filters.cxx

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Apply filters on Collisions, Tracks, and TPhi
struct Filters {
  Configurable<float> ptLow{"ptLow", 0.5f, ""};
  Configurable<float> ptUp{"ptUp", 2.0f, ""};
  Filter ptFilterA = aod::track::pt > ptLow;
  Filter ptFilterB = aod::track::pt < ptUp;

  Configurable<float> etaLow{"etaLow", -1.0f, ""};
  Configurable<float> etaUp{"etaUp", 1.0f, ""};
  Filter etafilter = (aod::track::eta < etaUp) && (aod::track::eta > etaLow);

  Configurable<float> phiLow{"phiLow", 1.0f, "Phi lower limit"};
  Configurable<float> phiUp{"phiUp", 2.0f, "Phi upper limit"};

  Configurable<float> vtxZ{"vtxZ", 10.f, ""};
  Filter posZfilter = nabs(aod::collision::posZ) < vtxZ;
  Filter bitwiseFilter = (aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor);

  // it is now possible to set filters as strings
  // note that column designators need the full prefix, i.e. o2::aod::
  // configurables can be used with ncfg(type, value, name)
  // where value is the default value
  // name is the full name in JSON, with prefix if there is any
  Configurable<std::string> extraFilter{"extraFilter", "(o2::aod::track::phi < ncfg(float,2.0,phiUp)) && (o2::aod::track::phi > ncfg(float,1.0,phiLow))", "extra filter string"};
  Filter extraF;

  void init(InitContext&)
  {
    if (!extraFilter->empty()) {
      // string-based filters need to be assigned in init()
      extraF = Parser::parse(extraFilter);
    }
  }

  // process only collisions and tracks which pass all defined filter criteria
  void process(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra>> const& tracks)
  {
    LOGF(info, "Collision: %d [N = %d out of %d], -%.1f < %.3f < %.1f",
         collision.globalIndex(), tracks.size(), tracks.tableSize(), (float)vtxZ, collision.posZ(), (float)vtxZ);
    for (auto const& track : tracks) {
      LOGP(info, "id = {}; eta:  {} < {} < {}; phi: {} < {} < {}; pt: {} < {} < {}",
           track.collisionId(), (float)etaLow, track.eta(), (float)etaUp, (float)phiLow, track.phi(), (float)phiUp, (float)ptLow, track.pt(), (float)ptUp);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<Filters>(cfgc)};
}
