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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "tables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ProduceExression {
  Defines<aod::TracksE> te;
  // bin lower bounds + the rightmost upper bound
  Configurable<std::vector<float>> binning{"binning", {0, o2::constants::math::PIHalf, o2::constants::math::PI, o2::constants::math::PI + o2::constants::math::PIHalf, o2::constants::math::TwoPI}, "Phi binning"};
  // parameters for the function as a flat array, grouped by parameter
  Configurable<std::vector<float>> parameters{"parameters", {1.0, 1.1, 1.2, 1.3,  // par 0
                                                             2.0, 2.1, 2.2, 2.3,  // par 1
                                                             3.0, 3.1, 3.2, 3.3}, // par 2
                                              "Function parameters for each phi bin"};

  void init(InitContext&)
  {
    // dummy function for benchmark (equivalent to dynamic-column-func.cxx)
    te.projectors[0] = binned((std::vector<float>)binning,
                              (std::vector<float>)parameters,
                              aod::track::phi, nsqrt(par(0) * aod::track::x * aod::track::x + par(1) * aod::track::y * aod::track::y + par(2) * aod::track::z * aod::track::z),
                              LiteralNode{-1.f});
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<ProduceExression>(cfgc)};
}
