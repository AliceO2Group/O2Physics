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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct v0mlscoresconverter {
  Produces<aod::V0GammaMLScores> gammaMLSelections;           // gamma scores
  Produces<aod::V0LambdaMLScores> lambdaMLSelections;         // lambda scores
  Produces<aod::V0AntiLambdaMLScores> antiLambdaMLSelections; // AntiLambda scores
  Produces<aod::V0K0ShortMLScores> k0ShortMLSelections;       // K0Short scores

  void process(aod::V0Cores const& v0cores)
  {
    for (int64_t i = 0; i < v0cores.size(); ++i) {
      gammaMLSelections(-1);
      lambdaMLSelections(-1);
      antiLambdaMLSelections(-1);
      k0ShortMLSelections(-1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0mlscoresconverter>(cfgc)};
}
