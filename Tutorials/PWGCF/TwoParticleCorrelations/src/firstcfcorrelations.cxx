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
/// \brief Joined tables can be used as argument to the process function.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct firstcorrelations {
  void init(InitContext&)
  {
    LOGF(info, "Starting init");

    LOGF(info, "Finishing init");
  }

  template <typename TCollision, typename TTracks>
  void fillQA(TCollision collision, float centrality, TTracks tracks)
  {
  }

  template <typename TTarget, typename TCollision>
  bool fillCollision(TTarget /*target*/, TCollision /*collision*/, float /*centrality*/)
  {
    return true;
  }

  template <typename TTarget, typename TTracks>
  void fillCorrelations(TTarget target, TTracks tracks1, TTracks tracks2, float centrality, float posZ)
  {
  }

  void processSame()
  {
  }

  void processMixed()
  {
  }
  void process(aod::Collision const&, aod::Tracks const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<firstcorrelations>(cfgc),
  };
}
