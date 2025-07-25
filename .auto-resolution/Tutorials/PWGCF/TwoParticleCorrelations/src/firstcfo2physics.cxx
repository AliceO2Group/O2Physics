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

#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;

struct firsto2physics {

  void processRecoData()
  {
  }

  void processDetectorLevel()
  {
  }

  void processGen()
  {
  }

  void process(aod::Collision const&, aod::Tracks const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<firsto2physics>(cfgc),
  };
}
