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
/// \file   pidITS.cxx
/// \since  2024-11-12
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesco Mazzaschi francesco.mazzaschi@cern.ch
/// \brief  Task to produce PID tables for ITS split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include <utility>
#include <vector>
#include <string>

// O2 includes
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/BasicCCDBManager.h"
#include "TOFBase/EventTimeMaker.h"

// O2Physics includes
#include "Common/DataModel/PIDResponseITS.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TOF signal from the trackTime information
struct itsPid {

  void init(o2::framework::InitContext& initContext)
  {
  }

  /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
  void process(aod::BCs const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<itsPid>(cfgc)};
  return workflow;
}
