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
/// \file   pidTOFbeta.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce TOF beta and TOF mass tables
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "pidTOFBase.h"
#include "TableHelper.h"

using namespace o2;
using namespace o2::pid;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Legacy. No effect."}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

struct tofPidBeta {
  Produces<aod::pidTOFbeta> tablePIDBeta;
  Produces<aod::pidTOFmass> tablePIDTOFMass;
  Configurable<float> expreso{"tof-expreso", 80, "Expected resolution for the computation of the expected beta"};

  bool enableTableBeta = false;
  bool enableTableMass = false;
  void init(o2::framework::InitContext& initContext)
  {
    if (isTableRequiredInWorkflow(initContext, "pidTOFbeta")) {
      enableTableBeta = true;
    }
    if (isTableRequiredInWorkflow(initContext, "pidTOFmass")) {
      enableTableMass = true;
    }
    responseBeta.mExpectedResolution = expreso.value;
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  o2::pid::tof::Beta<Trks::iterator> responseBeta;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void process(Trks const& tracks)
  {
    if (!enableTableBeta && !enableTableMass) {
      return;
    }
    float beta = 0.f;
    tablePIDBeta.reserve(tracks.size());
    for (auto const& trk : tracks) {
      beta = responseBeta.GetBeta(trk);
      if (enableTableBeta) {
        tablePIDBeta(beta,
                     responseBeta.GetExpectedSigma(trk),
                     responseBeta.GetExpectedSignal<o2::track::PID::Electron>(trk),
                     responseBeta.GetExpectedSigma(trk),
                     responseBeta.GetSeparation<o2::track::PID::Electron>(trk));
      }
      if (enableTableMass) {
        tablePIDTOFMass(o2::pid::tof::TOFMass<Trks::iterator>::GetTOFMass(trk, beta));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofPidBeta>(cfgc)};
}
