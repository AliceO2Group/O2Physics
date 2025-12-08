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

#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

void o2::analysis::PWGCF::registerConfiguration(PWGCF::FilterAndAnalysisFramework* fFilterFramework)
{
  fFilterFramework->RegisterConfiguration();
}

struct TwoParticleCorrelationsRegisterSkimming {
  bool registered = false;
  PWGCF::FilterAndAnalysisFramework* fFilterFramework = nullptr;

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h" // NOLINT

  void init(InitContext const&)
  {

    LOGF(info, "TwoParticleCorrelationsRegisterSkimming::init()");

    /* collision filtering configuration */
    PWGCF::EventSelectionConfigurable eventsel(eventfilter.bfield, eventfilter.centmultsel, {}, eventfilter.zvtxsel, eventfilter.pileuprej);
    /* track filtering configuration */
    PWGCF::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
#ifdef INCORPORATEBAYESIANPID
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
#else
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr);
#endif

    fFilterFramework = new PWGCF::FilterAndAnalysisFramework(filterccdb.ccdburl.value, filterccdb.ccdbpath.value, filterccdb.filterdate.value);
    fFilterFramework->SetConfiguration(eventsel, trksel, pidsel, PWGCF::SelectionFilterAndAnalysis::kFilter);
    registered = false;
  }

  void process(aod::Collisions const&)
  {
    if (!registered) {
      LOGF(info, "Registering filter configuration");
      PWGCF::registerConfiguration(fFilterFramework);
      registered = true;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelationsRegisterSkimming>(cfgc)};
  return workflow;
}
