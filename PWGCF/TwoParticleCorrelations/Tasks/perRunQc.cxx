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
//
// Minimal example to run this task:
// o2-analysis-centrality-table -b --configuration json://configuration.json | o2-analysis-timestamp -b --configuration json://configuration.json | o2-analysis-event-selection -b --configuration json://configuration.json | o2-analysis-multiplicity-table -b --configuration json://configuration.json | o2-analysis-lf-zdcsp -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json

#include <array>
#include <cmath>

#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

namespace perrunqatask
{
std::unordered_map<int, TH2*> gHadronicRate;

TH2* gCurrentHadronicRate;
} // namespace perrunqatask

struct DptDptPerRunQa {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber{-1};
  uint64_t mSOR{0};
  double mMinSeconds{-1.};

  ctpRateFetcher mRateFetcher;
  HistogramRegistry mHistos{"PerRunQaHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    using namespace perrunqatask;
    using namespace analysis::dptdptfilter;

    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHadronicRate.find(mRunNumber) == gHadronicRate.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);                /// round tsSOR to the highest integer lower than tsSOR
      double maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>(maxSec - mMinSeconds), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      gHadronicRate[mRunNumber] = mHistos.add<TH2>(Form("%i/hadronicRate", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"}, {1010, 0., 1010.}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  void process(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>::iterator const& collision, aod::BCsWithTimestamps const&)
  {
    using namespace perrunqatask;
    using namespace analysis::dptdptfilter;

    if (!collision.collisionaccepted()) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "T0VTX") * 1.e-3; //
    double seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
    gCurrentHadronicRate->Fill(seconds, hadronicRate);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DptDptPerRunQa>(cfgc)};
}
