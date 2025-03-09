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

/// \file dptDptPerRunQc.cxx
/// \brief basic per run check of the ITS dead chips and of the hadronic interaction rate
/// \author victor.gonzalez.sebastian@gmail.com

#include <array>
#include <cmath>
#include <unordered_map>
#include <memory>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsITSMFT/NoiseMap.h" // missing include in TimeDeadMap.h
#include "DataFormatsITSMFT/TimeDeadMap.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"

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

namespace perrunqctask
{
static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
std::unordered_map<int, TH2*> gHadronicRate;
std::unordered_map<int, std::shared_ptr<TH1>> gCollisionOrbitBefore;
std::unordered_map<int, std::shared_ptr<TH1>> gCollisionOrbitAfter;

TH2* gCurrentHadronicRate;
std::shared_ptr<TH1> gCurrentCollisionOrbitBefore;
std::shared_ptr<TH1> gCurrentCollisionOrbitAfter;
} // namespace perrunqctask

struct DptDptPerRunQc {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber{-1};
  double mMinSeconds{-1.};

  ctpRateFetcher mRateFetcher;
  HistogramRegistry mHistos{"PerRunQaHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    using namespace perrunqctask;
    using namespace analysis::dptdptfilter;

    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHadronicRate.find(mRunNumber) == gHadronicRate.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      uint64_t mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);                /// round tsSOR to the highest integer lower than tsSOR
      double maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>(maxSec - mMinSeconds), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      gHadronicRate[mRunNumber] = mHistos.add<TH2>(Form("%i/hadronicRate", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"}, {1010, 0., 1010.}}).get();

      /* initializing the ITS chips dead map orbit axis*/
      /* inspired in DPG/Tasks/AOTEvent/eventSelectionQa.cxx */
      auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), mRunNumber);
      int64_t tsSOR = runInfo.sor;
      int64_t tsEOR = runInfo.eor;
      o2::itsmft::TimeDeadMap* itsDeadMap = ccdb->getForTimeStamp<o2::itsmft::TimeDeadMap>("ITS/Calib/TimeDeadMap", (tsSOR + tsEOR) / 2);
      auto itsDeadMapOrbits = itsDeadMap->getEvolvingMapKeys(); // roughly every second, ~350 TFs = 350x32 orbits
      if (itsDeadMapOrbits.size() > 0) {
        std::vector<double> itsDeadMapOrbitsDouble(itsDeadMapOrbits.begin(), itsDeadMapOrbits.end());
        const AxisSpec axisItsDeadMapOrbits{itsDeadMapOrbitsDouble};
        gCollisionOrbitBefore[mRunNumber] = mHistos.add<TH1>(TString::Format("%d/Before/hCollisionOrbitB", mRunNumber), "Collision orbit before; orbit", kTH1I, {axisItsDeadMapOrbits});
        gCollisionOrbitAfter[mRunNumber] = mHistos.add<TH1>(TString::Format("%d/After/hCollisionOrbitA", mRunNumber), "Collision orbit; orbit", kTH1I, {axisItsDeadMapOrbits});
        gCurrentCollisionOrbitBefore = gCollisionOrbitBefore[mRunNumber];
        gCurrentCollisionOrbitAfter = gCollisionOrbitAfter[mRunNumber];
      } else {
        gCurrentCollisionOrbitBefore = nullptr;
        gCurrentCollisionOrbitAfter = nullptr;
      }
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
    LOGF(info, "Getting out");
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::DptDptCFCollisionsInfo>::iterator const& collision, aod::BCsWithTimestamps const&)
  {
    using namespace perrunqctask;
    using namespace analysis::dptdptfilter;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    int64_t orbit = bc.globalBC() / nBCsPerOrbit;
    gCurrentCollisionOrbitBefore->Fill(orbit);

    if (!collision.collisionaccepted()) {
      return;
    }

    double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "T0VTX") * 1.e-3; //
    double seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
    gCurrentHadronicRate->Fill(seconds, hadronicRate);
    gCurrentCollisionOrbitAfter->Fill(orbit);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DptDptPerRunQc>(cfgc)};
}
