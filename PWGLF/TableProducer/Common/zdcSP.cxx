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

#include "PWGLF/DataModel/LFzdcSPtables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsTPC/BetheBlochAleph.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <string>
#include <unordered_map>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

constexpr double kVeryNegative = -1.e12;

using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

namespace
{
std::unordered_map<int, std::array<TH2*, 2>> gCentroidC;
std::unordered_map<int, std::array<TH2*, 2>> gCentroidA;
std::unordered_map<int, TH2*> gHadronicRate;
std::array<TH2*, 2> gCurrentCentroidC;
std::array<TH2*, 2> gCurrentCentroidA;
TH2* gCurrentHadronicRate;
} // namespace

struct zdcSP {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Produces<aod::ZdcSPTable> zdcSPTable;
  int mRunNumber{-1};
  uint64_t mSOR{0};
  double mMinSeconds{-1.};
  Configurable<float> cfgCalibrationMaxCentFT0C{"cfgCalibrationMaxCentFT0C", 90.f, "Maximum centrality FT0C"};
  Configurable<float> cfgCalibrationDownscaling{"cfgCalibrationDownscaling", 1.f, "Percentage of events to be processed"};

  ctpRateFetcher mRateFetcher;
  HistogramRegistry mHistos{"qaHists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.centFT0C() < cfgCalibrationMaxCentFT0C;
  }

  void initCCDB(BCsRun3::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gCentroidC.find(mRunNumber) == gCentroidC.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);                /// round tsSOR to the highest integer lower than tsSOR
      double maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>(maxSec - mMinSeconds), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      std::array<TH2*, 2>& zncH = gCentroidC[mRunNumber];
      std::array<TH2*, 2>& znaH = gCentroidA[mRunNumber];
      std::string xy[2] = {"X", "Y"};
      for (int i{0}; i < 2; ++i) {
        zncH[i] = mHistos.add<TH2>(Form("%i/centroidC%s", mRunNumber, xy[i].data()), Form(";Time since SOR (s);ZNC centroid %s", xy[i].data()), kTH2D, {axisSeconds, {50, -1.5, 1.5}}).get();
        znaH[i] = mHistos.add<TH2>(Form("%i/centroidA%s", mRunNumber, xy[i].data()), Form(";Time since SOR (s);ZNA centroid %s", xy[i].data()), kTH2D, {axisSeconds, {50, -1.5, 1.5}}).get();
      }
      gHadronicRate[mRunNumber] = mHistos.add<TH2>(Form("%i/hadronicRate", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"}, {510, 0., 51.}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
    gCurrentCentroidC = gCentroidC[mRunNumber];
    gCurrentCentroidA = gCentroidA[mRunNumber];
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  Preslice<aod::Zdcs> zdcPerCollision = aod::collision::bcId;
  void processCalibrationData(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>::iterator const& collision, BCsRun3 const& bcs, aod::Zdcs const&)
  {
    if (bcs.size() != 0) {
      gRandom->SetSeed(bcs.iteratorAt(0).globalBC());
    }
    if (!eventSelection(collision)) {
      return;
    }

    auto bc = collision.foundBC_as<BCsRun3>();
    if (!bc.has_zdc()) {
      return;
    }
    initCCDB(bc);
    double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; //
    double seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
    gCurrentHadronicRate->Fill(seconds, hadronicRate);
    auto zdc = bc.zdc();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergy = zdc.energySectorZNA();
    for (int i = 0; i < 4; i++) { // avoid std::numeric_limits<float>::infinity() in the table
      if (zncEnergy[i] < kVeryNegative) {
        zncEnergy[i] = -1.f;
      }
      if (znaEnergy[i] < kVeryNegative) {
        znaEnergy[i] = -1.f;
      }
    }
    float znaCommon = zdc.energyCommonZNA() < 0 ? -1.f : zdc.energyCommonZNA();
    float zncCommon = zdc.energyCommonZNC() < 0 ? -1.f : zdc.energyCommonZNC();
    if (gRandom->Uniform() < cfgCalibrationDownscaling) {
      zdcSPTable(bc.timestamp() - mSOR, bc.runNumber(), hadronicRate, collision.posX(), collision.posY(), collision.posZ(), collision.centFT0C(),
                 znaCommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3],
                 zncCommon, zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
    }

    ///          | 3 | 4 |
    /// beam out | 1 | 2 | --> x       (ZNC)
    ///
    ///          | 3 | 4 |
    ///    <-- x | 1 | 2 |   beam in   (ZNA)

    constexpr float beamEne = 5.36 * 0.5;
    constexpr float x[4] = {-1.75, 1.75, -1.75, 1.75};
    constexpr float y[4] = {-1.75, -1.75, 1.75, 1.75};
    constexpr float alpha = 0.395;
    float numXZNC = 0., numYZNC = 0., denZNC = 0.;
    float numXZNA = 0., numYZNA = 0., denZNA = 0.;
    //
    for (int i = 0; i < 4; i++) {
      if (zncEnergy[i] > 0.) {
        float wZNC = std::pow(zncEnergy[i], alpha);
        numXZNC -= x[i] * wZNC;
        numYZNC += y[i] * wZNC;
        denZNC += wZNC;
      }
      if (znaEnergy[i] > 0.) {
        float wZNA = std::pow(znaEnergy[i], alpha);
        numXZNA += x[i] * wZNA;
        numYZNA += y[i] * wZNA;
        denZNA += wZNA;
      }
    }
    //
    float centrZNC[2], centrZNA[2];
    if (denZNC != 0.) {
      float nSpecnC = zncCommon / beamEne;
      float cZNC = 1.89358 - 0.71262 / (nSpecnC + 0.71789);
      centrZNC[0] = cZNC * numXZNC / denZNC;
      centrZNC[1] = cZNC * numYZNC / denZNC;
    } else {
      centrZNC[0] = centrZNC[1] = 999.;
    }
    //
    if (denZNA != 0.) {
      float nSpecnA = znaCommon / beamEne;
      float cZNA = 1.89358 - 0.71262 / (nSpecnA + 0.71789);
      centrZNA[0] = cZNA * numXZNA / denZNA;
      centrZNA[1] = cZNA * numYZNA / denZNA;
    } else {
      centrZNA[0] = centrZNA[1] = 999.;
    }
    for (int i{0}; i < 2; ++i) {
      gCurrentCentroidC[i]->Fill(seconds, centrZNC[i]);
      gCurrentCentroidA[i]->Fill(seconds, centrZNA[i]);
    }
  }
  PROCESS_SWITCH(zdcSP, processCalibrationData, "Data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcSP>(cfgc, TaskName{"zdc-sp"})};
}
