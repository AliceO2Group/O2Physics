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

#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/LFzdcSPtables.h"

#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

constexpr double kVeryNegative = -1.e12;

using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

namespace
{
std::unordered_map<int, std::array<TH2*, 5>> gHistosC;
std::unordered_map<int, std::array<TH2*, 5>> gHistosA;
std::array<TH2*, 5> gCurrentHistC;
std::array<TH2*, 5> gCurrentHistA;
} // namespace

struct zdcSP {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Produces<aod::ZdcSPTable> zdcSPTable;
  int mRunNumber{-1};
  Configurable<float> cfgCalibrationMaxCentFT0C{"cfgCalibrationMaxCentFT0C", 90.f, "Maximum centrality FT0C"};
  Configurable<float> cfgCalibrationDownscaling{"cfgCalibrationDownscaling", 1.f, "Percentage of events to be processed"};

  ctpRateFetcher mRateFetcher;
  HistogramRegistry mHistos{"qaHists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.centFT0C() < cfgCalibrationMaxCentFT0C && gRandom->Uniform() < cfgCalibrationDownscaling;
  }

  void initCCDB(BCsRun3::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHistosC.find(mRunNumber) == gHistosC.end()) {
      std::array<TH2*, 5>& zncH = gHistosC[mRunNumber];
      std::array<TH2*, 5>& znaH = gHistosA[mRunNumber];
      zncH[0] = mHistos.add<TH2>(Form("%i/znc_common", mRunNumber), ";Hadronic rate (kHz);ZNC common signal", kTH2D, {{51, 0, 51}, {400, 0., 400.}}).get();
      znaH[0] = mHistos.add<TH2>(Form("%i/zna_common", mRunNumber), ";Hadronic rate (kHz);ZNA common signal", kTH2D, {{51, 0, 51}, {400, 0., 400.}}).get();
      for (int i{0}; i < 4; ++i) {
        zncH[i + 1] = mHistos.add<TH2>(Form("%i/znc_%i", mRunNumber, i + 1), Form(";Hadronic rate (kHz);ZNC %i signal", i + 1), kTH2D, {{51, 0, 51}, {400, 0., 400.}}).get();
        znaH[i + 1] = mHistos.add<TH2>(Form("%i/zna_%i", mRunNumber, i + 1), Form(";Hadronic rate (kHz);ZNA %i signal", i + 1), kTH2D, {{51, 0, 51}, {400, 0., 400.}}).get();
      }
    }
    gCurrentHistC = gHistosC[mRunNumber];
    gCurrentHistA = gHistosA[mRunNumber];
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
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
    auto zdc = bc.zdc();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergy = zdc.energySectorZNA();
    bool goodC = true, goodA = true;
    for (int i = 0; i < 4; i++) { // avoid std::numeric_limits<float>::infinity() in the table
      if (zncEnergy[i] < kVeryNegative) {
        zncEnergy[i] = kVeryNegative;
        goodC = false;
      }
      if (znaEnergy[i] < kVeryNegative) {
        znaEnergy[i] = kVeryNegative;
        goodA = false;
      }
    }
    float znaCommon = zdc.energyCommonZNA() < 0 ? kVeryNegative : zdc.energyCommonZNA();
    float zncCommon = zdc.energyCommonZNC() < 0 ? kVeryNegative : zdc.energyCommonZNC();
    zdcSPTable(bc.globalBC(), bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.centFT0C(),
               znaCommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3],
               zncCommon, zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
    if (!goodC && !goodA) {
      return;
    }
    initCCDB(bc);
    double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic"); //
    if (goodC) {
      gCurrentHistC[0]->Fill(hadronicRate, zncCommon);
      for (int i{0}; i < 4; ++i) {
        gCurrentHistC[i + 1]->Fill(hadronicRate, zncEnergy[i]);
      }
    }
    if (goodA) {
      gCurrentHistA[0]->Fill(hadronicRate, znaCommon);
      for (int i{0}; i < 4; ++i) {
        gCurrentHistA[i + 1]->Fill(hadronicRate, znaEnergy[i]);
      }
    }
  }
  PROCESS_SWITCH(zdcSP, processCalibrationData, "Data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcSP>(cfgc, TaskName{"zdc-sp"})};
}
