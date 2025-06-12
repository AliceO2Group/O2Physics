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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/EventSelection.h"
#include "DataFormatsPHOS/Cell.h"
#include "DataFormatsPHOS/CalibParams.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"

/// \struct PHOS QA
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since March, 2022
///
/// This task monitors simple cell quantities like
/// - Amplitude distribution
/// - Time distribution
/// - Count rate in 2D representation

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosCellQA {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  ConfigurableAxis amplitudeAxisLarge{"amplitude", {1000, 0., 10.}, "Amplutude (GeV)"};
  ConfigurableAxis timeAxisLarge{"celltime", {1000, -1500.e-9, 3500.e-9}, "cell time (ns)"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<double> mMinCellAmplitude{"minCellAmplitude", 0.5, "Minimum cell energy for histograms (GeV)"};
  Configurable<double> mMinCellTimeMain{"minCellTimeMain", -50, "Min. cell time of main bunch selection"};
  Configurable<double> mMaxCellTimeMain{"maxCellTimeMain", 100, "Max. cell time of main bunch selection"};
  Configurable<int> mVetoBCID{"vetoBCID", -1, "BC ID to be excluded"};
  Configurable<std::string> mCalibPath{"calibPath", "PHS/Calib/CalibParams", "path to Calibration snapshot"};

  o2::framework::HistogramRegistry mHistManager{"phosCallQAHistograms"};

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    LOG(info) << "Initializing PHOS QA monitor task ...";

    const o2Axis
      cellXAxis{64, 0., 64, "x", ""},
      cellZAxis{56, 0., 56, "z", ""},
      bcAxis{3501, -0.5, 3500.5};

    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cellBCAll", "Bunch crossing ID of cell (all cells)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cellBCSelected", "Bunch crossing ID of cell (selected cells)", o2HistType::kTH1F, {{bcAxis}});

    for (int i = 1; i < 5; ++i) {
      mHistManager.add(Form("cellOccM%d", i), Form("Cell occupancy for module %d", i),
                       o2HistType::kTH2F, {cellXAxis, cellZAxis});
      mHistManager.add(Form("cellAmpM%d", i), Form("Cell amplitude for module %d", i),
                       o2HistType::kTH2F, {cellXAxis, cellZAxis});
      mHistManager.add(Form("cellTimeM%d", i), Form("Cell time for module %d", i),
                       o2HistType::kTH2F, {cellXAxis, cellZAxis});
      mHistManager.add(Form("cellAmpTimeM%d", i), Form("Correlation between cell amplitude and time in module %d", i),
                       o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge});
    }
    LOG(info) << "Cell monitor task configured ...";
  }

  using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  /// \brief Process PHOS data
  void process(o2::aod::Calos const& cells, BCsWithBcSels const& bcs)
  {
    LOG(debug) << "Processing next event";

    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(mCalibPath, timestamp);

    for (const auto& bc : bcs) {
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(bc.globalBC());
      mHistManager.fill(HIST("eventsAll"), 1);
      mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
      if (mVetoBCID >= 0 && eventIR.bc == mVetoBCID)
        continue;
      mHistManager.fill(HIST("eventsSelected"), 1);
      mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
    }
    for (const auto& cell : cells) {
      if (cell.caloType() != 0)
        continue;
      o2::InteractionRecord cellIR;
      cellIR.setFromLong(cell.bc_as<BCsWithBcSels>().globalBC());
      mHistManager.fill(HIST("cellBCAll"), cellIR.bc);
      if (mVetoBCID >= 0 && cellIR.bc == mVetoBCID)
        continue;
      mHistManager.fill(HIST("cellBCSelected"), cellIR.bc);

      if (!cell.bc_as<BCsWithBcSels>().alias_bit(mEvSelTrig))
        continue;

      // bool isHighGain = cell.cellType();
      double energy = calibParams->getGain(cell.cellNumber()) * cell.amplitude();
      // if (isHighGain) {
      //    energy = calibParams->getGain(cell.cellNumber()) * cell.amplitude();
      // } else {
      //    energy = calibParams->getGain(cell.cellNumber()) * cell.amplitude() * calibParams->getHGLGRatio(cell.cellNumber());
      // }

      if (energy < mMinCellAmplitude)
        continue;
      char relid[3];
      o2::phos::Geometry::absToRelNumbering(cell.cellNumber(), relid);
      if (relid[0] == 1) {
        mHistManager.fill(HIST("cellOccM1"), relid[1] - 0.5, relid[2] - 0.5);
        mHistManager.fill(HIST("cellAmpM1"), relid[1] - 0.5, relid[2] - 0.5, energy);
        mHistManager.fill(HIST("cellTimeM1"), relid[1] - 0.5, relid[2] - 0.5, cell.time());
        if (cell.time() > mMinCellTimeMain && cell.time() < mMaxCellTimeMain) {
          mHistManager.fill(HIST("cellAmpTimeM1"), cell.time(), energy);
        }
      }
      if (relid[0] == 2) {
        mHistManager.fill(HIST("cellOccM2"), relid[1] - 0.5, relid[2] - 0.5);
        mHistManager.fill(HIST("cellAmpM2"), relid[1] - 0.5, relid[2] - 0.5, energy);
        mHistManager.fill(HIST("cellTimeM2"), relid[1] - 0.5, relid[2] - 0.5, cell.time());
        if (cell.time() > mMinCellTimeMain && cell.time() < mMaxCellTimeMain) {
          mHistManager.fill(HIST("cellAmpTimeM2"), cell.time(), energy);
        }
      }
      if (relid[0] == 3) {
        mHistManager.fill(HIST("cellOccM3"), relid[1] - 0.5, relid[2] - 0.5);
        mHistManager.fill(HIST("cellAmpM3"), relid[1] - 0.5, relid[2] - 0.5, energy);
        mHistManager.fill(HIST("cellTimeM3"), relid[1] - 0.5, relid[2] - 0.5, cell.time());
        if (cell.time() > mMinCellTimeMain && cell.time() < mMaxCellTimeMain) {
          mHistManager.fill(HIST("cellAmpTimeM3"), cell.time(), energy);
        }
      }
      if (relid[0] == 4) {
        mHistManager.fill(HIST("cellOccM4"), relid[1] - 0.5, relid[2] - 0.5);
        mHistManager.fill(HIST("cellAmpM4"), relid[1] - 0.5, relid[2] - 0.5, energy);
        mHistManager.fill(HIST("cellTimeM4"), relid[1] - 0.5, relid[2] - 0.5, cell.time());
        if (cell.time() > mMinCellTimeMain && cell.time() < mMaxCellTimeMain) {
          mHistManager.fill(HIST("cellAmpTimeM4"), cell.time(), energy);
        }
      }
    }
    LOG(debug) << "Processing event done";
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosCellQA>(cfgc)};
}
