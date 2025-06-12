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
#include <sstream>
#include <string>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "DataFormatsEMCAL/Constants.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "CommonDataFormat/InteractionRecord.h"

/// \struct CellMonitor
/// \brief Simple monitoring task for cell related quantities
/// \author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laoratory
/// \since Nov 12, 2021
///
/// This task monitors simple cell quantities like
/// - Amplitude distribution
/// - Time distribution
/// - Count rate in col-row space
/// - Integrated amplitude in col-row space
/// In addition, in case the task has access to reconstructed clusters,
/// it plots
/// - Frequency of a cell contributing to clusters
/// - Summed cell amplitude fraction in a cluster
///
/// The task is a direct port of PWG::EMCAL::AliEmcalCellMonitorTask from
/// AliPhysics
struct CellMonitor {

  o2::framework::Configurable<double> mMinCellAmplitude{"minCellAmplitude", 0., "Minimum cell amplitude for histograms."};
  o2::framework::Configurable<double> mMinCellTimeMain{"minCellTimeMain", -50, "Min. cell time of main bunch selection"};
  o2::framework::Configurable<double> mMaxCellTimeMain{"maxCellTimeMain", 100, "Max. cell time of main bunch selection"};
  o2::framework::Configurable<double> mMinCellAmplitudeTimeHists{"minCellAmplitudeTimeHists", 0, "Min. cell amplitude used for time distribution"};
  o2::framework::Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID to be excluded"};
  o2::framework::Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID to be included"};

  o2::framework::HistogramRegistry mHistManager{"CellMonitorHistograms"};

  // Require EMCAL cells (CALO type 1)
  o2::framework::expressions::Filter emccellfilter = o2::aod::calo::caloType == 1;

  o2::emcal::Geometry* mGeometry = nullptr;
  std::shared_ptr<o2::emcal::BadChannelMap> mBadChannels;
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  /// \brief Create output histograms and initialize geometry
  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    LOG(info) << "Initializing EMCAL cell monitor task ...";
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
    auto nCells = mGeometry->GetNCells();
    LOG(info) << "Setting up histos for " << nCells << " cells";
    int ntimeMain = static_cast<int>(mMaxCellTimeMain - mMinCellTimeMain);
    LOG(info) << "Creating histogram for main bunch time range " << mMinCellTimeMain << " ns to " << mMaxCellTimeMain << " ns with " << ntimeMain << " bins";

    const o2Axis cellAxis{nCells, -0.5, nCells - 0.5, "cellID", "Cell abs. ID"},
      amplitudeAxis{makeAmplitudeBinning(), "amplitude", "Amplitude (GeV)"},
      amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplutude (GeV)"},
      timeAxisLarge{1500, -600, 900, "celltime", "cell time (ns)"},
      timeAxisMainBunch{ntimeMain, mMinCellTimeMain, mMaxCellTimeMain, "celltimeMain", "time (ns)"},
      colAxis{48, -0.5, 47.5, "colum", "Column"},
      rowAxis{24, -0.5, 23.5, "row", "Row"},
      bcAxis{3501, -0.5, 3500.5};
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events (selected BCs)", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsTriggered", "Number of triggered events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCTriggered", "Bunch crossing ID of event (triggered events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cellBCAll", "Bunch crossing ID of cell (all cells)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cellBCSelected", "Bunch crossing ID of cell (selected cells)", o2HistType::kTH1F, {{bcAxis}});
    mHistManager.add("cellMasking", "Monitoring for masked cells", o2HistType::kTH1F, {cellAxis});
    mHistManager.add("cellFrequency", "Frequency of cell firing", o2HistType::kTH1F, {cellAxis});
    mHistManager.add("cellAmplitude", "Energy distribution per cell", o2HistType::kTH2F, {amplitudeAxis, cellAxis});
    mHistManager.add("cellAmplitudeHG", "Energy distribution per cell high gain", o2HistType::kTH2F, {amplitudeAxis, cellAxis});
    mHistManager.add("cellAmplitudeLG", "Energy distribution per cell low gain", o2HistType::kTH2F, {amplitudeAxis, cellAxis});
    mHistManager.add("cellAmplitudeCut", "Energy distribution per cell cut", o2HistType::kTH2F, {amplitudeAxis, cellAxis});
    mHistManager.add("cellTime", "Time distribution per cell", o2HistType::kTH2F, {timeAxisLarge, cellAxis});
    mHistManager.add("cellTimeMain", "Time distribution per cell for the main bunch", o2HistType::kTH2F, {timeAxisMainBunch, cellAxis});
    mHistManager.add("cellAmplitudeBC", "Cell amplitude vs. bunch crossing ID", o2HistType::kTH2F, {bcAxis, amplitudeAxis});
    mHistManager.add("celTimeBC", "Cell time vs. bunch crossing ID", o2HistType::kTH2F, {bcAxis, timeAxisLarge});
    // mCellClusterOccurrency.setObject(new TH1F("cellClusterOccurrency", "Occurrency of a cell in clusters", nCells, -0.5, nCells - 0.5));
    // mCellAmplitudeFractionCluster.setObject(new TH2F("cellAmplitudeFractionCluster", "Summed cell amplitude fraction in a cluster", nCells, -0.5, nCells - 0.5, 200, 0., 200.));

    for (int ism = 0; ism < 20; ++ism) {
      mHistManager.add(Form("cellAmplitudeSM/cellAmpSM%d", ism), Form("Integrated cell amplitudes for SM %d", ism), o2HistType::kTH2F, {colAxis, rowAxis});
      mHistManager.add(Form("cellCountSM/cellCountSM%d", ism), Form("Count rate per cell for SM %d; col; row", ism), o2HistType::kTH2F, {colAxis, rowAxis});
      mHistManager.add(Form("cellAmplitudeTime/cellAmpTimeCorrSM%d", ism), Form("Correlation between cell amplitude and time in Supermodule %d", ism), o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge});
    }
    for (int ibc = 0; ibc < 4; ++ibc) {
      mHistManager.add(Form("cellTimeBC/cellTimeBC%d", ibc), Form("Time distribution for BC mod %d", ibc), o2HistType::kTH2F, {timeAxisLarge, cellAxis});
    }
    if (mVetoBCID->length()) {
      std::stringstream parser(mVetoBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Veto BCID " << bcid;
        mVetoBCIDs.push_back(bcid);
      }
    }
    if (mSelectBCID.value != "all") {
      std::stringstream parser(mSelectBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Select BCID " << bcid;
        mSelectBCIDs.push_back(bcid);
      }
    }
    LOG(info) << "Cell monitor task configured ...";
  }

  /// \brief Process EMCAL cells
  void process(o2::aod::BC const& bc, o2::soa::Filtered<o2::aod::Calos> const& cells)
  // void process(o2::aod::BC const& bc, o2::aod::Calos const& cells)
  {
    LOG(debug) << "Processing next event";
    o2::InteractionRecord eventIR;
    eventIR.setFromLong(bc.globalBC());
    mHistManager.fill(HIST("eventsAll"), 1);
    mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
    auto bcMod4 = eventIR.bc % 4;
    if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), eventIR.bc) != mVetoBCIDs.end()) {
      return;
    }
    if (mSelectBCIDs.size() > 0 && (std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), eventIR.bc) == mSelectBCIDs.end())) {
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
    if (cells.size()) {
      // BC not necessarily triggered - determine trigger payload based presence of cell payload
      mHistManager.fill(HIST("eventsTriggered"), 1);
      mHistManager.fill(HIST("eventBCTriggered"), eventIR.bc);
    }
    for (const auto& cell : cells) {
      // cells expected to be filtered -> only EMCAL cells
      // cell.cellNumber(),
      // cell.amplitude(),
      // cell.time(),
      if (isCellMasked(cell.cellNumber())) {
        continue;
      }
      o2::InteractionRecord cellIR;
      cellIR.setFromLong(cell.bc().globalBC());
      mHistManager.fill(HIST("cellBCAll"), cellIR.bc);
      mHistManager.fill(HIST("cellBCSelected"), cellIR.bc);
      mHistManager.fill(HIST("cellAmplitude"), cell.amplitude(), cell.cellNumber());
      if (cell.cellType() == o2::emcal::channelTypeToInt(o2::emcal::ChannelType_t::HIGH_GAIN)) {
        mHistManager.fill(HIST("cellAmplitudeHG"), cell.amplitude(), cell.cellNumber());
      } else if (cell.cellType() == o2::emcal::channelTypeToInt(o2::emcal::ChannelType_t::LOW_GAIN)) {
        mHistManager.fill(HIST("cellAmplitudeLG"), cell.amplitude(), cell.cellNumber());
      }

      if (cell.amplitude() < mMinCellAmplitude)
        continue;
      mHistManager.fill(HIST("cellAmplitudeCut"), cell.amplitude(), cell.cellNumber());
      mHistManager.fill(HIST("cellFrequency"), cell.cellNumber());
      if (cell.amplitude() >= mMinCellAmplitudeTimeHists) {
        mHistManager.fill(HIST("cellTime"), cell.time(), cell.cellNumber());
        if (cell.time() > mMinCellTimeMain && cell.time() < mMaxCellTimeMain) {
          mHistManager.fill(HIST("cellTimeMain"), cell.time(), cell.cellNumber());
        }
        mHistManager.fill(HIST("celTimeBC"), cellIR.bc, cell.time());
        fillHistTimeBCMod(bcMod4, cell.cellNumber(), cell.time());
      }

      mHistManager.fill(HIST("cellAmplitudeBC"), cellIR.bc, cell.amplitude());

      // Get Cell index in eta-phi of sm
      auto [supermodule, module, phiInModule, etaInModule] = mGeometry->GetCellIndex(cell.cellNumber());
      auto [row, col] = mGeometry->GetCellPhiEtaIndexInSModule(supermodule, module, phiInModule, etaInModule);

      fillSupermoduleHistograms(supermodule, col, row, cell.amplitude(), cell.time());
    }
    LOG(debug) << "Processing event done";
  }

  void fillSupermoduleHistograms(int supermoduleID, int col, int row, double amplitude, double celltime)
  {
    // workaround to have the histogram names per supermodule
    // handled as constexpr
    switch (supermoduleID) {
      case 0:
        supermoduleHistHelper<0>(col, row, amplitude, celltime);
        break;
      case 1:
        supermoduleHistHelper<1>(col, row, amplitude, celltime);
        break;
      case 2:
        supermoduleHistHelper<2>(col, row, amplitude, celltime);
        break;
      case 3:
        supermoduleHistHelper<3>(col, row, amplitude, celltime);
        break;
      case 4:
        supermoduleHistHelper<4>(col, row, amplitude, celltime);
        break;
      case 5:
        supermoduleHistHelper<5>(col, row, amplitude, celltime);
        break;
      case 6:
        supermoduleHistHelper<6>(col, row, amplitude, celltime);
        break;
      case 7:
        supermoduleHistHelper<7>(col, row, amplitude, celltime);
        break;
      case 8:
        supermoduleHistHelper<8>(col, row, amplitude, celltime);
        break;
      case 9:
        supermoduleHistHelper<9>(col, row, amplitude, celltime);
        break;
      case 10:
        supermoduleHistHelper<10>(col, row, amplitude, celltime);
        break;
      case 11:
        supermoduleHistHelper<11>(col, row, amplitude, celltime);
        break;
      case 12:
        supermoduleHistHelper<12>(col, row, amplitude, celltime);
        break;
      case 13:
        supermoduleHistHelper<13>(col, row, amplitude, celltime);
        break;
      case 14:
        supermoduleHistHelper<14>(col, row, amplitude, celltime);
        break;
      case 15:
        supermoduleHistHelper<15>(col, row, amplitude, celltime);
        break;
      case 16:
        supermoduleHistHelper<16>(col, row, amplitude, celltime);
        break;
      case 17:
        supermoduleHistHelper<17>(col, row, amplitude, celltime);
        break;
      case 18:
        supermoduleHistHelper<18>(col, row, amplitude, celltime);
        break;
      case 19:
        supermoduleHistHelper<19>(col, row, amplitude, celltime);
        break;
      default:
        break;
    }
  }

  void fillHistTimeBCMod(int bcMod4, int cellAbsID, double cellTime)
  {
    switch (bcMod4) {
      case 0:
        bcmodHistHelper<0>(cellAbsID, cellTime);
        break;
      case 1:
        bcmodHistHelper<1>(cellAbsID, cellTime);
        break;
      case 2:
        bcmodHistHelper<2>(cellAbsID, cellTime);
        break;
      case 3:
        bcmodHistHelper<3>(cellAbsID, cellTime);
        break;
    }
  }

  template <int supermoduleID>
  void supermoduleHistHelper(int col, int row, double amplitude, double celltime)
  {
    static constexpr std::string_view cellAmpHistSM[20] = {"cellAmplitudeSM/cellAmpSM0", "cellAmplitudeSM/cellAmpSM1", "cellAmplitudeSM/cellAmpSM2", "cellAmplitudeSM/cellAmpSM3", "cellAmplitudeSM/cellAmpSM4", "cellAmplitudeSM/cellAmpSM5", "cellAmplitudeSM/cellAmpSM6", "cellAmplitudeSM/cellAmpSM7", "cellAmplitudeSM/cellAmpSM8", "cellAmplitudeSM/cellAmpSM9", "cellAmplitudeSM/cellAmpSM10", "cellAmplitudeSM/cellAmpSM11", "cellAmplitudeSM/cellAmpSM12", "cellAmplitudeSM/cellAmpSM13", "cellAmplitudeSM/cellAmpSM14", "cellAmplitudeSM/cellAmpSM15", "cellAmplitudeSM/cellAmpSM16", "cellAmplitudeSM/cellAmpSM17", "cellAmplitudeSM/cellAmpSM18", "cellAmplitudeSM/cellAmpSM19"};
    static constexpr std::string_view cellCountHistSM[20] = {"cellCountSM/cellCountSM0", "cellCountSM/cellCountSM1", "cellCountSM/cellCountSM2", "cellCountSM/cellCountSM3", "cellCountSM/cellCountSM4", "cellCountSM/cellCountSM5", "cellCountSM/cellCountSM6", "cellCountSM/cellCountSM7", "cellCountSM/cellCountSM8", "cellCountSM/cellCountSM9", "cellCountSM/cellCountSM10", "cellCountSM/cellCountSM11", "cellCountSM/cellCountSM12", "cellCountSM/cellCountSM13", "cellCountSM/cellCountSM14", "cellCountSM/cellCountSM15", "cellCountSM/cellCountSM16", "cellCountSM/cellCountSM17", "cellCountSM/cellCountSM18", "cellCountSM/cellCountSM19"};
    static constexpr std::string_view cellAmpTimeHist[20] = {"cellAmplitudeTime/cellAmpTimeCorrSM0", "cellAmplitudeTime/cellAmpTimeCorrSM1", "cellAmplitudeTime/cellAmpTimeCorrSM2", "cellAmplitudeTime/cellAmpTimeCorrSM3", "cellAmplitudeTime/cellAmpTimeCorrSM4", "cellAmplitudeTime/cellAmpTimeCorrSM5", "cellAmplitudeTime/cellAmpTimeCorrSM6", "cellAmplitudeTime/cellAmpTimeCorrSM7", "cellAmplitudeTime/cellAmpTimeCorrSM8", "cellAmplitudeTime/cellAmpTimeCorrSM9", "cellAmplitudeTime/cellAmpTimeCorrSM10", "cellAmplitudeTime/cellAmpTimeCorrSM11", "cellAmplitudeTime/cellAmpTimeCorrSM12", "cellAmplitudeTime/cellAmpTimeCorrSM13", "cellAmplitudeTime/cellAmpTimeCorrSM14", "cellAmplitudeTime/cellAmpTimeCorrSM15", "cellAmplitudeTime/cellAmpTimeCorrSM16", "cellAmplitudeTime/cellAmpTimeCorrSM17", "cellAmplitudeTime/cellAmpTimeCorrSM18", "cellAmplitudeTime/cellAmpTimeCorrSM19"};
    mHistManager.fill(HIST(cellAmpHistSM[supermoduleID]), col, row);
    mHistManager.fill(HIST(cellCountHistSM[supermoduleID]), col, row, amplitude);
    mHistManager.fill(HIST(cellAmpTimeHist[supermoduleID]), celltime, amplitude);
  }

  template <int bcMod>
  void bcmodHistHelper(int cellAbsID, double celltime)
  {
    static constexpr std::string_view timeHistBC[4] = {"cellTimeBC/cellTimeBC0", "cellTimeBC/cellTimeBC1", "cellTimeBC/cellTimeBC2", "cellTimeBC/cellTimeBC3"};
    mHistManager.fill(HIST(timeHistBC[bcMod]), celltime, cellAbsID);
  }

  /// \brief Check if a cell is masked
  /// \param towerID ID of of the tower to be checked
  /// \return true if the tower is masked, false if the tower is considered good
  ///
  /// Requires presence of the bad channel map from the CCDB
  bool isCellMasked(int towerID)
  {
    bool masked = false;
    if (mBadChannels) {
      auto maskStatus = mBadChannels->getChannelStatus(towerID);
      masked = (maskStatus != o2::emcal::BadChannelMap::MaskType_t::GOOD_CELL);
    }
    return masked;
  }

  /// \brief Create binning for cell amplitude axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double> makeAmplitudeBinning() const
  {
    auto fillBinLimits = [](std::vector<double>& binlimits, double max, double binwidth) {
      auto current = *binlimits.rbegin();
      while (current < max) {
        current += binwidth;
        binlimits.emplace_back(current);
      }
    };
    std::vector<double> result = {0.};
    fillBinLimits(result, 2., 0.1);
    fillBinLimits(result, 5., 0.2);
    fillBinLimits(result, 10., 0.5);
    fillBinLimits(result, 20., 1.);
    fillBinLimits(result, 50., 2.);
    fillBinLimits(result, 100., 5.);
    fillBinLimits(result, 200., 10.);
    return result;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<CellMonitor>(cfgc)};
}
