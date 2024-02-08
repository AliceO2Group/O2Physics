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
#include "TFile.h"
#include "TGrid.h"
#include "TLorentzVector.h"

#include "Common/DataModel/EventSelection.h"
#include "DataFormatsPHOS/Cell.h"
#include "DataFormatsPHOS/Cluster.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "DataFormatsPHOS/BadChannelsMap.h"
#include "DataFormatsPHOS/CalibParams.h"
#include "PHOSBase/Geometry.h"
#include "PHOSReconstruction/Clusterer.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \struct PHOS Calibration
/// \brief Collectes hitograms to evaluate L1pahse, time and energy calibration, non-linearity and check
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///
struct phosCalibration {

  class photon : public TLorentzVector
  {
   public:
    photon() = default;
    explicit photon(const photon& p) = default;
    explicit photon(const TLorentzVector& p) : TLorentzVector(p), fabsId(0), fBad(false) {}
    explicit photon(double px, double py, double pz, double e) : TLorentzVector(px, py, pz, e), fabsId(0), fBad(false) {}
    void setAbsId(int16_t id) { fabsId = id; }
    void setBad() { fBad = true; }
    void setBC(uint64_t bc) { fBC = bc; }
    int16_t getAbsId() { return fabsId; }
    bool isBad() { return fBad; }
    uint64_t getBC() { return fBC; }

   protected:
    uint64_t fBC;
    int16_t fabsId;
    bool fBad;
  };

  ConfigurableAxis timeAxis{"t", {200, -50.e-9, 150.e-9}, "t (s)"};
  ConfigurableAxis timeAxisLarge{"celltime", {1000, -1500.e-9, 3500.e-9}, "cell time (ns)"};
  ConfigurableAxis timeAxisRunStart{"timeRunStart", {50000, 0., 50000.}, "Time from start of the run (s)"};
  ConfigurableAxis amplitudeAxisLarge{"amplitude", {1000, 0., 100.}, "Amplutude (GeV)"};
  ConfigurableAxis mggAxis{"mgg", {250, 0., 1.}, "m_{#gamma#gamma} (GeV/c^{2})"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<double> mMinCellAmplitude{"minCellAmplitude", 0.3, "Minimum cell amplitude for histograms."};
  Configurable<double> mCellTimeMinE{"minCellTimeAmp", 100, "Minimum cell amplitude for time histograms (ADC)"};
  Configurable<double> mMinCellTimeMain{"minCellTimeMain", -50.e-9, "Min. cell time of main bunch selection"};
  Configurable<double> mMaxCellTimeMain{"maxCellTimeMain", 100.e-9, "Max. cell time of main bunch selection"};
  Configurable<int> mMixedEvents{"mixedEvents", 10, "number of events to mix"};
  Configurable<bool> mSkipL1phase{"skipL1phase", false, "do not correct L1 phase from CCDB"};
  Configurable<std::string> rctPath{"rct-path", "RCT/Info/RunInformation", "path to the ccdb RCT objects for the SOR timestamps"};
  Configurable<std::string> mBadMapPath{"badmapPath", "PHS/Calib/BadMap", "path to BadMap snapshot"};
  Configurable<std::string> mCalibPath{"calibPath", "PHS/Calib/CalibParams", "path to Calibration snapshot"};
  Configurable<std::string> mL1PhasePath{"L1phasePath", "PHS/Calib/L1phase", "path to L1phase snapshot"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdb_api;

  HistogramRegistry mHistManager{"phosCallQAHistograms"};

  o2::phos::Geometry* geom; // singleton
  std::unique_ptr<o2::phos::Clusterer> clusterizer;
  std::vector<o2::phos::Cell> phosCells;
  std::vector<o2::phos::TriggerRecord> phosCellTRs;
  std::vector<o2::phos::CluElement> phosCluElements;
  std::vector<o2::phos::Cluster> phosClusters;
  std::vector<o2::phos::TriggerRecord> phosClusterTrigRecs;
  std::vector<photon> event;
  int mL1 = 0;

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    LOG(info) << "Initializing PHOS calibration task ...";

    const AxisSpec
      bcAxis{4, 0., 4., "bc%4"},
      modAxis{4, 1., 5., "module", "Module"},
      timeDdlAxis{200, -200.e-9, 400.e-9, "t (s)"},
      ddlAxis{14, 0., 14., "ddl"},
      cellXAxis{64, 0., 64, "x", ""},
      cellZAxis{56, 0., 56, "z", ""},
      absIdAxis{12544, 1793, 14336, "absId", "absId"};

    // Number of events
    mHistManager.add("eventsAll", "Number of events", HistType::kTH1F, {{2, 0., 2.}});
    mHistManager.add("eventsTrig", "Number of trigger events", HistType::kTH1F, {{2, 0., 2.}});

    // Cells
    mHistManager.add("cellOcc", "Cell occupancy per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("cellAmp", "Cell amplitude per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("cellTimeHG", "Time per cell, High Gain", HistType::kTH2F, {absIdAxis, timeAxis});
    mHistManager.add("cellTimeLG", "Time per cell, Low Gain", HistType::kTH2F, {absIdAxis, timeAxis});
    mHistManager.add("timeDDL", "time vs bc for DDL", HistType::kTH3F, {timeDdlAxis, bcAxis, ddlAxis});
    mHistManager.add("cellTimeFromRunStart", "time in cells vs time from start of the run", HistType::kTH3F, {absIdAxis, timeAxis, timeAxisRunStart});

    // Clusters
    mHistManager.add("hSoftClu", "Soft clu occupancy per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("hSoftCluGood", "Soft clu occupancy per module after bad map", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("hHardClu", "Hard clu occupancy per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("hSpClu", "Spectra", HistType::kTH2F, {amplitudeAxisLarge, modAxis});
    mHistManager.add("hTimeEClu", "Time vs E vs DDL", HistType::kTH3F, {ddlAxis, amplitudeAxisLarge, timeAxis});
    mHistManager.add("hTimeDdlCorr", "Time vs DDL", HistType::kTH3F, {ddlAxis, timeAxis, bcAxis});
    mHistManager.add("hRemgg", "Real m_{#gamma#gamma}", HistType::kTH2F, {absIdAxis, mggAxis});
    mHistManager.add("hMimgg", "Mixed m_{#gamma#gamma}", HistType::kTH2F, {absIdAxis, mggAxis});
    mHistManager.add("hResum", "Real m_{#gamma#gamma}", HistType::kTH2F, {mggAxis, amplitudeAxisLarge});
    mHistManager.add("hMisum", "Mixed m_{#gamma#gamma}", HistType::kTH2F, {mggAxis, amplitudeAxisLarge});

    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb_api.init(o2::base::NameConf::getCCDBServer());

    geom = o2::phos::Geometry::GetInstance("Run3");
    LOG(info) << "Calibration configured ...";
  }

  using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  /// \brief Process PHOS data
  void process(BCsWithBcSels const& bcs,
               o2::aod::Calos const& cells,
               o2::aod::CaloTriggers const& ctrs)
  {
    // Fill cell histograms
    // clusterize
    // Fill clusters histograms

    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }

    if (!clusterizer) {
      clusterizer = std::make_unique<o2::phos::Clusterer>();
      clusterizer->initialize();
    }

    const o2::phos::BadChannelsMap* badMap = ccdb->getForTimeStamp<o2::phos::BadChannelsMap>(mBadMapPath, timestamp);
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(mCalibPath, timestamp);

    if (badMap) {
      clusterizer->setBadMap(badMap);
    } else {
      LOG(fatal) << "Can not get PHOS Bad Map";
    }
    if (calibParams) {
      clusterizer->setCalibration(calibParams);
    } else {
      LOG(fatal) << "Can not get PHOS calibration";
    }

    if (!mSkipL1phase && mL1 == 0) { // should be read, but not read yet
      const std::vector<int>* vec = ccdb->getForTimeStamp<std::vector<int>>(mL1PhasePath, timestamp);
      if (vec) {
        clusterizer->setL1phase((*vec)[0]);
        mL1 = (*vec)[0];
        LOG(info) << "Got L1phase=" << mL1;
      } else {
        LOG(fatal) << "Can not get PHOS L1phase calibration";
      }
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    phosCluElements.clear();
    phosClusters.clear();
    phosClusterTrigRecs.clear();
    event.clear();

    InteractionRecord ir;
    for (auto& c : cells) {
      if (c.caloType() != 0) { // PHOS
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc_as<BCsWithBcSels>().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc_as<BCsWithBcSels>().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc_as<BCsWithBcSels>().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
        mHistManager.fill(HIST("eventsAll"), 1.);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));

      if (!c.bc_as<BCsWithBcSels>().alias_bit(mEvSelTrig))
        continue;
      mHistManager.fill(HIST("eventsTrig"), 1.);

      // Fill calibraiton histos
      if (c.amplitude() < mMinCellAmplitude)
        continue;
      char relid[3];
      o2::phos::Geometry::absToRelNumbering(c.cellNumber(), relid);
      mHistManager.fill(HIST("cellOcc"), relid[0], relid[1], relid[2]);
      mHistManager.fill(HIST("cellAmp"), relid[0], relid[1], relid[2], c.amplitude());

      int ddl = (relid[0] - 1) * 4 + (relid[1] - 1) / 16 - 2;
      uint64_t bc = c.bc_as<BCsWithBcSels>().globalBC();
      float tcorr = c.time();
      if (c.cellType() == o2::phos::HIGH_GAIN) {
        tcorr -= calibParams->getHGTimeCalib(c.cellNumber());
      } else {
        tcorr -= calibParams->getLGTimeCalib(c.cellNumber());
      }

      if (!mSkipL1phase) {
        int shift = (mL1 >> (ddl * 2)) & 3; // extract 2 bits corresponding to this ddl
        shift = bc % 4 - shift;
        if (shift < 0) {
          shift += 4;
        }
        tcorr -= shift * 25.e-9;
      }

      int runNumber = c.bc_as<BCsWithBcSels>().runNumber();
      std::map<std::string, std::string> metadata, headers;
      const std::string run_path = Form("%s/%i", rctPath.value.data(), runNumber);
      headers = ccdb_api.retrieveHeaders(run_path, metadata, -1);
      if (headers.count("SOR") == 0)
        LOGF(fatal, "Cannot find start-of-run timestamp for run number in path '%s'.", run_path.data());
      int64_t sorTimestamp = atol(headers["SOR"].c_str()); // timestamp of the SOR in ms

      if (c.amplitude() > mCellTimeMinE) {
        mHistManager.fill(HIST("timeDDL"), tcorr, bc % 4, ddl);
        mHistManager.fill(HIST("cellTimeFromRunStart"), c.cellNumber(), tcorr, (timestamp - sorTimestamp) / 1000.);
        if (c.cellType() == o2::phos::HIGH_GAIN) {
          mHistManager.fill(HIST("cellTimeHG"), c.cellNumber(), tcorr);
        } else {
          if (c.cellType() == o2::phos::LOW_GAIN) {
            mHistManager.fill(HIST("cellTimeLG"), c.cellNumber(), tcorr);
          }
        }
      }
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
      mHistManager.fill(HIST("eventsAll"), phosCellTRs.size());
    } else {
      return;
    }

    // Clusterize
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> dummyMC;
    clusterizer->processCells(phosCells, phosCellTRs, nullptr, phosClusters, phosCluElements, phosClusterTrigRecs, dummyMC);

    // Calibrate from clusters
    for (auto& tr : phosClusterTrigRecs) {
      int firstClusterInEvent = tr.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + tr.getNumberOfObjects();

      // primary vertex
      TVector3 vtx = {0., 0., 0.};

      for (int i = firstClusterInEvent; i < lastClusterInEvent; i++) {
        o2::phos::Cluster& clu = phosClusters[i];
        float e = Nonlinearity(clu.getCoreEnergy());
        if (e < 0.2) {
          continue;
        }

        int mod = clu.module();
        float x, z;
        clu.getLocalPosition(x, z);

        // find most energetic cluElement
        float maxE = 0;
        int16_t absId = 0;
        for (uint32_t iClEl = clu.getFirstCluEl(); iClEl < clu.getLastCluEl(); iClEl++) {
          auto& el = phosCluElements[iClEl];
          if (el.energy > maxE) {
            maxE = el.energy;
            absId = el.absId;
          }
        }

        char relid[3];
        phos::Geometry::absToRelNumbering(absId, relid);
        int ddl = (relid[0] - 1) * 4 + (relid[1] - 1) / 16 - 2;

        mHistManager.fill(HIST("hTimeEClu"), ddl, e, clu.getTime());
        mHistManager.fill(HIST("hSpClu"), e, mod);
        if (e > 0.5) {
          mHistManager.fill(HIST("hTimeDdlCorr"), static_cast<float>(ddl), clu.getTime(), tr.getBCData().toLong() % 4);
        }

        if (e > 0.5) {
          mHistManager.fill(HIST("hSoftClu"), mod, relid[1], relid[2]);
          if (badMap->isChannelGood(absId)) {
            mHistManager.fill(HIST("hSoftCluGood"), mod, relid[1], relid[2]);
          }
        }
        if (e > 1.5) {
          mHistManager.fill(HIST("hHardClu"), mod, relid[1], relid[2]);
        }

        if (clu.getTime() < mMinCellTimeMain || clu.getTime() > mMaxCellTimeMain) {
          continue;
        }

        TVector3 globaPos;
        // Correction for the depth of the shower starting point (TDR p 127)
        const float para = 0.925;
        const float parb = 6.52;
        float depth = para * TMath::Log(e) + parb;
        x -= x * depth / 460.;
        z -= (z - vtx.Z()) * depth / 460.;

        geom->local2Global(mod, x, z, globaPos);
        globaPos -= vtx;

        double sc = e / globaPos.Mag();
        event.emplace_back(sc * globaPos.X(), sc * globaPos.Y(), sc * globaPos.Z(), e);
        event.back().setAbsId(absId);
        event.back().setBC(tr.getBCData().toLong());
        // SetBadMap
        if (!badMap->isChannelGood(absId)) {
          event.back().setBad();
        }
      }
    }

    // Make Real and Mixed
    if (event.size() == 0) { // to avoid overflow at uint type event.size() - 1
      return;
    }
    for (size_t i = 0; i < event.size() - 1; i++) {
      int absId1 = event[i].getAbsId();
      int nMix = mMixedEvents; // Number of events to mix
      uint64_t bcurrent = 0;
      for (size_t j = i + 1; nMix && (j < event.size()); j++) {
        TLorentzVector s = event[i] + event[j];
        int absId2 = event[j].getAbsId();
        if (event[j].getBC() != bcurrent) {
          --nMix;
          bcurrent = event[j].getBC();
        }
        if (s.Pt() > 2.) {
          if (!event[j].isBad()) {
            if (event[i].getBC() == event[j].getBC()) {
              mHistManager.fill(HIST("hRemgg"), absId1, s.M());
            } else {
              mHistManager.fill(HIST("hMimgg"), absId1, s.M());
            }
          }
          if (!event[i].isBad()) {
            if (event[i].getBC() == event[j].getBC()) {
              mHistManager.fill(HIST("hRemgg"), absId2, s.M());
            } else {
              mHistManager.fill(HIST("hMimgg"), absId2, s.M());
            }
          }
        }
        if (!event[i].isBad() && !event[j].isBad()) {
          if (event[i].getBC() == event[j].getBC()) {
            mHistManager.fill(HIST("hResum"), s.M(), s.Pt());
          } else {
            mHistManager.fill(HIST("hMisum"), s.M(), s.Pt());
          }
        }
      }
    }
  }

  float Nonlinearity(float en)
  {
    // Correct for non-linearity
    // Parameters to be read from ccdb
    const double a = 9.34913e-01;
    const double b = 2.33e-03;
    const double c = -8.10e-05;
    const double d = 3.2e-02;
    const double f = -8.0e-03;
    const double g = 1.e-01;
    const double h = 2.e-01;
    const double k = -1.48e-04;
    const double l = 0.194;
    const double m = 0.0025;

    return en * (a + b * en + c * en * en + d / en + f / ((en - g) * (en - g) + h) + k / ((en - l) * (en - l) + m));
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{o2::framework::adaptAnalysisTask<phosCalibration>(cfgc)};
}
