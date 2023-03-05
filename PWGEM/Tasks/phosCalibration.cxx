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

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2;
using namespace o2::framework;

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

  Configurable<double> mMinCellAmplitude{"minCellAmplitude", 0.3, "Minimum cell amplitude for histograms."};
  Configurable<double> mCellTimeMinE{"minCellTimeAmp", 100, "Minimum cell amplitude for time histograms (ADC)"};
  Configurable<double> mMinCellTimeMain{"minCellTimeMain", -50.e-9, "Min. cell time of main bunch selection"};
  Configurable<double> mMaxCellTimeMain{"maxCellTimeMain", 100.e-9, "Max. cell time of main bunch selection"};
  Configurable<int> mMixedEvents{"mixedEvents", 10, "number of events to mix"};
  Configurable<uint32_t> mL1{"L1", 0, "L1 phase"};
  Configurable<std::string> mBadMapPath{"badmapPath", "alien:///alice/cern.ch/user/p/prsnko/Calib/BadMap/snapshot.root", "path to BadMap snapshot"};
  Configurable<std::string> mCalibPath{"calibPath", "alien:///alice/cern.ch/user/p/prsnko/Calib/CalibParams/snapshot.root", "path to Calibration snapshot"};

  HistogramRegistry mHistManager{"phosCallQAHistograms"};

  o2::phos::Geometry* geom; // singleton
  std::unique_ptr<o2::phos::Clusterer> clusterizer;
  std::vector<o2::phos::Cell> phosCells;
  std::vector<o2::phos::TriggerRecord> phosCellTRs;
  std::vector<o2::phos::CluElement> phosCluElements;
  std::vector<o2::phos::Cluster> phosClusters;
  std::vector<o2::phos::TriggerRecord> phosClusterTrigRecs;
  std::vector<photon> event;

  // calibration will be set on first processing
  std::unique_ptr<const o2::phos::BadChannelsMap> badMap;   // = ccdb->get<o2::phos::BadChannelsMap>("PHS/Calib/BadMap");
  std::unique_ptr<const o2::phos::CalibParams> calibParams; // = ccdb->get<o2::phos::CalibParams>("PHS/Calib/CalibParams");

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    LOG(info) << "Initializing PHOS calibration task ...";

    const AxisSpec
      bcAxis{4, 0., 4., "bc%4"},
      modAxis{4, 1., 5., "module", "Module"},
      timeAxis{200, -50.e-9, 150.e-9, "t (s)"},
      timeDdlAxis{200, -200.e-9, 400.e-9, "t (s)"},
      timeAxisLarge{1000, -1500.e-9, 3500.e-9, "celltime", "cell time (ns)"},
      ddlAxis{14, 0., 14., "ddl"},
      cellXAxis{64, 0., 64, "x", ""},
      cellZAxis{56, 0., 56, "z", ""},
      absIdAxis{12544, 1793, 14336, "absId", "absId"},
      amplitudeAxisLarge{200, 0., 40., "amplitude", "Amplutude (GeV)"},
      mggAxis{200, 0., 0.5, "m_{#gamma#gamma} (GeV/c^{2})"};

    // Number of events
    mHistManager.add("eventsAll", "Number of events", HistType::kTH1F, {{2, 0., 2.}});

    // Cells
    mHistManager.add("cellOcc", "Cell occupancy per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("cellAmp", "Cell amplitude per module", HistType::kTH3F, {modAxis, cellXAxis, cellZAxis});
    mHistManager.add("cellTimeHG", "Time per cell, High Gain", HistType::kTH2F, {absIdAxis, timeAxis});
    mHistManager.add("cellTimeLG", "Time per cell, Low Gain", HistType::kTH2F, {absIdAxis, timeAxis});
    mHistManager.add("timeDDL", "time vs bc for DDL", HistType::kTH3F, {timeDdlAxis, bcAxis, ddlAxis});

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

    geom = o2::phos::Geometry::GetInstance("Run3");
    LOG(info) << "Calibration configured ...";
  }

  /// \brief Process PHOS data
  void process(o2::aod::BCs const& bcs,
               o2::aod::Collisions const& collisions,
               o2::aod::Calos const& cells,
               o2::aod::CaloTriggers const& ctrs)
  {
    // Fill cell histograms
    // clusterize
    // Fill clusters histograms

    if (!clusterizer) {
      clusterizer = std::make_unique<o2::phos::Clusterer>();
      clusterizer->initialize();
    }

    if (!badMap) {
      LOG(info) << "Reading BadMap from: " << mBadMapPath.value;
      TFile* fBadMap = TFile::Open(mBadMapPath.value.data());
      if (fBadMap == nullptr) { // probably, TGrid not connected yet?
        TGrid::Connect("alien");
        fBadMap = TFile::Open(mBadMapPath.value.data());
      }
      o2::phos::BadChannelsMap* bm1 = (o2::phos::BadChannelsMap*)fBadMap->Get("ccdb_object");
      badMap.reset(bm1);
      fBadMap->Close();
      clusterizer->setBadMap(badMap.get());
      LOG(info) << "Read bad map";
    }
    if (!calibParams) {
      LOG(info) << "Reading Calibration from: " << mCalibPath.value;
      TFile* fCalib = TFile::Open(mCalibPath.value.data());
      o2::phos::CalibParams* calib1 = (o2::phos::CalibParams*)fCalib->Get("ccdb_object");
      calibParams.reset(calib1);
      fCalib->Close();
      clusterizer->setCalibration(calibParams.get());
      clusterizer->setL1phase(mL1);
      LOG(info) << "Read calibration";
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    phosCluElements.clear();
    phosClusters.clear();
    phosClusterTrigRecs.clear();
    event.clear();

    // Make map between collision and BC tables
    //  map: (bcId_long,collision index)
    std::map<int64_t, int> bcMap;
    int collId = 0;
    for (auto cl : collisions) {
      bcMap[cl.bc().globalBC()] = collId;
      collId++;
    }

    InteractionRecord ir;
    for (auto& c : cells) {
      if (c.caloType() != 0) { // PHOS
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
        mHistManager.fill(HIST("eventsAll"), 1.);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));

      // Fill calibraiton histos
      if (c.amplitude() < mMinCellAmplitude)
        continue;
      char relid[3];
      o2::phos::Geometry::absToRelNumbering(c.cellNumber(), relid);
      mHistManager.fill(HIST("cellOcc"), relid[0], relid[1], relid[2]);
      mHistManager.fill(HIST("cellAmp"), relid[0], relid[1], relid[2], c.amplitude());

      int ddl = (relid[0] - 1) * 4 + (relid[1] - 1) / 16 - 2;
      uint64_t bc = c.bc().globalBC();
      int shift = (mL1 >> (ddl * 2)) & 3; // extract 2 bits corresponding to this ddl
      shift = bc % 4 - shift;
      if (shift < 0) {
        shift += 4;
      }
      float tcorr = c.time() - shift * 25.e-9;
      if (c.amplitude() > mCellTimeMinE) {
        mHistManager.fill(HIST("timeDDL"), tcorr, bc % 4, ddl);
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

      // find collision corresponding to current clu
      auto clvtx = collisions.begin() + (bcMap[tr.getBCData().toLong()]);

      // Extract primary vertex
      TVector3 vtx = {clvtx.posX(), clvtx.posY(), clvtx.posZ()};

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
