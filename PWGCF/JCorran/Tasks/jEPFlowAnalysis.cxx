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

/// \author Maxim Virta (maxim.virta@cern.ch)
/// \brief flow measurement with q-vectors
/// \file jEPFlowAnalysis.cxx
/// \since Jul 2024

#include "FlowJHistManager.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <TDatabasePDG.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct jEPFlowAnalysis {
  HistogramRegistry epFlowHistograms{"EPFlow", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  EventPlaneHelper helperEP;
  FlowJHistManager histManager;
  bool debug = kFALSE;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Set Configurables here
  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track seletion."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "event selection"};
  Configurable<int> cfgEvtSel{"cfgEvtSel", 0, "Event selection flags\n0: Sel8\n1: Sel8+kIsGoodZvtxFT0vsPV+kNoSameBunchPileup\n2: Sel8+kIsGoodZvtxFT0vsPV+kNoSameBunchPileup+kNoCollInTimeRangeStandard\n3: Sel8+kNoSameBunchPileup"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<bool> cfgEffCor{"cfgEffCor", false, "flag for efficiency correction"};
  Configurable<std::string> cfgEffCorDir{"cfgEffCorDir", "Users/n/nmallick/Run3OO/Eff/LHC25h3b_FT0C", "path for efficiency correction"};

  Configurable<bool> cfgSystStudy{"cfgSystStudy", false, "flag for syst study"};
  Configurable<int> cfgITSNCls{"cfgITSNCls", 5, "minimum number of its clusters"};
  Configurable<int> cfgTPCNclsCR{"cfgTPCNclsCR", 70, "minimum number of tpc cluster crossed rows"};
  Configurable<float> cfgTPCChi2{"cfgTPCChi2", 4.0, "maximum TPC chi2"};
  Configurable<float> cfgITSChi2{"cfgITSChi2", 36.0, "maximum ITS chi2"};
  Configurable<float> cfgdcaZ{"cfgdcaZ", 2.0, "maximum dca z"};
  Configurable<float> cfgdcaXY0{"cfgdcaXY0", 0.0105, "maximum constant dca xy"};
  Configurable<float> cfgdcaXY1{"cfgdcaXY1", 0.035, "maximum pt deepdent dca xy"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "Total number of detectors in qVectorsTable"};
  Configurable<int> cfgnMode{"cfgnMode", 1, "the number of modulations"};

  Configurable<bool> cfgManShiftCorr{"cfgManShiftCorr", false, "additional shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};
  Configurable<float> cfgVertexZ{"cfgVertexZ", 10.0, "Maximum vertex Z selection"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {100, 0, 100}, ""};
  ConfigurableAxis cfgAxisPt{"cfgAxisPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 30.0, 50.0, 70.0, 100.0}, ""};
  ConfigurableAxis cfgAxisCos{"cfgAxisCos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgAxisQvec{"cfgAxisQvec", {200, -5.0, 5.0}, ""};

  ConfigurableAxis cfgAxisCentMC{"cfgAxisCentMC", {5, 0, 100}, ""};
  ConfigurableAxis cfgAxisVtxZMC{"cfgAxisVtxZMC", {20, -10, 10}, ""};
  ConfigurableAxis cfgAxisEtaMC{"cfgAxisEtaMC", {20, -1, 1}, ""};
  ConfigurableAxis cfgAxisPhiMC{"cfgAxisPhiMC", {36, 0, constants::math::PI * 2.0}, ""};
  ConfigurableAxis cfgAxisPtMC{"cfgAxisPtMC", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 30.0, 50.0, 70.0, 100.0}, ""};

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using MyCollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>;
  using MyTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;

  int detId;
  int refAId;
  int refBId;
  int harmInd;

  int currentRunNumber = -999;
  int lastRunNumber = -999;

  std::vector<TProfile3D*> shiftprofile{};
  std::string fullCCDBShiftCorrPath;

  THn* effMap = nullptr;

  template <typename T>
  int getdetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCPos") {
      return 4;
    } else if (name.value == "TPCNeg") {
      return 5;
    } else if (name.value == "TPCTot") {
      return 6;
    } else {
      return 0;
    }
  }

  template <typename Trk>
  uint8_t trackSel(const Trk& track)
  {
    uint8_t tracksel = 0;
    if (!track.isGlobalTrack()) {
      tracksel += 1;
    }
    if (track.itsNCls() <= cfgITSNCls && cfgSystStudy) {
      tracksel += 2;
    }
    if (track.tpcNClsCrossedRows() <= cfgTPCNclsCR && cfgSystStudy) {
      tracksel += 4;
    }
    if (track.tpcChi2NCl() >= cfgTPCChi2 && cfgSystStudy) {
      tracksel += 8;
    }
    if (track.itsChi2NCl() >= cfgITSChi2 && cfgSystStudy) {
      tracksel += 16;
    }
    if (std::abs(track.dcaZ()) >= cfgdcaZ && cfgSystStudy) {
      tracksel += 32;
    }
    if (std::abs(track.dcaXY()) >= cfgdcaXY0 + cfgdcaXY1 / std::pow(track.pt(), 1.1) && cfgSystStudy) {
      tracksel += 64;
    }
    if (track.pt() <= cfgTrackCuts.cfgPtMin) {
      tracksel += 128;
    }
    if (std::abs(track.eta()) >= cfgTrackCuts.cfgEtaMax) {
      tracksel += 256;
    }

    return tracksel;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);

    return eff->GetBinContent(effVars);
  }

  void init(InitContext const&)
  {
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    detId = getdetId(cfgDetName);
    refAId = getdetId(cfgRefAName);
    refBId = getdetId(cfgRefBName);

    AxisSpec axisMod{cfgnMode, 2., cfgnMode + 2.};
    AxisSpec axisEvtPl{360, -constants::math::PI * 1.1, constants::math::PI * 1.1};
    AxisSpec axisVertex{150, -12.5, 12.5};

    AxisSpec axisCent{cfgAxisCent, "cent"};
    AxisSpec axisPt{cfgAxisPt, "pT"};
    AxisSpec axisCos{cfgAxisCos, "cos"};
    AxisSpec axisQvec{cfgAxisQvec, "Qvec"};

    AxisSpec axisCentMC{cfgAxisCentMC, "cent"};
    AxisSpec axisVtxZMC{cfgAxisVtxZMC, "vtxz"};
    AxisSpec axisEtaMC{cfgAxisEtaMC, "eta"};
    AxisSpec axisPhiMC{cfgAxisPhiMC, "phi"};
    AxisSpec axisPtMC{cfgAxisPtMC, "pt"};

    epFlowHistograms.add("EpDet", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpRefA", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpRefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});

    epFlowHistograms.add("EpResDetRefA", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpResDetRefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpResRefARefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});

    epFlowHistograms.add("vncos", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisCos}});
    epFlowHistograms.add("vnsin", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisCos}});

    epFlowHistograms.add("EpResQvecDetRefAxx", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});
    epFlowHistograms.add("EpResQvecDetRefAxy", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});
    epFlowHistograms.add("EpResQvecDetRefBxx", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});
    epFlowHistograms.add("EpResQvecDetRefBxy", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});
    epFlowHistograms.add("EpResQvecRefARefBxx", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});
    epFlowHistograms.add("EpResQvecRefARefBxy", "", {HistType::kTH3F, {axisMod, axisCent, axisQvec}});

    epFlowHistograms.add("SPvnxx", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisQvec}});
    epFlowHistograms.add("SPvnxy", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisQvec}});

    epFlowHistograms.add("hCentrality", "", {HistType::kTH1F, {axisCent}});
    epFlowHistograms.add("hVertex", "", {HistType::kTH1F, {axisVertex}});

    epFlowHistograms.add("MC/hPartGen", "", {kTHnSparseF, {cfgAxisCentMC, cfgAxisVtxZMC, cfgAxisEtaMC, cfgAxisPhiMC, cfgAxisPtMC}});
    epFlowHistograms.add("MC/hPartRecPr", "", {kTHnSparseF, {cfgAxisCentMC, cfgAxisVtxZMC, cfgAxisEtaMC, cfgAxisPhiMC, cfgAxisPtMC}});
    epFlowHistograms.add("MC/hPartRec", "", {kTHnSparseF, {cfgAxisCentMC, cfgAxisVtxZMC, cfgAxisEtaMC, cfgAxisPhiMC, cfgAxisPtMC}});
  }

  void processDefault(MyCollisions::iterator const& coll, soa::Filtered<MyTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    if (cfgAddEvtSel) {
      if (std::abs(coll.posZ()) > cfgVertexZ)
        return;
      switch (cfgEvtSel) {
        case 0: // Sel8
          if (!coll.sel8())
            return;
          break;
        case 1: // PbPb standard
          if (!coll.sel8() || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 2: // PbPb with pileup
          if (!coll.sel8() || !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ||
              !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 3: // Small systems (OO, NeNe, pp)
          if (!coll.sel8() || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        default:
          LOGF(warning, "Event selection flag was not found, continuing without basic event selections!\n");
      }
      // Check occupancy
      if (coll.trackOccupancyInTimeRange() > cfgMaxOccupancy || coll.trackOccupancyInTimeRange() < cfgMinOccupancy)
        return;
    }

    if (cfgEffCor) {
      auto bc = coll.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        effMap = ccdb->getForTimeStamp<THnT<float>>(cfgEffCorDir, bc.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    float cent = coll.cent();
    epFlowHistograms.fill(HIST("hCentrality"), cent);
    epFlowHistograms.fill(HIST("hVertex"), coll.posZ());
    float eps[3] = {0.};
    float qx_shifted[3] = {0.};
    float qy_shifted[3] = {0.};

    if (cfgManShiftCorr) {
      auto bc = coll.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        shiftprofile.clear();
        for (int i = 0; i < cfgnMode; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i + 2);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }

    if (coll.qvecAmp()[detId] < 1e-5 || coll.qvecAmp()[refAId] < 1e-5 || coll.qvecAmp()[refBId] < 1e-5)
      return;

    for (int i = 0; i < cfgnMode; i++) {       // loop over different harmonic orders
      harmInd = cfgnTotalSystem * 4 * (i) + 3; // harmonic index to access corresponding Q-vector as all Q-vectors are in same vector
      eps[0] = helperEP.GetEventPlane(coll.qvecRe()[4 * detId + harmInd], coll.qvecIm()[4 * detId + harmInd], i + 2);
      eps[1] = helperEP.GetEventPlane(coll.qvecRe()[4 * refAId + harmInd], coll.qvecIm()[4 * refAId + harmInd], i + 2);
      eps[2] = helperEP.GetEventPlane(coll.qvecRe()[4 * refBId + harmInd], coll.qvecIm()[4 * refBId + harmInd], i + 2);

      auto deltapsiDet = 0.0;
      auto deltapsiRefA = 0.0;
      auto deltapsiRefB = 0.0;

      float weight = 1.0;

      qx_shifted[0] = coll.qvecRe()[4 * detId + harmInd];
      qy_shifted[0] = coll.qvecIm()[4 * detId + harmInd];
      qx_shifted[1] = coll.qvecRe()[4 * refAId + harmInd];
      qy_shifted[1] = coll.qvecIm()[4 * refAId + harmInd];
      qx_shifted[2] = coll.qvecRe()[4 * refBId + harmInd];
      qy_shifted[2] = coll.qvecIm()[4 * refBId + harmInd];

      if (cfgManShiftCorr) {
        constexpr int kShiftBins = 10;
        for (int ishift = 1; ishift <= kShiftBins; ishift++) {
          auto coeffshiftxDet = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * detId + 0.5, ishift - 0.5));
          auto coeffshiftyDet = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * detId + 1.5, ishift - 0.5));
          auto coeffshiftxRefA = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * refAId + 0.5, ishift - 0.5));
          auto coeffshiftyRefA = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * refAId + 1.5, ishift - 0.5));
          auto coeffshiftxRefB = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * refBId + 0.5, ishift - 0.5));
          auto coeffshiftyRefB = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.0 * refBId + 1.5, ishift - 0.5));

          deltapsiDet += ((2. / (1.0 * ishift)) * (-coeffshiftxDet * std::cos(ishift * static_cast<float>(i + 2) * eps[0]) + coeffshiftyDet * std::sin(ishift * static_cast<float>(i + 2) * eps[0]))) / static_cast<float>(i + 2);
          deltapsiRefA += ((2. / (1.0 * ishift)) * (-coeffshiftxRefA * std::cos(ishift * static_cast<float>(i + 2) * eps[1]) + coeffshiftyRefA * std::sin(ishift * static_cast<float>(i + 2) * eps[1]))) / static_cast<float>(i + 2);
          deltapsiRefB += ((2. / (1.0 * ishift)) * (-coeffshiftxRefB * std::cos(ishift * static_cast<float>(i + 2) * eps[2]) + coeffshiftyRefB * std::sin(ishift * static_cast<float>(i + 2) * eps[2]))) / static_cast<float>(i + 2);
        }

        eps[0] += deltapsiDet;
        eps[1] += deltapsiRefA;
        eps[2] += deltapsiRefB;

        qx_shifted[0] = coll.qvecRe()[4 * detId + harmInd] * std::cos(deltapsiDet) - coll.qvecIm()[4 * detId + harmInd] * std::sin(deltapsiDet);
        qy_shifted[0] = coll.qvecRe()[4 * detId + harmInd] * std::sin(deltapsiDet) + coll.qvecIm()[4 * detId + harmInd] * std::cos(deltapsiDet);
        qx_shifted[1] = coll.qvecRe()[4 * refAId + harmInd] * std::cos(deltapsiRefA) - coll.qvecIm()[4 * refAId + harmInd] * std::sin(deltapsiRefA);
        qy_shifted[1] = coll.qvecRe()[4 * refAId + harmInd] * std::sin(deltapsiRefA) + coll.qvecIm()[4 * refAId + harmInd] * std::cos(deltapsiRefA);
        qx_shifted[2] = coll.qvecRe()[4 * refBId + harmInd] * std::cos(deltapsiRefB) - coll.qvecIm()[4 * refBId + harmInd] * std::sin(deltapsiRefB);
        qy_shifted[2] = coll.qvecRe()[4 * refBId + harmInd] * std::sin(deltapsiRefB) + coll.qvecIm()[4 * refBId + harmInd] * std::cos(deltapsiRefB);
      }

      float resNumA = helperEP.GetResolution(eps[0], eps[1], i + 2);
      float resNumB = helperEP.GetResolution(eps[0], eps[2], i + 2);
      float resDenom = helperEP.GetResolution(eps[1], eps[2], i + 2);

      epFlowHistograms.fill(HIST("EpDet"), i + 2, cent, eps[0]);
      epFlowHistograms.fill(HIST("EpRefA"), i + 2, cent, eps[1]);
      epFlowHistograms.fill(HIST("EpRefB"), i + 2, cent, eps[2]);

      epFlowHistograms.fill(HIST("EpResDetRefA"), i + 2, cent, resNumA);
      epFlowHistograms.fill(HIST("EpResDetRefB"), i + 2, cent, resNumB);
      epFlowHistograms.fill(HIST("EpResRefARefB"), i + 2, cent, resDenom);

      epFlowHistograms.fill(HIST("EpResQvecDetRefAxx"), i + 2, cent, qx_shifted[0] * qx_shifted[1] + qy_shifted[0] * qy_shifted[1]);
      epFlowHistograms.fill(HIST("EpResQvecDetRefAxy"), i + 2, cent, qx_shifted[1] * qy_shifted[0] - qx_shifted[0] * qy_shifted[1]);
      epFlowHistograms.fill(HIST("EpResQvecDetRefBxx"), i + 2, cent, qx_shifted[0] * qx_shifted[2] + qy_shifted[0] * qy_shifted[2]);
      epFlowHistograms.fill(HIST("EpResQvecDetRefBxy"), i + 2, cent, qx_shifted[2] * qy_shifted[0] - qx_shifted[0] * qy_shifted[2]);
      epFlowHistograms.fill(HIST("EpResQvecRefARefBxx"), i + 2, cent, qx_shifted[1] * qx_shifted[2] + qy_shifted[1] * qy_shifted[2]);
      epFlowHistograms.fill(HIST("EpResQvecRefARefBxy"), i + 2, cent, qx_shifted[2] * qy_shifted[1] - qx_shifted[1] * qy_shifted[2]);

      for (const auto& track : tracks) {
        if (trackSel(track))
          continue;

        if (cfgEffCor) {
          weight = getEfficiencyCorrection(effMap, track.eta(), track.pt(), cent, coll.posZ());
        }

        float vn = std::cos((i + 2) * (track.phi() - eps[0]));
        float vnSin = std::sin((i + 2) * (track.phi() - eps[0]));

        epFlowHistograms.fill(HIST("vncos"), i + 2, cent, track.pt(), vn, weight);
        epFlowHistograms.fill(HIST("vnsin"), i + 2, cent, track.pt(), vnSin, weight);

        epFlowHistograms.fill(HIST("SPvnxx"), i + 2, cent, track.pt(), (std::cos(track.phi() * static_cast<float>(i + 2)) * qx_shifted[0] + std::sin(track.phi() * static_cast<float>(i + 2)) * qy_shifted[0]), weight);
        epFlowHistograms.fill(HIST("SPvnxy"), i + 2, cent, track.pt(), (std::sin(track.phi() * static_cast<float>(i + 2)) * qx_shifted[0] - std::cos(track.phi() * static_cast<float>(i + 2)) * qy_shifted[0]), weight);
      }
    }
  }
  PROCESS_SWITCH(jEPFlowAnalysis, processDefault, "default process", true);

  void processMCRec(MyCollisionsMC::iterator const& coll, MyTracksMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!coll.has_mcCollision()) {
      return;
    }

    if (cfgAddEvtSel) {
      if (std::abs(coll.posZ()) > cfgVertexZ)
        return;
      switch (cfgEvtSel) {
        case 0: // Sel8
          if (!coll.sel8())
            return;
          break;
        case 1: // PbPb standard
          if (!coll.sel8() || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 2: // PbPb with pileup
          if (!coll.sel8() || !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ||
              !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 3: // Small systems (OO, NeNe, pp)
          if (!coll.sel8() || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        default:
          LOGF(warning, "Event selection flag was not found, continuing without basic event selections!\n");
      }
      // Check occupancy
      if (coll.trackOccupancyInTimeRange() > cfgMaxOccupancy || coll.trackOccupancyInTimeRange() < cfgMinOccupancy)
        return;
    }

    float cent = coll.centFT0C();

    if (cfgEffCor) {
      auto bc = coll.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        effMap = ccdb->getForTimeStamp<THnT<float>>(cfgEffCorDir, bc.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    for (auto trk : tracks) {
      if (!trk.has_mcParticle()) {
        continue;
      }

      if (trackSel(trk)) {
        continue;
      }

      epFlowHistograms.fill(HIST("MC/hPartRec"), cent, coll.posZ(), trk.eta(), trk.phi(), trk.pt());
      auto mctrk = trk.mcParticle();
      if (mctrk.isPhysicalPrimary()) {
        epFlowHistograms.fill(HIST("MChPartRecPr"), cent, coll.posZ(), trk.eta(), trk.phi(), trk.pt());
      }
    }
  }
  PROCESS_SWITCH(jEPFlowAnalysis, processMCRec, "process for MC", false);

  void processMCGen(MyCollisionsMC::iterator const& coll, aod::McParticles const& mcParticles, aod::McCollisions const&)
  {
    if (!coll.has_mcCollision())
      return;
    const auto mcColl = coll.mcCollision();

    if (cfgAddEvtSel) {
      if (std::abs(mcColl.posZ()) > cfgVertexZ) {
        return;
      }
    }

    float cent = coll.centFT0C();

    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > cfgTrackCuts.cfgEtaMax)
        continue;

      auto* p = TDatabasePDG::Instance()->GetParticle(mcParticle.pdgCode());
      if (p) {
        if (std::abs(p->Charge()) < 1e-1) {
          continue;
        }
      }

      if (!mcParticle.isPhysicalPrimary())
        continue;

      epFlowHistograms.fill(HIST("MC/hPartGen"), cent, mcColl.posZ(), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
    }
  }
  PROCESS_SWITCH(jEPFlowAnalysis, processMCGen, "process for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jEPFlowAnalysis>(cfgc)};
}
