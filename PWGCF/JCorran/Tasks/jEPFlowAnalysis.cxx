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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
using MyTracks = aod::Tracks;

struct jEPFlowAnalysis {

  HistogramRegistry epFlowHistograms{"EPFlow", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  EventPlaneHelper helperEP;
  FlowJHistManager histManager;
  bool debug = kFALSE;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Set Configurables here
  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "Use event selection"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "Total number of detectors in qVectorsTable"};
  Configurable<int> cfgnMode{"cfgnMode", 1, "the number of modulations"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};
  Configurable<bool> cfgSPmethod{"cfgSPmethod", false, "flag for scalar product"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {100, 0, 100}, ""};
  ConfigurableAxis cfgAxisPt{"cfgAxisPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 30.0, 50.0, 70.0, 100.0}, ""};
  ConfigurableAxis cfgAxisCos{"cfgAxisCos", {102, -1.02, 1.02}, ""};

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);

  int detId;
  int refAId;
  int refBId;
  int harmInd;

  int currentRunNumber = -999;
  int lastRunNumber = -999;

  std::vector<TProfile3D*> shiftprofile{};
  std::string fullCCDBShiftCorrPath;

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

  void init(InitContext const&)
  {
    detId = getdetId(cfgDetName);
    refAId = getdetId(cfgRefAName);
    refBId = getdetId(cfgRefBName);

    AxisSpec axisMod{cfgnMode, 2., cfgnMode + 2.};
    AxisSpec axisEvtPl{360, -constants::math::PI * 1.1, constants::math::PI * 1.1};

    AxisSpec axisCent{cfgAxisCent, "cent"};
    AxisSpec axisPt{cfgAxisPt, "pT"};
    AxisSpec axisCos{cfgAxisCos, "cos"};

    epFlowHistograms.add("EpDet", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpRefA", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpRefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});

    epFlowHistograms.add("EpResDetRefA", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpResDetRefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});
    epFlowHistograms.add("EpResRefARefB", "", {HistType::kTH3F, {axisMod, axisCent, axisEvtPl}});

    epFlowHistograms.add("vncos", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisCos}});
    epFlowHistograms.add("vnsin", "", {HistType::kTHnSparseF, {axisMod, axisCent, axisPt, axisCos}});
  }

  void process(MyCollisions::iterator const& coll, soa::Filtered<MyTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    if (cfgAddEvtSel && (!coll.sel8() || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup)))
      return;

    float cent = coll.cent();
    epFlowHistograms.fill(HIST("FullCentrality"), cent);
    float eps[3] = {0.};

    if (cfgShiftCorr) {
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

    for (int i = 0; i < cfgnMode; i++) {       // loop over different harmonic orders
      harmInd = cfgnTotalSystem * 4 * (i) + 3; // harmonic index to access corresponding Q-vector as all Q-vectors are in same vector
      eps[0] = helperEP.GetEventPlane(coll.qvecRe()[detId + harmInd], coll.qvecIm()[detId + harmInd], i + 2);
      eps[1] = helperEP.GetEventPlane(coll.qvecRe()[refAId + harmInd], coll.qvecIm()[refAId + harmInd], i + 2);
      eps[2] = helperEP.GetEventPlane(coll.qvecRe()[refBId + harmInd], coll.qvecIm()[refBId + harmInd], i + 2);

      auto deltapsiDet = 0.0;
      auto deltapsiRefA = 0.0;
      auto deltapsiRefB = 0.0;

      float weight = 1.0;

      if (cfgShiftCorr) {
        constexpr int kShiftBins = 10;
        for (int ishift = 1; ishift <= kShiftBins; ishift++) {
          auto coeffshiftxDet = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 0.5, ishift - 0.5));
          auto coeffshiftyDet = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 1.5, ishift - 0.5));
          auto coeffshiftxRefA = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 2.5, ishift - 0.5));
          auto coeffshiftyRefA = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 3.5, ishift - 0.5));
          auto coeffshiftxRefB = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 4.5, ishift - 0.5));
          auto coeffshiftyRefB = shiftprofile.at(i)->GetBinContent(shiftprofile.at(i)->FindBin(cent, 5.5, ishift - 0.5)); // currently only FT0C/TPCpos/TPCneg

          deltapsiDet += ((1 / (1.0 * ishift)) * (-coeffshiftxDet * std::cos(ishift * static_cast<float>(i + 2) * eps[0]) + coeffshiftyDet * std::sin(ishift * static_cast<float>(i + 2) * eps[0])));
          deltapsiRefA += ((1 / (1.0 * ishift)) * (-coeffshiftxRefA * std::cos(ishift * static_cast<float>(i + 2) * eps[1]) + coeffshiftyRefA * std::sin(ishift * static_cast<float>(i + 2) * eps[1])));
          deltapsiRefB += ((1 / (1.0 * ishift)) * (-coeffshiftxRefB * std::cos(ishift * static_cast<float>(i + 2) * eps[2]) + coeffshiftyRefB * std::sin(ishift * static_cast<float>(i + 2) * eps[2])));
        }

        eps[0] += deltapsiDet;
        eps[1] += deltapsiRefA;
        eps[2] += deltapsiRefB;
      }

      if (cfgSPmethod)
        weight *= std::sqrt(std::pow(coll.qvecRe()[detId + harmInd], 2) + std::pow(coll.qvecIm()[detId + harmInd], 2));

      float resNumA = helperEP.GetResolution(eps[0], eps[1], i + 2);
      float resNumB = helperEP.GetResolution(eps[0], eps[2], i + 2);
      float resDenom = helperEP.GetResolution(eps[1], eps[2], i + 2);

      epFlowHistograms.fill(HIST("EpDet"), i + 2, cent, eps[0]);
      epFlowHistograms.fill(HIST("EpRefA"), i + 2, cent, eps[1]);
      epFlowHistograms.fill(HIST("EpRefB"), i + 2, cent, eps[2]);

      epFlowHistograms.fill(HIST("EpResDetRefA"), i + 2, cent, resNumA);
      epFlowHistograms.fill(HIST("EpResDetRefB"), i + 2, cent, resNumB);
      epFlowHistograms.fill(HIST("EpResRefARefB"), i + 2, cent, resDenom);

      for (int j = 0; j < cfgnMode; j++) { // loop over detectors used
        for (const auto& track : tracks) {
          float vn = std::cos((i + 2) * (track.phi() - eps[j]));
          float vnSin = std::sin((i + 2) * (track.phi() - eps[j]));

          epFlowHistograms.fill(HIST("vncos"), i + 2, cent, track.pt(), vn * weight);
          epFlowHistograms.fill(HIST("vnsin"), i + 2, cent, track.pt(), vnSin * weight);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jEPFlowAnalysis>(cfgc)};
}
