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
/// \since Jul 2024
/// \brief task for flow measurement
/// \file jEPFlowAnalysis.cxx

#include "JEPFlowAnalysis.h"

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

struct JEPFlowAnalysis {

  HistogramRegistry epFlowHistograms{"EPFlow", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  JEPFlowAnalysis epAnalysis;
  EventPlaneHelper helperEP;
  FlowJHistManager histManager;
  bool debug = kFALSE;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Set Configurables here
  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "Use event selection"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "Total number of detectors in qVectorsTable"};

  Configurable<int> cfgNmode{"cfgNmode", 1, "the number of modulations to be calculated"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};
  Configurable<bool> cfgSPmethod{"cfgSPmethod", false, "flag for scalar product"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);

  int detId;
  int refAId;
  int refBId;
  int harmInd;

  int currentRunNumber = -999;
  int lastRunNumber = -999;

  int initflow = 2;
  int nshift = 10;
  int nsystem = 3;
  int nqvec = 4;
  int nstep = 3;

  std::vector<TProfile3D*> shiftprofile{};
  std::string fullCCDBShiftCorrPath;

  template <typename T>
  int getDetId(const T& name)
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
    detId = getDetId(cfgDetName);
    refAId = getDetId(cfgRefAName);
    refBId = getDetId(cfgRefBName);

    epAnalysis.SetHistRegistry(&epFlowHistograms);
    epAnalysis.CreateHistograms();
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
        for (int i = initflow; i < initflow + cfgNmode; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }

    for (int i = initflow; i < initflow + cfgNmode; i++) {                  // loop over different harmonic orders
      harmInd = cfgnTotalSystem * nqvec * (i - initflow) + nstep; // harmonic index to access corresponding Q-vector as all Q-vectors are in same vector
      eps[0] = helperEP.GetEventPlane(coll.qvecRe()[detId + harmInd], coll.qvecIm()[detId + harmInd], i);
      eps[1] = helperEP.GetEventPlane(coll.qvecRe()[refAId + harmInd], coll.qvecIm()[refAId + harmInd], i);
      eps[2] = helperEP.GetEventPlane(coll.qvecRe()[refBId + harmInd], coll.qvecIm()[refBId + harmInd], i);

      auto deltapsiDet = 0.0;
      auto deltapsiRefA = 0.0;
      auto deltapsiRefB = 0.0;

      float weight = 1.0;

      if (cfgShiftCorr) {
        for (int ishift = 1; ishift <= nshift; ishift++) {
          auto coeffshiftxDet = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 0.5, ishift - 0.5));
          auto coeffshiftyDet = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 1.5, ishift - 0.5));
          auto coeffshiftxRefA = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 2.5, ishift - 0.5));
          auto coeffshiftyRefA = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 3.5, ishift - 0.5));
          auto coeffshiftxRefB = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 4.5, ishift - 0.5));
          auto coeffshiftyRefB = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 5.5, ishift - 0.5)); // currently only FT0C/TPCpos/TPCneg

          deltapsiDet += ((1 / (1.0 * ishift)) * (-coeffshiftxDet * std::cos(ishift * static_cast<float>(i) * eps[0]) + coeffshiftyDet * std::sin(ishift * static_cast<float>(i) * eps[0])));
          deltapsiRefA += ((1 / (1.0 * ishift)) * (-coeffshiftxRefA * std::cos(ishift * static_cast<float>(i) * eps[1]) + coeffshiftyRefA * std::sin(ishift * static_cast<float>(i) * eps[1])));
          deltapsiRefB += ((1 / (1.0 * ishift)) * (-coeffshiftxRefB * std::cos(ishift * static_cast<float>(i) * eps[2]) + coeffshiftyRefB * std::sin(ishift * static_cast<float>(i) * eps[2])));
        }

        eps[0] += deltapsiDet;
        eps[1] += deltapsiRefA;
        eps[2] += deltapsiRefB;
      }

      if (cfgSPmethod)
        weight *= std::sqrt(std::pow(coll.qvecRe()[detId + harmInd], 2) + std::pow(coll.qvecIm()[detId + harmInd], 2));

      float resNumA = helperEP.GetResolution(eps[0], eps[1], i);
      float resNumB = helperEP.GetResolution(eps[0], eps[2], i);
      float resDenom = helperEP.GetResolution(eps[1], eps[2], i);
      epAnalysis.FillResolutionHistograms(cent, static_cast<float>(i), resNumA, resNumB, resDenom);

      for (uint j = 0; j < nsystem; j++) { // loop over detectors used
        for (const auto& track : tracks) {
          float vn = std::cos((i) * (track.phi() - eps[j]));
          float vn_sin = std::sin((i) * (track.phi() - eps[j]));
          epAnalysis.FillVnHistograms(i, cent, static_cast<float>(j + 1), track.pt(), vn * weight, vn_sin * weight);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JEPFlowAnalysis>(cfgc)};
}
