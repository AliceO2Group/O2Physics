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

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "FlowJHistManager.h"
#include "JEPFlowAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
using MyTracks = aod::Tracks;

struct jEPFlowAnalysis {

  HistogramRegistry EPFlowHistograms{"EPFlow", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  JEPFlowAnalysis epAnalysis;
  EventPlaneHelper helperEP;
  FlowJHistManager histManager;
  Bool_t debug = kFALSE;
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

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};
  Configurable<bool> cfgSPmethod{"cfgSPmethod", false, "flag for scalar product"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);

  int DetId;
  int RefAId;
  int RefBId;
  int harmInd;

  int currentRunNumber = -999;
  int lastRunNumber = -999;

  std::vector<TProfile3D*> shiftprofile{};
  std::string fullCCDBShiftCorrPath;

  template <typename T>
  int GetDetId(const T& name)
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
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);

    epAnalysis.SetHistRegistry(&EPFlowHistograms);
    epAnalysis.CreateHistograms();
  }

  void process(MyCollisions::iterator const& coll, soa::Filtered<MyTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    if (cfgAddEvtSel && (!coll.sel8() || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !coll.selection_bit(aod::evsel::kNoSameBunchPileup)))
      return;

    Float_t cent = coll.cent();
    EPFlowHistograms.fill(HIST("FullCentrality"), cent);
    Float_t EPs[3] = {0.};

    if (cfgShiftCorr) {
      auto bc = coll.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        shiftprofile.clear();
        for (int i = 2; i < 5; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }

    for (int i = 2; i < 5; i++) {                  // loop over different harmonic orders
      harmInd = cfgnTotalSystem * 4 * (i - 2) + 3; // harmonic index to access corresponding Q-vector as all Q-vectors are in same vector
      EPs[0] = helperEP.GetEventPlane(coll.qvecRe()[DetId + harmInd], coll.qvecIm()[DetId + harmInd], i);
      EPs[1] = helperEP.GetEventPlane(coll.qvecRe()[RefAId + harmInd], coll.qvecIm()[RefAId + harmInd], i);
      EPs[2] = helperEP.GetEventPlane(coll.qvecRe()[RefBId + harmInd], coll.qvecIm()[RefBId + harmInd], i);

      auto deltapsiDet = 0.0;
      auto deltapsiRefA = 0.0;
      auto deltapsiRefB = 0.0;

      float weight = 1.0;

      if (cfgShiftCorr) {
        for (int ishift = 1; ishift <= 10; ishift++) {
          auto coeffshiftxDet = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 0.5, ishift - 0.5));
          auto coeffshiftyDet = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 1.5, ishift - 0.5));
          auto coeffshiftxRefA = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 2.5, ishift - 0.5));
          auto coeffshiftyRefA = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 3.5, ishift - 0.5));
          auto coeffshiftxRefB = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 4.5, ishift - 0.5));
          auto coeffshiftyRefB = shiftprofile.at(i - 2)->GetBinContent(shiftprofile.at(i - 2)->FindBin(cent, 5.5, ishift - 0.5)); //currently only FT0C/TPCpos/TPCneg

          deltapsiDet += ((1 / (1.0 * ishift)) * (-coeffshiftxDet * TMath::Cos(ishift * static_cast<float>(i) * EPs[0]) + coeffshiftyDet * TMath::Sin(ishift * static_cast<float>(i) * EPs[0])));
          deltapsiRefA += ((1 / (1.0 * ishift)) * (-coeffshiftxRefA * TMath::Cos(ishift * static_cast<float>(i) * EPs[1]) + coeffshiftyRefA * TMath::Sin(ishift * static_cast<float>(i) * EPs[1])));
          deltapsiRefB += ((1 / (1.0 * ishift)) * (-coeffshiftxRefB * TMath::Cos(ishift * static_cast<float>(i) * EPs[2]) + coeffshiftyRefB * TMath::Sin(ishift * static_cast<float>(i) * EPs[2])));
        }

        EPs[0] += deltapsiDet;
        EPs[1] += deltapsiRefA;
        EPs[2] += deltapsiRefB;
      }

      if (cfgSPmethod) weight *= sqrt(pow(coll.qvecRe()[DetId + harmInd], 2) + pow(coll.qvecIm()[DetId + harmInd], 2));

      Float_t resNumA = helperEP.GetResolution(EPs[0], EPs[1], i);
      Float_t resNumB = helperEP.GetResolution(EPs[0], EPs[2], i);
      Float_t resDenom = helperEP.GetResolution(EPs[1], EPs[2], i);
      epAnalysis.FillResolutionHistograms(cent, static_cast<float>(i), resNumA, resNumB, resDenom);
      for (uint j = 0; j < 3; j++) { // loop over detectors used
        for (auto& track : tracks) {
          Float_t vn = TMath::Cos((i) * (track.phi() - EPs[j]));
          Float_t vn_sin = TMath::Sin((i) * (track.phi() - EPs[j]));
          epAnalysis.FillVnHistograms(i, cent, static_cast<float>(j + 1), track.pt(), vn * weight, vn_sin * weight);
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
