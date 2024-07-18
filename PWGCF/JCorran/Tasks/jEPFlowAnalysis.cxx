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

#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
// #include "Common/DataModel/Centrality.h"
// #include "ReconstructionDataFormats/V0.h"

// #include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/EventPlaneHelper.h"

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

  // Set Configurables here
  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "Use event selection"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);

  int DetId;
  int RefAId;
  int RefBId;

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


  void init(InitContext const&) {
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);

    epAnalysis.SetHistRegistry(&EPFlowHistograms);
    epAnalysis.CreateHistograms();
  }

  Float_t ResolutionByEP(const Float_t EP_A, const Float_t EP_B, const Float_t EP_C, const Int_t ind) {
    Float_t resNumA = helperEP.GetResolution(EP_A, EP_B, ind);
    Float_t resNumB = helperEP.GetResolution(EP_A, EP_C, ind);
    Float_t resDen = helperEP.GetResolution(EP_B, EP_C, ind);
    if (debug )printf("EP_A: %.5f, EP_B: %.5f, ResAB: %.5f, Ind: %d\n", EP_A, EP_B, resNumA, ind);
    Float_t resolution = TMath::Sqrt(resNumA*resNumB/resDen);
    return resolution;
  }

  void process(MyCollisions::iterator const& coll, soa::Filtered<MyTracks> const& tracks)
  {
    if (cfgAddEvtSel && (!coll.sel8() || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)
        || !coll.selection_bit(aod::evsel::kNoSameBunchPileup))) return;
    
    Int_t cBin = histManager.GetCentBin(coll.cent());
    Float_t EPs[3] = {0.};
    Float_t vn[3][3] = {{0.}};
    for (uint i = 0; i<3; i++) {
      EPs[0] = helperEP.GetEventPlane(coll.qvecRe()[DetId+3], coll.qvecRe()[DetId+3], i+2); 
      EPs[1] = helperEP.GetEventPlane(coll.qvecRe()[RefAId+3], coll.qvecRe()[RefAId+3], i+2);
      EPs[2] = helperEP.GetEventPlane(coll.qvecRe()[RefBId+3], coll.qvecRe()[RefBId+3], i+2);

      for (uint j=0; j<3; j++) {
        Float_t resolution = ResolutionByEP(EPs[j%3],EPs[(j+1)%3],EPs[(j+2)%3], i+2);
      
        if (debug) printf("Ind: %u, Res: %.5f\n", j, resolution);
        Float_t sumCos=0;
        for (auto& track: tracks) {
          sumCos += TMath::Cos((i+2)*(track.phi()-EPs[j]));
        }

        vn[i][j] = sumCos/(float(tracks.size())*resolution);
      }

    }
    for (uint j=0; j<3; j++) epAnalysis.FillHistograms(cBin, float(j+1), vn[0][j], vn[1][j], vn[2][j]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jEPFlowAnalysis>(cfgc)};
}
