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
// Test hypertriton task
// =====================
//
// First rough code for hypertriton checks
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

#using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr>;
using TracksCompleteIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, , aod::pidTPCPr>;

struct hypertritonAnalysis {

  HistogramRegistry registry{
    "registry",
    {
      {"h2dMassHypertriton", "h2dMassHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h2dMassAntiHypertriton", "h2dMassAntiHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hEventSelection", "hEventSelection", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"V0loopFiltersCounts", "V0loopFiltersCounts", {HistType::kTH1F, {{11, 0.0f, 11.0f}}}},
    },
  };
  
  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(1, "V0 Candidates");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(2, "V0Radius and CosPA");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(4, "Lambda Rapidity");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(5, "Lambda lifetime cut");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(6, "Lambda TPC PID cut");
    
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "Sel8 cut");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    
    if (doprocessRun3 && doprocessRun2) {
      LOGF(fatal, "processRun3 and processRun2 are both set to true; try again with only one of them set to true");
    }
    if (!doprocessRun3 && !doprocessRun2) {
      LOGF(fatal, "processRun3 nor processRun2 are both set to false; try again with only one of them set to false");
    }
  }
  
  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  
  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};
  
  const double lHypertritonMass = 2.99131;
  const double l3HeMass = 2.80923;
  
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;
  
  // void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s) //for now CentV0M info is not available for run 3 pp
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, TracksCompleteIU const& tracks)
  {
    registry.fill(HIST("hEventSelection"), 0.5);
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventSelection"), 1.5);
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    registry.fill(HIST("hEventSelection"), 2.5);
    
    for (auto& v0 : fullV0s) {
      // FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("V0loopFiltersCounts"), 0.5);
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("V0loopFiltersCounts"), 1.5);
        if (TMath::Abs(v0.yHypertriton()) < rapidity) {
          registry.fill(HIST("V0loopFiltersCounts"), 3.5);
          // Hypertriton
          if (v0.posTrack_as<TracksCompleteIU>().tpcNSigmaPr() > TpcPidNsigmaCut) { // NEEDS FIXING and no cut on K0S
            registry.fill(HIST("V0loopFiltersCounts"), 5.5);
            registry.fill(HIST("h2dMassHypertriton"), v0.ptHypertriton(), v0.mHypertriton());
          }
        }
        if (TMath::Abs(v0.yAntiHypertriton()) < rapidity) {
          if (v0.negTrack_as<TracksCompleteIU>().tpcNSigmaPr() > TpcPidNsigmaCut) { // NEEDS FIXING
            registry.fill(HIST("V0loopFiltersCounts"), 5.5);
            registry.fill(HIST("h2dMassAntiHypertriton"), v0.ptAntiHypertriton(), v0.mAntiHypertriton());
          }
          
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysis>(cfgc)};
}
