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
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"

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

//using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr>;
using TracksCompleteIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCLfPi, aod::pidTPCLfDe>;
using V0full = soa::Join<aod::V0Datas, aod::V0Covs>;

struct hypertritonAnalysis {

  HistogramRegistry registry{
    "registry",
    {
      {"h2dMassHypertriton", "h2dMassHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h2dMassAntiHypertriton", "h2dMassAntiHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},

      //Basic QA histogram
      {"h3dPtVsMassHyVsDCAxy", "h3dPtVsMassHyVsDCAxy", {HistType::kTH3F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 2.900f, 3.100f, "Inv. Mass (GeV/c^{2})"}, {50, 0.0f, 1.0f, "abs(dcaxy)"}}}},
      {"h3dPtVsMassAHyVsDCAxy", "h3dPtVsMassAHyVsDCAxy", {HistType::kTH3F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 2.900f, 3.100f, "Inv. Mass (GeV/c^{2})"}, {50, 0.0f, 1.0f, "abs(dcaxy)"}}}},
      
      {"hEventSelection", "hEventSelection", {HistType::kTH1F, {{3, -0.5f, 2.5f}}}},
      {"V0loopFiltersCounts", "V0loopFiltersCounts", {HistType::kTH1F, {{10, -0.5f, 9.5f}}}},
    },
  };
  
  void init(InitContext const&)
  {
    resetHistos();
    
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "Sel8 cut");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
  }
  
  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dca3hetopv{"dca3hetopv", .0, "DCA Neg To PV"};
  Configurable<float> dcapiontopv{"dcapiontopv", .2, "DCA Pos To PV"};
  Configurable<float> dcahyptopv{"dcahyptopv", .2, "DCA Hypertriton To PV"};
  Configurable<float> v0radius{"v0radius", 1.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> deuteronUpperTPC{"deuteronUpperTPC", 8, "deuteronUpperTPC"};
  Configurable<float> piondEdx{"piondEdx", 5, "piondEdx"};
  
  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};
  
  // Material correction to use when propagating hypertriton: none
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::track::TrackParCov lHyTrack;
  
  Filter preFilterV0 = aod::v0data::dcaV0daughters < dcav0dau;
  
  enum anastep { kHypAll = 0,
                 kHypRadiusCosPA,
                 kHypRapidity,
                 kHypTPCdEdx,
                 kHypDauDCAtoPV,
                 kHypDCAtoPV,
                 kHypAllSteps};
  
  enum evselstep { kEvSelAll = 0,
                   kEvSelBool,
                   kEvSelVtxZ,
                   kEvSelAllSteps};
  
  // Helper to do bookkeeping and late filling of QA histos
  std::array<long, kHypAllSteps> stats;
  std::array<long, kEvSelAllSteps> evselstats;

  void resetHistos()
  {
    for (Int_t ii = 0; ii < kHypAllSteps; ii++)
      stats[ii] = 0;
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      evselstats[ii] = 0;
  }

  void fillHistos()
  {
    for (Int_t ii = 0; ii < kHypAllSteps; ii++)
      registry.fill(HIST("V0loopFiltersCounts"), ii, stats[ii]);
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      registry.fill(HIST("hEventSelection"), ii, evselstats[ii]);
  }
  
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<V0full> const& fullV0s, TracksCompleteIU const& tracks)
  {
    const double lHypertritonMass = 2.99131;
    const double l3HeMass = 2.80923;
    
    gpu::gpustd::array<float, 2> dcaInfo;
    
    evselstats[kEvSelAll]++;
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    evselstats[kEvSelBool]++;
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    evselstats[kEvSelVtxZ]++;
    
    for (auto& v0 : fullV0s) {
      stats[kHypAll]++;
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        stats[kHypRadiusCosPA]++;
        //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
        // Hypertriton
        if (TMath::Abs(v0.yHypertriton()) < rapidity) {
          stats[kHypRapidity]++;
          if (v0.posTrack_as<TracksCompleteIU>().tpcNSigmaDe() > deuteronUpperTPC&&
              TMath::Abs(v0.negTrack_as<TracksCompleteIU>().tpcNSigmaPi()) > piondEdx ) {
            stats[kHypTPCdEdx]++;
            if( TMath::Abs( v0.posTrack_as<TracksCompleteIU>().dcaXY() ) < dca3hetopv&&
               TMath::Abs( v0.negTrack_as<TracksCompleteIU>().dcaXY() ) < dcapiontopv){
              stats[kHypDauDCAtoPV]++;
              // Set up covariance matrices (should in fact be optional)
              std::array<float, 21> covV = {0.};
              constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
              for (int i = 0; i < 6; i++) {
                covV[MomInd[i]] = v0.momentumCovMat()[i]; //fixme: rigidity <-> momentum not taken into account in cov
                covV[i] = v0.positionCovMat()[i];
              }
              lHyTrack = o2::track::TrackParCov(
                {v0.x(), v0.y(), v0.z()},
                {2.0f*v0.pxpos() + v0.pxneg(), 2.0f*v0.pypos() + v0.pyneg(), 2.0f*v0.pzpos() + v0.pzneg()},
                covV, +1, true);

              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
              
              if( TMath::Abs( dcaInfo[0] ) < dcahyptopv ){
                stats[kHypDCAtoPV]++;
                registry.fill(HIST("h2dMassHypertriton"), v0.ptHypertriton(), v0.mHypertriton());
                registry.fill(HIST("h3dPtVsMassHyVsDCAxy"), v0.ptHypertriton(), v0.mHypertriton(), TMath::Abs( dcaInfo[0] ) );
              }
            }
          }
        }
        //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
        // AntiHypertriton
        if (TMath::Abs(v0.yAntiHypertriton()) < rapidity) {
          if (v0.negTrack_as<TracksCompleteIU>().tpcNSigmaDe() > deuteronUpperTPC&&
              TMath::Abs(v0.posTrack_as<TracksCompleteIU>().tpcNSigmaPi()) > piondEdx ) {
            stats[kHypTPCdEdx]++;
            if( TMath::Abs( v0.posTrack_as<TracksCompleteIU>().dcaXY() ) < dcapiontopv&&
               TMath::Abs( v0.negTrack_as<TracksCompleteIU>().dcaXY() ) < dca3hetopv){
              stats[kHypDauDCAtoPV]++;
              // Set up covariance matrices (should in fact be optional)
              std::array<float, 21> covV = {0.};
              constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
              for (int i = 0; i < 6; i++) {
                covV[MomInd[i]] = v0.momentumCovMat()[i]; //fixme: rigidity <-> momentum not taken into account in cov
                covV[i] = v0.positionCovMat()[i];
              }
              lHyTrack = o2::track::TrackParCov(
                {v0.x(), v0.y(), v0.z()},
                {v0.pxpos() + 2.0f*v0.pxneg(), v0.pypos() + 2.0f*v0.pyneg(), v0.pzpos() + 2.0f*v0.pzneg()},
                covV, -1, true);

              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
              
              if( TMath::Abs( dcaInfo[0] ) < dcahyptopv ){
                stats[kHypDCAtoPV]++;
                registry.fill(HIST("h2dMassAntiHypertriton"), v0.ptAntiHypertriton(), v0.mAntiHypertriton());
                registry.fill(HIST("h3dPtVsMassAHyVsDCAxy"), v0.ptAntiHypertriton(), v0.mAntiHypertriton(), TMath::Abs( dcaInfo[0] ) );
              }
            }
          }
        }
      }
    }//end v0 loop
    fillHistos();
    resetHistos();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysis>(cfgc)};
}
