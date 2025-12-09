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
// Strangeness reconstruction QA
// =============================
//
// Dedicated task to understand reconstruction
// Special emphasis on PV reconstruction when strangeness is present
// Tested privately, meant to be used on central MC productions now
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr>;
using TracksCompleteIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using TracksCompleteIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::McTrackLabels>;
using V0DataLabeled = soa::Join<aod::V0Datas, aod::McV0Labels>;
using CascMC = soa::Join<aod::CascDataExt, aod::McCascLabels>;
using TraCascMC = soa::Join<aod::TraCascDatas, aod::McTraCascLabels>;
using RecoedMCCollisions = soa::Join<aod::McCollisions, aod::McCollsExtra>;
using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;

struct preProcessMCcollisions {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  struct collisionStats {
    float nContribsWithTPC = 0;
    float nContribsWithTRD = 0;
    float nContribsWithTRDNoTOF = 0;
    float nContribsWithTOF = 0;
    float nContribsWithITS = 0;
    float covTrace = 0;
    float deltaXY = -999; // positive def
    float deltaZ = -999;
    float deltaT = -999;
    float time = -999;
    float timeResolution = -999;
    float globalBC = -999;
    float positionX = -999;
    float positionY = -999;
    float positionZ = -999;
  };

  template <typename T>
  std::vector<std::size_t> sort_indices(const std::vector<T>& v)
  {
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&v](std::size_t i1, std::size_t i2) { return v[i1] > v[i2]; });
    return idx;
  }

  void init(InitContext const&)
  {
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};
    const AxisSpec axisTrackCount{(int)50, -0.5f, +49.5f, ""};
    const AxisSpec axisContributors{(int)200, -0.5f, +199.5f, ""};
    const AxisSpec axisContributorsTRD{(int)200, -0.5f, +199.5f, "N_{contribs}^{TRD}"};
    const AxisSpec axisContributorsTOF{(int)200, -0.5f, +199.5f, "N_{contribs}^{TOF}"};
    const AxisSpec axisCovariance{(int)400, 0.0f, +0.1f, ""};
    const AxisSpec axisCovarianceTest{(int)400, -0.05f, +0.05f, ""};
    const AxisSpec axisTwenty{(int)20, -0.5f, +19.5f, ""};
    const AxisSpec axisTwentyWithNegative{(int)20, -10.5, +9.5f, ""};
    const AxisSpec axisTwentyFloatWithNegative{(int)100, -10.0f, +10.0f, ""};

    const AxisSpec axisPVResolutionXY{(int)400, 0.0f, +0.05f, ""};
    const AxisSpec axisPVResolutionZ{(int)400, -0.1f, +0.1f, ""};
    const AxisSpec axisPVResolutionT{(int)500, -20.0f, +20.0f, ""};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});
    histos.add("hNTimesCollWithXiRecoed", "hNTimesCollWithXiRecoed", kTH1F, {axisNTimesCollRecoed});

    // A trick to store more information, please
    histos.add("h2dTrackCounter", "h2dTrackCounter", kTH2D, {axisTrackCount, axisNTimesCollRecoed});
    histos.add("h2dTrackCounterWithXi", "h2dTrackCounterWithXi", kTH2D, {axisTrackCount, axisNTimesCollRecoed});

    // Number of contributor distributions - Y offset controls exact case
    histos.add("h2dNContributors", "h2dNContributors", kTH2D, {axisContributors, axisTwenty});
    histos.add("h2dNContributorsWithXi", "h2dNContributorsWithXi", kTH2D, {axisContributors, axisTwenty});

    // PV uncertainty estimate: trace of PV covariance matrix
    histos.add("hCyyTest", "hCyyTest", kTH1F, {axisCovarianceTest});
    histos.add("h2dCovarianceTrace", "h2dCovarianceTrace", kTH2D, {axisCovariance, axisTwenty});
    histos.add("h2dCovarianceTraceWithXi", "h2dCovarianceTraceWithXi", kTH2D, {axisCovariance, axisTwenty});

    // PV true error: from reco vs MC
    histos.add("h2dPVResolutionXY", "h2dPVResolutionXY", kTH2D, {axisPVResolutionXY, axisTwenty});
    histos.add("h2dPVResolutionZ", "h2dPVResolutionZ", kTH2D, {axisPVResolutionZ, axisTwenty});
    histos.add("h2dPVResolutionT", "h2dPVResolutionT", kTH2D, {axisPVResolutionT, axisTwenty});
    histos.add("h2dPVResolutionXYWithXi", "h2dPVResolutionXYWithXi", kTH2D, {axisPVResolutionXY, axisTwenty});
    histos.add("h2dPVResolutionZWithXi", "h2dPVResolutionZWithXi", kTH2D, {axisPVResolutionZ, axisTwenty});
    histos.add("h2dPVResolutionTWithXi", "h2dPVResolutionTWithXi", kTH2D, {axisPVResolutionT, axisTwenty});

    // Helper to decipher this histogram
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(1, "Recoed 1 time, 1st PV");      // size 1 = 0
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(2, "Recoed 2 times, Biggest PV"); // size 2 = 1
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(3, "Recoed 2 times, Smallest PV");
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(4, "Recoed 3 times, Biggest PV"); // size 3 = 3
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(5, "Recoed 3 times, Intermediate PV");
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(6, "Recoed 3 times, Smallest PV");
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(7, "Recoed 4 times, Biggest PV"); // size 4 = 6
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(8, "Recoed 4 times, 2nd Biggest PV");
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(9, "Recoed 4 times, 3rd Biggest PV");
    histos.get<TH2>(HIST("h2dNContributors"))->GetYaxis()->SetBinLabel(10, "Recoed 4 times, Smallest PV");

    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(1, "Recoed 1 time, 1st PV");      // size 1 = 0
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(2, "Recoed 2 times, Biggest PV"); // size 2 = 1
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(3, "Recoed 2 times, Smallest PV");
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(4, "Recoed 3 times, Biggest PV"); // size 3 = 3
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(5, "Recoed 3 times, Intermediate PV");
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(6, "Recoed 3 times, Smallest PV");
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(7, "Recoed 4 times, Biggest PV"); // size 4 = 6
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(8, "Recoed 4 times, 2nd Biggest PV");
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(9, "Recoed 4 times, 3rd Biggest PV");
    histos.get<TH2>(HIST("h2dNContributorsWithXi"))->GetYaxis()->SetBinLabel(10, "Recoed 4 times, Smallest PV");

    // Deltas between split vertices
    histos.add("h2dDeltaTimeVsDeltaZ", "h2dDeltaTimeVsDeltaZ", kTH2D, {axisTwentyFloatWithNegative, axisPVResolutionT});
    histos.add("h2dDeltaBCIndexVsDeltaZ", "h2dPVDeltaBCIndexVsDeltaZ", kTH2D, {axisTwentyFloatWithNegative, axisTwentyWithNegative});

    // Contributors correlation
    histos.add("h2dNContribCorrAll", "h2dNContribCorrAll", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribCorr1", "h2dNContribCorr1", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribCorr2a", "h2dNContribCorr2a", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribCorr2b", "h2dNContribCorr2b", kTH2D, {axisContributorsTRD, axisContributorsTOF});

    // Contributors correlation - TRD without TOF vs TOF
    histos.add("h2dNContribSpecialCorrAll", "h2dNContribSpecialCorrAll", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribSpecialCorr1", "h2dNContribSpecialCorr1", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribSpecialCorr2a", "h2dNContribSpecialCorr2a", kTH2D, {axisContributorsTRD, axisContributorsTOF});
    histos.add("h2dNContribSpecialCorr2b", "h2dNContribSpecialCorr2b", kTH2D, {axisContributorsTRD, axisContributorsTOF});
  }

  void processData(aod::Collisions const&)
  {
    // Dummy process
  }
  PROCESS_SWITCH(preProcessMCcollisions, processData, "Real data dummy", false);

  void processMC(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, TracksCompleteIUMC const& tracks, aod::McParticles const& mcParticles, aod::BCs const&)
  {
    int lNumberOfXi = 0;
    for (auto& mcp : mcParticles) {
      // mimic triggering strategy precisely
      if (TMath::Abs(mcp.eta()) < 0.8 && mcp.pdgCode() == 3312)
        lNumberOfXi++;
    }

    std::vector<collisionStats> collisionStatAggregator(collisions.size());
    std::vector<int> collisionNContribs;

    histos.fill(HIST("hNTimesCollRecoed"), collisions.size());
    if (lNumberOfXi > 0)
      histos.fill(HIST("hNTimesCollWithXiRecoed"), collisions.size());
    int lCollisionIndex = 0;
    for (auto& collision : collisions) {
      float cyy = TMath::Sign(1, collision.covYY()) * TMath::Sqrt(TMath::Abs(collision.covYY()));
      histos.fill(HIST("hCyyTest"), cyy); // check for bug
      collisionNContribs.emplace_back(collision.numContrib());
      collisionStatAggregator[lCollisionIndex].covTrace = TMath::Sqrt(TMath::Abs(collision.covXX()) + TMath::Abs(collision.covYY()) + TMath::Abs(collision.covZZ()));
      collisionStatAggregator[lCollisionIndex].deltaXY = TMath::Sqrt(TMath::Power(collision.posX() - mcCollision.posX(), 2) + TMath::Power(collision.posY() - mcCollision.posY(), 2));
      collisionStatAggregator[lCollisionIndex].deltaZ = collision.posZ() - mcCollision.posZ();
      collisionStatAggregator[lCollisionIndex].deltaT = (collision.collisionTime() - mcCollision.t()) / 1000;
      collisionStatAggregator[lCollisionIndex].time = collision.collisionTime();
      collisionStatAggregator[lCollisionIndex].timeResolution = collision.collisionTimeRes();
      collisionStatAggregator[lCollisionIndex].globalBC = collision.bc().globalBC();

      collisionStatAggregator[lCollisionIndex].positionX = collision.posX();
      collisionStatAggregator[lCollisionIndex].positionY = collision.posY();
      collisionStatAggregator[lCollisionIndex].positionZ = collision.posZ();

      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto& track : groupedTracks) {
        if (track.isPVContributor()) {
          if (track.hasITS())
            collisionStatAggregator[lCollisionIndex].nContribsWithITS++;
          if (track.hasTPC())
            collisionStatAggregator[lCollisionIndex].nContribsWithTPC++;
          if (track.hasTRD())
            collisionStatAggregator[lCollisionIndex].nContribsWithTRD++;
          if (track.hasTRD() && !track.hasTOF())
            collisionStatAggregator[lCollisionIndex].nContribsWithTRDNoTOF++;
          if (track.hasTOF())
            collisionStatAggregator[lCollisionIndex].nContribsWithTOF++;
        }
      }
      // Increment counter
      lCollisionIndex++;
    }
    // Collisions now exist, loop over them in NContribs order please
    lCollisionIndex = 0;
    auto sortedIndices = sort_indices(collisionNContribs);
    int lYAxisOffset = 0.5 * collisions.size() * (collisions.size() - 1);
    for (auto ic : sortedIndices) {
      histos.fill(HIST("h2dNContribCorrAll"), collisionStatAggregator[ic].nContribsWithTRD, collisionStatAggregator[ic].nContribsWithTOF);
      histos.fill(HIST("h2dNContribSpecialCorrAll"), collisionStatAggregator[ic].nContribsWithTRDNoTOF, collisionStatAggregator[ic].nContribsWithTOF);
      int lIndexBin = 7 * lCollisionIndex; // use offset to make plot much easier to read
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 0, collisions.size());
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 1, collisions.size(), collisionNContribs[ic]);
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 2, collisions.size(), collisionStatAggregator[ic].nContribsWithITS);
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 3, collisions.size(), collisionStatAggregator[ic].nContribsWithTPC);
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 4, collisions.size(), collisionStatAggregator[ic].nContribsWithTRD);
      histos.fill(HIST("h2dTrackCounter"), lIndexBin + 5, collisions.size(), collisionStatAggregator[ic].nContribsWithTOF);
      histos.fill(HIST("h2dNContributors"), collisionNContribs[ic], lYAxisOffset + lCollisionIndex);
      histos.fill(HIST("h2dCovarianceTrace"), collisionStatAggregator[ic].covTrace, lYAxisOffset + lCollisionIndex);
      histos.fill(HIST("h2dPVResolutionXY"), collisionStatAggregator[ic].deltaXY, lYAxisOffset + lCollisionIndex);
      histos.fill(HIST("h2dPVResolutionZ"), collisionStatAggregator[ic].deltaZ, lYAxisOffset + lCollisionIndex);
      histos.fill(HIST("h2dPVResolutionT"), collisionStatAggregator[ic].deltaT, lYAxisOffset + lCollisionIndex);
      if (lNumberOfXi > 0) {
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 0, collisions.size());
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 1, collisions.size(), collisionNContribs[ic]);
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 2, collisions.size(), collisionStatAggregator[ic].nContribsWithITS);
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 3, collisions.size(), collisionStatAggregator[ic].nContribsWithTPC);
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 4, collisions.size(), collisionStatAggregator[ic].nContribsWithTRD);
        histos.fill(HIST("h2dTrackCounterWithXi"), lIndexBin + 5, collisions.size(), collisionStatAggregator[ic].nContribsWithTOF);
        histos.fill(HIST("h2dNContributorsWithXi"), collisionNContribs[ic], lYAxisOffset + lCollisionIndex);
        histos.fill(HIST("h2dCovarianceTraceWithXi"), collisionStatAggregator[ic].covTrace, lYAxisOffset + lCollisionIndex);
        histos.fill(HIST("h2dPVResolutionXYWithXi"), collisionStatAggregator[ic].deltaXY, lYAxisOffset + lCollisionIndex);
        histos.fill(HIST("h2dPVResolutionZWithXi"), collisionStatAggregator[ic].deltaZ, lYAxisOffset + lCollisionIndex);
        histos.fill(HIST("h2dPVResolutionTWithXi"), collisionStatAggregator[ic].deltaT, lYAxisOffset + lCollisionIndex);
      }
      lCollisionIndex++;
    }
    // Now do investigation of split-into-two vertices, please - ordered
    if (collisions.size() == 1) {
      histos.fill(HIST("h2dNContribCorr1"), collisionStatAggregator[sortedIndices[0]].nContribsWithTRD, collisionStatAggregator[sortedIndices[0]].nContribsWithTOF);
      histos.fill(HIST("h2dNContribSpecialCorr1"), collisionStatAggregator[sortedIndices[0]].nContribsWithTRDNoTOF, collisionStatAggregator[sortedIndices[0]].nContribsWithTOF);
    }
    if (collisions.size() == 2) {
      float lDeltaTime = collisionStatAggregator[sortedIndices[0]].time - collisionStatAggregator[sortedIndices[1]].time;
      float lDeltaNBC = collisionStatAggregator[sortedIndices[0]].globalBC - collisionStatAggregator[sortedIndices[1]].globalBC;
      float lDeltaZ = collisionStatAggregator[sortedIndices[0]].positionZ - collisionStatAggregator[sortedIndices[1]].positionZ;
      histos.fill(HIST("h2dDeltaTimeVsDeltaZ"), lDeltaZ, lDeltaTime);
      histos.fill(HIST("h2dDeltaBCIndexVsDeltaZ"), lDeltaZ, lDeltaNBC);

      // Correaltion between number of contributors
      histos.fill(HIST("h2dNContribCorr2a"), collisionStatAggregator[sortedIndices[0]].nContribsWithTRD, collisionStatAggregator[sortedIndices[0]].nContribsWithTOF);
      histos.fill(HIST("h2dNContribCorr2b"), collisionStatAggregator[sortedIndices[1]].nContribsWithTRD, collisionStatAggregator[sortedIndices[1]].nContribsWithTOF);
      histos.fill(HIST("h2dNContribSpecialCorr2a"), collisionStatAggregator[sortedIndices[0]].nContribsWithTRDNoTOF, collisionStatAggregator[sortedIndices[0]].nContribsWithTOF);
      histos.fill(HIST("h2dNContribSpecialCorr2b"), collisionStatAggregator[sortedIndices[1]].nContribsWithTRDNoTOF, collisionStatAggregator[sortedIndices[1]].nContribsWithTOF);
    }
  }
  PROCESS_SWITCH(preProcessMCcollisions, processMC, "MC prepare", true);
};

struct straRecoStudy {
  // one to hold them all
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // binning matters
  Configurable<float> MaxPt{"MaxPt", 10, "maximum pT"};
  Configurable<int> NBinsPt{"NBinsPt", 100, "N bins"};
  Configurable<int> NBinsPtCoarse{"NBinsPtCoarse", 10, "N bins, coarse"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Selection criteria - compatible with core wagon autodetect
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.95, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  Configurable<float> maxV0Radius{"maxV0Radius", 200, "maxV0Radius"};
  Configurable<float> maxCascRadius{"maxCascRadius", 200, "maxCascRadius"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  Configurable<int> tpcmincrossedrows{"mincrossedrows", 70, "Minimum crossed rows"};
  Configurable<int> itsminclusters{"itsminclusters", 4, "Minimum ITS clusters"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  enum evselstep { kEvSelAll = 0,
                   kEvSelBool,
                   kEvSelVtxZ,
                   kEvSelAllSteps };

  std::array<long, kEvSelAllSteps> evselstats;

  void resetCounters()
  {
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      evselstats[ii] = 0;
  }

  void fillHistos()
  {
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      histos.fill(HIST("hEventSelection"), ii, evselstats[ii]);
  }

  Filter preFilterMcCollisions = aod::mccollisionprop::numRecoCollision > 0;

  Filter preFilterCascade =
    nabs(aod::cascdata::dcapostopv) > v0setting_dcapostopv&& nabs(aod::cascdata::dcanegtopv) > v0setting_dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < v0setting_dcav0dau&& aod::cascdata::dcacascdaughters<cascadesetting_dcacascdau && aod::mccasclabel::mcParticleId> - 1;

  Filter preFilterV0 =
    aod::mcv0label::mcParticleId > -1 && nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv&& nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv&& aod::v0data::dcaV0daughters < v0setting_dcav0dau;

  void init(InitContext const&)
  {
    const AxisSpec axisEventSelection{(int)10, -0.5f, +9.5f, ""};
    histos.add("hEventSelection", "hEventSelection", kTH1F, {axisEventSelection});

    // Creation of axes
    const AxisSpec axisVsPt{(int)NBinsPt, 0, MaxPt, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisVsPtCoarse{(int)NBinsPtCoarse, 0, MaxPt, "#it{p}_{T} (GeV/c)"};

    const AxisSpec axisGammaMass{400, 0.000f, 0.400f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisK0ShortMass{400, 0.400f, 0.600f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisLambdaMass{400, 1.01f, 1.21f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisXiMass{400, 1.22f, 1.42f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisOmegaMass{400, 1.57f, 1.77f, "Inv. Mass (GeV/c^{2})"};

    const AxisSpec axisV0Radius{200, 0.0f, 50.0f, "V0 decay radius (cm)"};
    const AxisSpec axisCascRadius{200, 0.0f, 50.0f, "Cascade decay radius (cm)"};
    const AxisSpec axisDCA{200, -2.0f, +2.0f, "DCA single-track to PV (cm)"};
    const AxisSpec axisDCADaughters{200, 0.0f, +2.0f, "DCA between daughters (cm)"};
    const AxisSpec axisDCAWD{200, 0.0f, +2.0f, "DCA to PV (cm)"};
    const AxisSpec axisPA{200, 0.0f, +1.0f, "Pointing angle (rad)"};

    const AxisSpec axisITSClu{10, -0.5f, +9.5f, "ITS clusters"};
    const AxisSpec axisTPCCroRo{160, -0.5f, +159.5f, "TPC crossed rows"};

    // bit packed ITS cluster map
    const AxisSpec axisITSCluMap{(int)128, -0.5f, +127.5f, "Packed ITS map"};
    const AxisSpec axisRadius{(int)160, 0.0f, +80.0f, "Radius (cm)"};

    TString lSpecies[] = {"Gamma", "K0Short", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
    const AxisSpec lMassAxis[] = {axisGammaMass, axisK0ShortMass, axisLambdaMass, axisLambdaMass, axisXiMass, axisXiMass, axisOmegaMass, axisOmegaMass};

    // Creation of histograms: MC generated
    for (Int_t i = 0; i < 8; i++)
      histos.add(Form("hGen%s", lSpecies[i].Data()), Form("hGen%s", lSpecies[i].Data()), kTH1F, {axisVsPt});
    for (Int_t i = 0; i < 8; i++)
      histos.add(Form("hGenWithPV%s", lSpecies[i].Data()), Form("hGenWithPV%s", lSpecies[i].Data()), kTH1F, {axisVsPt});

    // Creation of histograms: mass affairs
    for (Int_t i = 0; i < 8; i++)
      histos.add(Form("h2dMass%s", lSpecies[i].Data()), Form("h2dMass%s", lSpecies[i].Data()), kTH2F, {axisVsPt, lMassAxis[i]});

    // Topo sel QA
    histos.add("h2dK0ShortQAV0Radius", "h2dK0ShortQAV0Radius", kTH2F, {axisVsPtCoarse, axisV0Radius});
    histos.add("h2dK0ShortQADCAV0Dau", "h2dK0ShortQADCAV0Dau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dK0ShortQADCAPosToPV", "h2dK0ShortQADCAPosToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dK0ShortQADCANegToPV", "h2dK0ShortQADCANegToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dK0ShortQADCAToPV", "h2dK0ShortQADCAToPV", kTH2F, {axisVsPtCoarse, axisDCAWD});
    histos.add("h2dK0ShortQAPointingAngle", "h2dK0ShortQAPointingAngle", kTH2F, {axisVsPtCoarse, axisPA});

    histos.add("h2dLambdaQAV0Radius", "h2dLambdaQAV0Radius", kTH2F, {axisVsPtCoarse, axisV0Radius});
    histos.add("h2dLambdaQADCAV0Dau", "h2dLambdaQADCAV0Dau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dLambdaQADCAPosToPV", "h2dLambdaQADCAPosToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dLambdaQADCANegToPV", "h2dLambdaQADCANegToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dLambdaQADCAToPV", "h2dLambdaQADCAToPV", kTH2F, {axisVsPtCoarse, axisDCAWD});
    histos.add("h2dLambdaQAPointingAngle", "h2dLambdaQAPointingAngle", kTH2F, {axisVsPtCoarse, axisPA});

    histos.add("h2dXiMinusQAV0Radius", "h2dXiMinusQAV0Radius", kTH2F, {axisVsPtCoarse, axisV0Radius});
    histos.add("h2dXiMinusQACascadeRadius", "h2dXiMinusQACascadeRadius", kTH2F, {axisVsPtCoarse, axisCascRadius});
    histos.add("h2dXiMinusQADCAV0Dau", "h2dXiMinusQADCAV0Dau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dXiMinusQADCACascDau", "h2dXiMinusQADCACascDau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dXiMinusQADCAPosToPV", "h2dXiMinusQADCAPosToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dXiMinusQADCANegToPV", "h2dXiMinusQADCANegToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dXiMinusQADCABachToPV", "h2dXiMinusQADCABachToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dXiMinusQADCACascToPV", "h2dXiMinusQADCACascToPV", kTH2F, {axisVsPtCoarse, axisDCAWD});
    histos.add("h2dXiMinusQAPointingAngle", "h2dXiMinusQAPointingAngle", kTH2F, {axisVsPtCoarse, axisPA});

    histos.add("h2dOmegaMinusQAV0Radius", "h2dOmegaMinusQAV0Radius", kTH2F, {axisVsPtCoarse, axisV0Radius});
    histos.add("h2dOmegaMinusQACascadeRadius", "h2dOmegaMinusQACascadeRadius", kTH2F, {axisVsPtCoarse, axisCascRadius});
    histos.add("h2dOmegaMinusQADCAV0Dau", "h2dOmegaMinusQADCAV0Dau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dOmegaMinusQADCACascDau", "h2dOmegaMinusQADCACascDau", kTH2F, {axisVsPtCoarse, axisDCADaughters});
    histos.add("h2dOmegaMinusQADCAPosToPV", "h2dOmegaMinusQADCAPosToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dOmegaMinusQADCANegToPV", "h2dOmegaMinusQADCANegToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dOmegaMinusQADCABachToPV", "h2dOmegaMinusQADCABachToPV", kTH2F, {axisVsPtCoarse, axisDCA});
    histos.add("h2dOmegaMinusQADCACascToPV", "h2dOmegaMinusQADCACascToPV", kTH2F, {axisVsPtCoarse, axisDCAWD});
    histos.add("h2dOmegaMinusQAPointingAngle", "h2dOmegaMinusQAPointingAngle", kTH2F, {axisVsPtCoarse, axisPA});

    histos.add("h2dITSCluMap_V0Positive", "h2dITSCluMap_V0Positive", kTH2D, {axisITSCluMap, axisRadius});
    histos.add("h2dITSCluMap_V0Negative", "h2dITSCluMap_V0Negative", kTH2D, {axisITSCluMap, axisRadius});
    histos.add("h2dITSCluMap_CascPositive", "h2dITSCluMap_CascPositive", kTH2D, {axisITSCluMap, axisRadius});
    histos.add("h2dITSCluMap_CascNegative", "h2dITSCluMap_CascNegative", kTH2D, {axisITSCluMap, axisRadius});
    histos.add("h2dITSCluMap_CascBachelor", "h2dITSCluMap_CascBachelor", kTH2D, {axisITSCluMap, axisRadius});

    // Track quality tests
    histos.add("h3dTrackPtsK0ShortP", "h3dTrackPtsK0ShortP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsK0ShortN", "h3dTrackPtsK0ShortN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsLambdaP", "h3dTrackPtsLambdaP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsLambdaN", "h3dTrackPtsLambdaN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsAntiLambdaP", "h3dTrackPtsAntiLambdaP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsAntiLambdaN", "h3dTrackPtsAntiLambdaN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiMinusP", "h3dTrackPtsXiMinusP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiMinusN", "h3dTrackPtsXiMinusN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiMinusB", "h3dTrackPtsXiMinusB", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiPlusP", "h3dTrackPtsXiPlusP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiPlusN", "h3dTrackPtsXiPlusN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsXiPlusB", "h3dTrackPtsXiPlusB", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaMinusP", "h3dTrackPtsOmegaMinusP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaMinusN", "h3dTrackPtsOmegaMinusN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaMinusB", "h3dTrackPtsOmegaMinusB", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaPlusP", "h3dTrackPtsOmegaPlusP", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaPlusN", "h3dTrackPtsOmegaPlusN", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});
    histos.add("h3dTrackPtsOmegaPlusB", "h3dTrackPtsOmegaPlusB", kTH3F, {axisVsPtCoarse, axisITSClu, axisTPCCroRo});

    resetCounters();

    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "Sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");

    histos.add("hDCAxyVsPt", "hDCAxyVsPt", kTH2F, {{20, 0.0f, 10.0f}, {10000, -0.5f, 0.5f}});
    histos.add("hDCAzVsPt", "hDCAzVsPt", kTH2F, {{20, 0.0f, 10.0f}, {10000, -0.5f, 0.5f}});
  }

  void processV0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<V0DataLabeled> const& fullV0s, soa::Filtered<CascMC> const&, TracksCompleteIUMC const& tracks, aod::McParticles const&, aod::V0sLinked const&)
  {
    evselstats[kEvSelAll]++;
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    evselstats[kEvSelBool]++;
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    evselstats[kEvSelVtxZ]++;
    for (auto const& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (mcParticle.isPhysicalPrimary() && track.itsNCls() >= 7)
          histos.fill(HIST("hDCAxyVsPt"), track.pt(), track.dcaXY());
        histos.fill(HIST("hDCAzVsPt"), track.pt(), track.dcaXY());
      }
    }

    for (auto& v0 : fullV0s) {
      // MC association
      auto posPartTrack = v0.posTrack_as<TracksCompleteIUMC>();
      auto negPartTrack = v0.negTrack_as<TracksCompleteIUMC>();
      if (!v0.has_mcParticle() || !posPartTrack.has_mcParticle() || !negPartTrack.has_mcParticle())
        continue;
      auto v0mc = v0.mcParticle();
      if (TMath::Abs(v0mc.y()) > 0.5)
        continue;

      // fill track quality
      if (v0mc.pdgCode() == 310) {
        histos.fill(HIST("h3dTrackPtsK0ShortP"), v0.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsK0ShortN"), v0.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      }
      if (v0mc.pdgCode() == 3122) {
        histos.fill(HIST("h3dTrackPtsLambdaP"), v0.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsLambdaN"), v0.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      }
      if (v0mc.pdgCode() == -3122) {
        histos.fill(HIST("h3dTrackPtsAntiLambdaP"), v0.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsAntiLambdaN"), v0.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      }

      if (posPartTrack.itsNCls() < itsminclusters || negPartTrack.itsNCls() < itsminclusters)
        continue;
      if (posPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || negPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows)
        continue;

      if (v0mc.pdgCode() == 310) {
        histos.fill(HIST("h2dK0ShortQAV0Radius"), v0.pt(), v0.v0radius());
        histos.fill(HIST("h2dK0ShortQADCAV0Dau"), v0.pt(), v0.dcaV0daughters());
        histos.fill(HIST("h2dK0ShortQADCAPosToPV"), v0.pt(), v0.posTrack_as<TracksCompleteIUMC>().dcaXY());
        histos.fill(HIST("h2dK0ShortQADCANegToPV"), v0.pt(), v0.negTrack_as<TracksCompleteIUMC>().dcaXY());
        histos.fill(HIST("h2dK0ShortQADCAToPV"), v0.pt(), v0.dcav0topv());
        histos.fill(HIST("h2dK0ShortQAPointingAngle"), v0.pt(), TMath::ACos(v0.v0cosPA()));
      }
      if (v0mc.pdgCode() == 3122) {
        histos.fill(HIST("h2dLambdaQAV0Radius"), v0.pt(), v0.v0radius());
        histos.fill(HIST("h2dLambdaQADCAV0Dau"), v0.pt(), v0.dcaV0daughters());
        histos.fill(HIST("h2dLambdaQADCAPosToPV"), v0.pt(), v0.posTrack_as<TracksCompleteIUMC>().dcaXY());
        histos.fill(HIST("h2dLambdaQADCANegToPV"), v0.pt(), v0.negTrack_as<TracksCompleteIUMC>().dcaXY());
        histos.fill(HIST("h2dLambdaQADCAToPV"), v0.pt(), v0.dcav0topv());
        histos.fill(HIST("h2dLambdaQAPointingAngle"), v0.pt(), TMath::ACos(v0.v0cosPA()));
      }

      if (v0.v0radius() > v0setting_radius && v0.v0radius() < maxV0Radius) {
        if (v0.v0cosPA() > v0setting_cospa) {
          if (v0.dcaV0daughters() < v0setting_dcav0dau) {
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // Fill invariant masses
            if (v0mc.pdgCode() == 310)
              histos.fill(HIST("h2dMassK0Short"), v0.pt(), v0.mK0Short());
            if (v0mc.pdgCode() == 3122)
              histos.fill(HIST("h2dMassLambda"), v0.pt(), v0.mLambda());
            if (v0mc.pdgCode() == -3122)
              histos.fill(HIST("h2dMassAntiLambda"), v0.pt(), v0.mAntiLambda());
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
          }
        }
      }
    } // end v0 loop
    fillHistos();
    resetCounters();
  }
  PROCESS_SWITCH(straRecoStudy, processV0, "Regular V0 analysis", true);

  void processV0RealData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& v0tables, aod::CascDataExt const&, TracksCompleteIU const&, aod::V0sLinked const&)
  {
    evselstats[kEvSelAll]++;
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    evselstats[kEvSelBool]++;
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    evselstats[kEvSelVtxZ]++;
    for (auto& v0 : v0tables) {
      // MC association
      auto posPartTrack = v0.posTrack_as<TracksCompleteIU>();
      auto negPartTrack = v0.negTrack_as<TracksCompleteIU>();
      histos.fill(HIST("h2dITSCluMap_V0Positive"), (float)posPartTrack.itsClusterMap(), v0.v0radius());
      histos.fill(HIST("h2dITSCluMap_V0Negative"), (float)negPartTrack.itsClusterMap(), v0.v0radius());

      if (posPartTrack.itsNCls() < itsminclusters || negPartTrack.itsNCls() < itsminclusters)
        continue;
      if (posPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || negPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows)
        continue;

      if (v0.v0radius() > v0setting_radius && v0.v0radius() < maxV0Radius) {
        if (v0.v0cosPA() > v0setting_cospa) {
          if (v0.dcaV0daughters() < v0setting_dcav0dau) {
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // Fill invariant masses
            histos.fill(HIST("h2dMassK0Short"), v0.pt(), v0.mK0Short());
            histos.fill(HIST("h2dMassLambda"), v0.pt(), v0.mLambda());
            histos.fill(HIST("h2dMassAntiLambda"), v0.pt(), v0.mAntiLambda());
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
          }
        }
      }
    } // end v0 loop
    fillHistos();
    resetCounters();
  }
  PROCESS_SWITCH(straRecoStudy, processV0RealData, "Regular V0 analysis in real data", false);

  template <typename TCascCandidate>
  void processCascadeCandidate(TCascCandidate casc)
  {
    // MC association
    if (!casc.has_collision()) {
      return;
    }
    auto collision = casc.template collision_as<CollisionsWithEvSels>();
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }

    if (!casc.has_mcParticle())
      return;
    auto cascmc = casc.mcParticle();
    if (TMath::Abs(cascmc.y()) > 0.5)
      return;

    auto bachPartTrack = casc.template bachelor_as<TracksCompleteIUMC>();
    auto posPartTrack = casc.template posTrack_as<TracksCompleteIUMC>();
    auto negPartTrack = casc.template negTrack_as<TracksCompleteIUMC>();

    if (cascmc.pdgCode() == 3312) {
      histos.fill(HIST("h3dTrackPtsXiMinusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsXiMinusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsXiMinusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
    }
    if (cascmc.pdgCode() == -3312) {
      histos.fill(HIST("h3dTrackPtsXiPlusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsXiPlusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsXiPlusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
    }
    if (cascmc.pdgCode() == 3334) {
      histos.fill(HIST("h3dTrackPtsOmegaMinusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsOmegaMinusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsOmegaMinusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
    }
    if (cascmc.pdgCode() == -3334) {
      histos.fill(HIST("h3dTrackPtsOmegaPlusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsOmegaPlusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
      histos.fill(HIST("h3dTrackPtsOmegaPlusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
    }

    if (posPartTrack.itsNCls() < itsminclusters || negPartTrack.itsNCls() < itsminclusters || bachPartTrack.itsNCls() < itsminclusters)
      return;
    if (posPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || negPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || bachPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows)
      return;

    if (cascmc.pdgCode() == 3312) {
      histos.fill(HIST("h2dXiMinusQAV0Radius"), casc.pt(), casc.v0radius());
      histos.fill(HIST("h2dXiMinusQACascadeRadius"), casc.pt(), casc.cascradius());
      histos.fill(HIST("h2dXiMinusQADCAV0Dau"), casc.pt(), casc.dcaV0daughters());
      histos.fill(HIST("h2dXiMinusQADCACascDau"), casc.pt(), casc.dcacascdaughters());
      histos.fill(HIST("h2dXiMinusQADCAPosToPV"), casc.pt(), casc.dcapostopv());
      histos.fill(HIST("h2dXiMinusQADCANegToPV"), casc.pt(), casc.dcanegtopv());
      histos.fill(HIST("h2dXiMinusQADCABachToPV"), casc.pt(), casc.dcabachtopv());
      histos.fill(HIST("h2dXiMinusQADCACascToPV"), casc.pt(), casc.dcaXYCascToPV());
      histos.fill(HIST("h2dXiMinusQAPointingAngle"), casc.pt(), TMath::ACos(casc.casccosPA(collision.posX(), collision.posY(), collision.posZ())));
    }
    if (cascmc.pdgCode() == 3334) {
      histos.fill(HIST("h2dOmegaMinusQAV0Radius"), casc.pt(), casc.v0radius());
      histos.fill(HIST("h2dOmegaMinusQACascadeRadius"), casc.pt(), casc.cascradius());
      histos.fill(HIST("h2dOmegaMinusQADCAV0Dau"), casc.pt(), casc.dcaV0daughters());
      histos.fill(HIST("h2dOmegaMinusQADCACascDau"), casc.pt(), casc.dcacascdaughters());
      histos.fill(HIST("h2dOmegaMinusQADCAPosToPV"), casc.pt(), casc.dcapostopv());
      histos.fill(HIST("h2dOmegaMinusQADCANegToPV"), casc.pt(), casc.dcanegtopv());
      histos.fill(HIST("h2dOmegaMinusQADCABachToPV"), casc.pt(), casc.dcabachtopv());
      histos.fill(HIST("h2dOmegaMinusQADCACascToPV"), casc.pt(), casc.dcaXYCascToPV());
      histos.fill(HIST("h2dOmegaMinusQAPointingAngle"), casc.pt(), TMath::ACos(casc.casccosPA(collision.posX(), collision.posY(), collision.posZ())));
    }

    if (casc.v0radius() < maxV0Radius && casc.cascradius() < maxCascRadius) {
      if (casc.v0radius() > v0setting_radius && casc.cascradius() > cascadesetting_cascradius) {
        if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0setting_cospa) {
          if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > cascadesetting_cospa) {
            if (casc.dcaV0daughters() < v0setting_dcav0dau) {
              //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
              // Fill invariant masses
              if (cascmc.pdgCode() == 3312)
                histos.fill(HIST("h2dMassXiMinus"), casc.pt(), casc.mXi());
              if (cascmc.pdgCode() == -3312)
                histos.fill(HIST("h2dMassXiPlus"), casc.pt(), casc.mXi());
              if (cascmc.pdgCode() == 3334)
                histos.fill(HIST("h2dMassOmegaMinus"), casc.pt(), casc.mOmega());
              if (cascmc.pdgCode() == -3334)
                histos.fill(HIST("h2dMassOmegaPlus"), casc.pt(), casc.mOmega());
              //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            }
          }
        }
      }
    }
  }

  void processCascade(CollisionsWithEvSels const&, aod::V0Datas const&, soa::Filtered<CascMC> const& Cascades, TracksCompleteIUMC const&, aod::McParticles const&)
  {
    for (auto& casc : Cascades) {
      processCascadeCandidate(casc);
    }
  }
  PROCESS_SWITCH(straRecoStudy, processCascade, "Regular cascade analysis", true);
  void processTrackedCascade(CollisionsWithEvSels const&, aod::V0Datas const&, TraCascMC const& Cascades, TracksCompleteIUMC const&, aod::McParticles const&)
  {
    for (auto& casc : Cascades) {
      processCascadeCandidate(casc);
    }
  }
  PROCESS_SWITCH(straRecoStudy, processTrackedCascade, "Tracked cascade analysis", true);

  void processCascadeRealData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const&, aod::CascDataExt const& Cascades, TracksCompleteIU const&, aod::V0sLinked const&)
  {
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    for (auto& casc : Cascades) {
      auto bachPartTrack = casc.bachelor_as<TracksCompleteIU>();
      auto posPartTrack = casc.posTrack_as<TracksCompleteIU>();
      auto negPartTrack = casc.negTrack_as<TracksCompleteIU>();

      histos.fill(HIST("h2dITSCluMap_CascPositive"), (float)posPartTrack.itsClusterMap(), casc.v0radius());
      histos.fill(HIST("h2dITSCluMap_CascNegative"), (float)negPartTrack.itsClusterMap(), casc.v0radius());
      histos.fill(HIST("h2dITSCluMap_CascBachelor"), (float)bachPartTrack.itsClusterMap(), casc.cascradius());

      if (casc.sign() < 0) {
        histos.fill(HIST("h3dTrackPtsXiMinusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsXiMinusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsXiMinusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
      }
      if (casc.sign() > 0) {
        histos.fill(HIST("h3dTrackPtsXiPlusP"), casc.pt(), posPartTrack.itsNCls(), posPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsXiPlusN"), casc.pt(), negPartTrack.itsNCls(), negPartTrack.tpcNClsCrossedRows());
        histos.fill(HIST("h3dTrackPtsXiPlusB"), casc.pt(), bachPartTrack.itsNCls(), bachPartTrack.tpcNClsCrossedRows());
      }

      if (posPartTrack.itsNCls() < itsminclusters || negPartTrack.itsNCls() < itsminclusters || bachPartTrack.itsNCls() < itsminclusters)
        continue;
      if (posPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || negPartTrack.tpcNClsCrossedRows() < tpcmincrossedrows || bachPartTrack.itsNCls() < tpcmincrossedrows)
        continue;

      if (casc.v0radius() < maxV0Radius && casc.cascradius() < maxCascRadius) {
        if (casc.v0radius() > v0setting_radius && casc.cascradius() > cascadesetting_cascradius) {
          if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0setting_cospa) {
            if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > cascadesetting_cospa) {
              if (casc.dcaV0daughters() < v0setting_dcav0dau) {
                //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
                // Fill invariant masses
                if (casc.sign() < 0)
                  histos.fill(HIST("h2dMassXiMinus"), casc.pt(), casc.mXi());
                if (casc.sign() > 0)
                  histos.fill(HIST("h2dMassXiPlus"), casc.pt(), casc.mXi());
                if (casc.sign() < 0)
                  histos.fill(HIST("h2dMassOmegaMinus"), casc.pt(), casc.mOmega());
                if (casc.sign() > 0)
                  histos.fill(HIST("h2dMassOmegaPlus"), casc.pt(), casc.mOmega());
                //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(straRecoStudy, processCascadeRealData, "Regular cascade analysis, real data", false);

  void processGeneratedReconstructible(soa::Filtered<RecoedMCCollisions>::iterator const&, aod::McParticles const& mcParticles)
  {
    // check if collision successfully reconstructed
    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.eta()) < 0.5) {
        if (mcp.pdgCode() == 310)
          histos.fill(HIST("hGenWithPVK0Short"), mcp.pt());
        if (mcp.pdgCode() == 3122)
          histos.fill(HIST("hGenWithPVLambda"), mcp.pt());
        if (mcp.pdgCode() == -3122)
          histos.fill(HIST("hGenWithPVAntiLambda"), mcp.pt());
        if (mcp.pdgCode() == 3312)
          histos.fill(HIST("hGenWithPVXiMinus"), mcp.pt());
        if (mcp.pdgCode() == -3312)
          histos.fill(HIST("hGenWithPVXiPlus"), mcp.pt());
        if (mcp.pdgCode() == 3334)
          histos.fill(HIST("hGenWithPVOmegaMinus"), mcp.pt());
        if (mcp.pdgCode() == -3334)
          histos.fill(HIST("hGenWithPVOmegaPlus"), mcp.pt());
      }
    }
  }
  PROCESS_SWITCH(straRecoStudy, processGeneratedReconstructible, "generated analysis in events with PV", true);

  void processPureGenerated(aod::McParticles const& mcParticles)
  {
    // check if collision successfully reconstructed
    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.eta()) < 0.5) {
        if (mcp.pdgCode() == 22)
          histos.fill(HIST("hGenGamma"), mcp.pt());
        if (mcp.pdgCode() == 310)
          histos.fill(HIST("hGenK0Short"), mcp.pt());
        if (mcp.pdgCode() == 3122)
          histos.fill(HIST("hGenLambda"), mcp.pt());
        if (mcp.pdgCode() == -3122)
          histos.fill(HIST("hGenAntiLambda"), mcp.pt());
        if (mcp.pdgCode() == 3312)
          histos.fill(HIST("hGenXiMinus"), mcp.pt());
        if (mcp.pdgCode() == -3312)
          histos.fill(HIST("hGenXiPlus"), mcp.pt());
        if (mcp.pdgCode() == 3334)
          histos.fill(HIST("hGenOmegaMinus"), mcp.pt());
        if (mcp.pdgCode() == -3334)
          histos.fill(HIST("hGenOmegaPlus"), mcp.pt());
      }
    }
  }
  PROCESS_SWITCH(straRecoStudy, processPureGenerated, "generated analysis in events with PV", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<preProcessMCcollisions>(cfgc),
    adaptAnalysisTask<straRecoStudy>(cfgc)};
}
