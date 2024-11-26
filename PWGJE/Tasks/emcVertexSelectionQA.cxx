// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Monitoring task for EMCAL vertex selection (number and source of contributors to the vertices)
//
/// \author Nicolas Strangmann <nicolas.strangmann@cern.ch>, Goethe University Frankfurt / Oak Ridge National Laoratory

#include <unordered_map>
#include <vector>
#include <TRobustEstimator.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using bcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using collEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA>;

struct EmcVertexSelectionQA {
  o2::framework::HistogramRegistry mHistManager{"EMCALVertexSelectionQAHistograms"};

  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 20.f, "max. Zvtx"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
  Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", false, "require FT0AND in event cut"};
  Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
  Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
  Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
  Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};

  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;

    o2Axis matchingAxis{4, -0.5, 3.5, "matchingStatus", "Matching status"}; // 0, no vertex,1 vertex found , 2 vertices found, 3 or more vertices found
    o2Axis qualityAxis{3, -0.5, 2.5, "VertexQuality", "Vertex quality"};    // 0: TPC/ITS contributor, 1: TRD contributor , 2: TOF contributor
    o2Axis DCAStdDevAxis{500, 0., 0.05, "StdDevDCA", "#sigma DCA"};
    o2Axis TimeStdDevAxis{1000, 0., 10000., "TimeStdDev", "#sigma t"};

    mHistManager.add("hCollisionMatching", "Collision Status", o2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingReadout", "Collision Status EMCAL Readout", o2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hnVerticeswithTOForTRDcontr", "Number of vertices vs number with TOF/TRD Contributors", o2HistType::kTH2F, {matchingAxis, matchingAxis});
    mHistManager.add("hnVertexContributors", "Number of Contributors to first and second vertex", o2HistType::kTH2F, {{50, 0, 50}, {50, 0, 50}});

    // The standard deviations of the DCA and times of tracks contributing to the vertex are investigated as proxies for the quality of a vertex

    // First using the simple RMS from TMath
    mHistManager.add("hVertexStdDevDCA", "StdDev DCA of Vertex vs its Quality", o2HistType::kTH2F, {DCAStdDevAxis, qualityAxis});
    mHistManager.add("hVertexTrackTimeStdDev", "Standard Deviation of Tracks Contributing to the Vertex vs its Quality", o2HistType::kTH2F, {TimeStdDevAxis, qualityAxis});

    // And secondly using the TRobustEstimator to exclude outliers
    mHistManager.add("hVertexRobustStdDevDCA", "Robust StdDev DCA of Vertex vs its Quality", o2HistType::kTH2F, {DCAStdDevAxis, qualityAxis});
    mHistManager.add("hVertexTrackTimeRobustStdDev", "Robust Standard Deviation of Tracks Contributing to the Vertex vs its Quality", o2HistType::kTH2F, {TimeStdDevAxis, qualityAxis});

    // The ratio of normal/robust RMS to see how much the outlier rejection has impacted the RMS
    mHistManager.add("hVertexRelDiffRobustStdDevDCA", "Relative Difference of Robust StdDev DCA to StdDev DCA of Vertex vs its Quality", o2HistType::kTH2F, {{200, 0, 2}, qualityAxis});
    mHistManager.add("hVertexRelDiffRobustStdDevTrackTime", "Relative Difference of Robust Standard Deviation to Standard Deviation of Tracks Contributing to the Vertex vs its Quality", o2HistType::kTH2F, {{200, 0, 2}, qualityAxis});

    // Z vertex positions of two vertices in the same BC (to look for two collisions in the same BC with similar z-vertex positions)
    mHistManager.add("hDoubleZVertex", "Z Vertex Position Correlation of two Collions in the same BC", o2HistType::kTH2F, {{300, -15., 15.}, {300, -15., 15.}});
    mHistManager.add("hZVertexDiff", "Difference Z Vertex Position of two Collions in the same BC", o2HistType::kTH1F, {{40000, -20., 20., "Z_{vtx 1}-Z_{vtx 2} (cm)"}});

    // Set axis labels and bin titles
    initVertexHistogram(mHistManager.get<TH1>(HIST("hCollisionMatching")).get());
    initVertexHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingReadout")).get());
    initVertexHistogram(mHistManager.get<TH2>(HIST("hnVerticeswithTOForTRDcontr")).get());

    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexStdDevDCA")).get());
    mHistManager.get<TH2>(HIST("hVertexStdDevDCA")).get()->GetXaxis()->SetTitle("<DCA> of vtx contributers");
    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexTrackTimeStdDev")).get());
    mHistManager.get<TH2>(HIST("hVertexTrackTimeStdDev")).get()->GetXaxis()->SetTitle("#sigma t of vtx contributers");
    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexRobustStdDevDCA")).get());
    mHistManager.get<TH2>(HIST("hVertexRobustStdDevDCA")).get()->GetXaxis()->SetTitle("Robust <DCA> of vtx contributers");
    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexTrackTimeRobustStdDev")).get());
    mHistManager.get<TH2>(HIST("hVertexTrackTimeRobustStdDev")).get()->GetXaxis()->SetTitle("Robust #sigma t of vtx contributers");
    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexRelDiffRobustStdDevDCA")).get());
    initVertexQualityAxis(mHistManager.get<TH2>(HIST("hVertexRelDiffRobustStdDevTrackTime")).get());

    mHistManager.get<TH2>(HIST("hnVertexContributors")).get()->GetXaxis()->SetTitle("N_{contr} to Vtx 1");
    mHistManager.get<TH2>(HIST("hnVertexContributors")).get()->GetYaxis()->SetTitle("N_{contr} to Vtx 2");
    mHistManager.get<TH2>(HIST("hDoubleZVertex")).get()->GetXaxis()->SetTitle("Z_{vtx 1} (cm)");
    mHistManager.get<TH2>(HIST("hDoubleZVertex")).get()->GetYaxis()->SetTitle("Z_{vtx 2} (cm)");
  }

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;
  PresliceUnsorted<collEventSels> perFoundBC = aod::evsel::foundBCId;
  void process(bcEvSels const& bcs, collEventSels const& collisions, FullTracksIU const& tracks)
  {
    for (const auto& bc : bcs) {
      bool isEMCALreadout = (bc.alias_bit(kTVXinEMC) || bc.alias_bit(kEMC7) || bc.alias_bit(kEG1) || bc.alias_bit(kEG2) || bc.alias_bit(kDG1) || bc.alias_bit(kDG2) || bc.alias_bit(kEJ1) || bc.alias_bit(kEJ2) || bc.alias_bit(kDJ1) || bc.alias_bit(kDJ2));

      auto colsinbc = collisions.sliceBy(perFoundBC, bc.globalIndex());

      bool isBCAccepted = true;
      for (auto& col : colsinbc) {
        if (cfgRequireSel8 && !col.sel8())
          isBCAccepted = false;
        if (cfgRequireFT0AND && !col.selection_bit(o2::aod::evsel::kIsTriggerTVX))
          isBCAccepted = false;
        if (col.posZ() < -cfgZvtxMax || col.posZ() > cfgZvtxMax)
          isBCAccepted = false;
        if (cfgRequireNoTFB && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
          isBCAccepted = false;
        if (cfgRequireNoITSROFB && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
          isBCAccepted = false;
        if (cfgRequireNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
          isBCAccepted = false;
        if (cfgRequireVertexITSTPC && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
          isBCAccepted = false;
        if (cfgRequireGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
          isBCAccepted = false;
      }
      if (!isBCAccepted)
        continue;

      int collisionStatus = -1;
      if (!colsinbc.size()) {
        collisionStatus = 0;
      } else if (colsinbc.size() == 1) {
        collisionStatus = 1;
      } else if (colsinbc.size() == 2) {
        collisionStatus = 2;
      } else {
        collisionStatus = 3;
      }
      if (collisionStatus >= 0) {
        mHistManager.fill(HIST("hCollisionMatching"), collisionStatus);
        if (isEMCALreadout) {
          mHistManager.fill(HIST("hCollisionMatchingReadout"), collisionStatus);
        }
      }
      if (collisionStatus > 0) {
        int nVtx = 0;
        int nVtxwithTOForTRDcontr = 0;
        std::vector<float> zVertexPositions;
        std::vector<int> nVtxContributors;
        for (auto& col : colsinbc) { // Loop over all collisions/vertices
          int ivtxquality = 0;       // 0: TPC/ITS contributor, 1: TRD contributor , 2: TOF contributor
          // int nITStracks = 0;
          // int nTPCtracks = 0;
          int nTOFtracks = 0;
          int nTRDtracks = 0;
          int nPVContributorTracks = 0;
          std::vector<double> TrackDCA;
          std::vector<double> TrackTime;
          auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
          for (auto& track : tracksGrouped) {
            if (!track.isPVContributor()) {
              continue;
            }
            // nITStracks += track.hasITS();
            // nTPCtracks += track.hasTPC();
            nTOFtracks += track.hasTOF();
            nTRDtracks += track.hasTRD() && !track.hasTOF();
            TrackDCA.push_back(track.dcaXY());
            TrackTime.push_back(track.trackTime());
            // LOGF(info, "TOF: %d, TRD: %d - track.dcaXY[%d] = %f, track.trackTime[%d] = %f, pT = %f, nTPCrows = %d", track.hasTOF(), track.hasTRD(), nPVContributorTracks, track.dcaXY(), nPVContributorTracks, track.trackTime(), track.pt(), track.tpcNClsCrossedRows());

            nPVContributorTracks++;
          }

          if (nTRDtracks > 0) {
            ivtxquality = 1;
          }
          if (nTOFtracks > 0) {
            ivtxquality = 2;
          }

          // Calculate the arithmetic RMS for all tracks contributing to the vertex
          double StdDevTrackDCA = TMath::RMS(TrackDCA.begin(), TrackDCA.end());
          double StdDevTrackTime = TMath::RMS(TrackTime.begin(), TrackTime.end());

          // Calculate the RMS using the robust estimator to reject outliers
          TRobustEstimator robustEstimator;
          double RobustMeanTrackDCA = 0, RobustStdDevTrackDCA = 0;
          double RobustMeanTrackTime = 0, RobustStdDevTrackTime = 0;
          robustEstimator.EvaluateUni(TrackDCA.size(), TrackDCA.data(), RobustMeanTrackDCA, RobustStdDevTrackDCA, 0.5 * TrackDCA.size());
          robustEstimator.EvaluateUni(TrackTime.size(), TrackTime.data(), RobustMeanTrackTime, RobustStdDevTrackTime, 0.5 * TrackTime.size());

          mHistManager.fill(HIST("hVertexStdDevDCA"), StdDevTrackDCA, ivtxquality);
          mHistManager.fill(HIST("hVertexTrackTimeStdDev"), StdDevTrackTime, ivtxquality);
          mHistManager.fill(HIST("hVertexRobustStdDevDCA"), RobustStdDevTrackDCA, ivtxquality);
          mHistManager.fill(HIST("hVertexTrackTimeRobustStdDev"), RobustStdDevTrackTime, ivtxquality);
          mHistManager.fill(HIST("hVertexRelDiffRobustStdDevDCA"), RobustStdDevTrackDCA / StdDevTrackDCA, ivtxquality);
          mHistManager.fill(HIST("hVertexRelDiffRobustStdDevTrackTime"), RobustStdDevTrackTime / StdDevTrackTime, ivtxquality);
          // LOGF(info, "StdDevTrackDCA = %f, StdDevTrackTime = %f, RobustStdDevTrackDCA = %f, RobustStdDevTrackTime = %f", StdDevTrackDCA, StdDevTrackTime, RobustStdDevTrackDCA, RobustStdDevTrackTime);

          nVtx++;
          if (nTOFtracks + nTRDtracks > 0) {
            nVtxwithTOForTRDcontr++;
          }
          nVtxContributors.push_back(nPVContributorTracks);
          zVertexPositions.push_back(col.posZ());
        }
        mHistManager.fill(HIST("hnVerticeswithTOForTRDcontr"), nVtx, nVtxwithTOForTRDcontr);
        if (collisionStatus == 2) {
          mHistManager.fill(HIST("hnVertexContributors"), nVtxContributors.at(0), nVtxContributors.at(1));
          mHistManager.fill(HIST("hDoubleZVertex"), zVertexPositions.at(0), zVertexPositions.at(1));
          mHistManager.fill(HIST("hZVertexDiff"), zVertexPositions.at(0) - zVertexPositions.at(1));
        }
        nVtxContributors.clear();
        zVertexPositions.clear();
      }
    }
  }

  // Beautify bc type lables
  void
    initVertexHistogram(TH1* hist)
  {
    hist->GetXaxis()->SetTitle("N_{Vtx}");
    hist->GetXaxis()->SetBinLabel(1, "No vertex");
    hist->GetXaxis()->SetBinLabel(2, "1 vertex");
    hist->GetXaxis()->SetBinLabel(3, "2 vertices");
    hist->GetXaxis()->SetBinLabel(4, "3+ vertices");
    hist->GetYaxis()->SetTitle("Number of BCs");
  }
  // Beautify bc type lables
  void initVertexHistogram(TH2* hist)
  {
    hist->GetXaxis()->SetTitle("N_{Vtx}");
    hist->GetXaxis()->SetBinLabel(1, "No vertex");
    hist->GetXaxis()->SetBinLabel(2, "1 vertex");
    hist->GetXaxis()->SetBinLabel(3, "2 vertices");
    hist->GetXaxis()->SetBinLabel(4, "3+ vertices");
    hist->GetYaxis()->SetTitle("N_{Vtx} with TOF/TRD Contributors");
    hist->GetYaxis()->SetBinLabel(1, "No vertex");
    hist->GetYaxis()->SetBinLabel(2, "1 vertex");
    hist->GetYaxis()->SetBinLabel(3, "2 vertices");
    hist->GetYaxis()->SetBinLabel(4, "3+ vertices");
    hist->GetZaxis()->SetTitle("Number of BCs");
  }
  // Beautify vertex quality lables
  void initVertexQualityAxis(TH2* hist)
  {
    hist->GetYaxis()->SetTitle("Highest contributer quality");
    hist->GetYaxis()->SetBinLabel(1, "TPC/ITS");
    hist->GetYaxis()->SetBinLabel(2, "TRD");
    hist->GetYaxis()->SetBinLabel(3, "TOF");
    hist->GetZaxis()->SetTitle("Number of BCs");
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<EmcVertexSelectionQA>(cfgc)};
}
