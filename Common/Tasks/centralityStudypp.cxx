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
// This task is dedicated to percentile calibration and event
// selection studies in pp collisions using derived data based on the
// multCentTable output

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TProfile.h>

#include <cstdint>
#include <format>
#include <map>
#include <memory>
#include <string>

using namespace o2;
using namespace o2::framework;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
#define getHist(type, name) std::get<std::shared_ptr<type>>(histPointers[name])

struct centralityStudypp {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;
  std::string histPath;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = -1;

  // vertex Z equalization
  TList* hCalibObjects = nullptr;
  TProfile* hVtxZFV0A = nullptr;
  TProfile* hVtxZFT0A = nullptr;
  TProfile* hVtxZFT0C = nullptr;
  TProfile* hVtxZNTracks = nullptr;
  TProfile* hVtxZNGlobals = nullptr;
  TProfile* hVtxZMFT = nullptr;
  TProfile* hVtxZFDDA = nullptr;
  TProfile* hVtxZFDDC = nullptr;

  // calibration histograms
  TH1* hCentralityFV0A = nullptr;
  TH1* hCentralityFT0A = nullptr;
  TH1* hCentralityFT0C = nullptr;
  TH1* hCentralityFT0M = nullptr;
  TH1* hCentralityFDDM = nullptr;
  TH1* hCentralityNTPV = nullptr;
  TH1* hCentralityNGlo = nullptr;
  TH1* hCentralityMFT = nullptr;

  // Configurables
  Configurable<bool> do2DPlots{"do2DPlots", true, "0 - no, 1 - yes"};
  // _______________________________________
  // event selection criteria
  Configurable<bool> applyVertexZEqualization{"applyVertexZEqualization", false, "0 - no, 1 - yes"};
  Configurable<bool> saveUnequalized{"saveUnequalized", false, "save unequalized raw: 0 - no, 1 - yes"};
  Configurable<bool> doCentralityQA{"doCentralityQA", false, "do centrality QA: 0 - no, 1 - yes"};
  Configurable<bool> applySel8{"applySel8", true, "0 - no, 1 - yes"};
  Configurable<bool> applyVtxZ{"applyVtxZ", true, "0 - no, 1 - yes"};
  Configurable<bool> requireINELgtZERO{"requireINELgtZERO", true, "0 no, 1 - yes"};
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", false, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", false, "reject events at TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", false, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions in case of pileup with another collision in the same foundBC"};

  Configurable<float> vertexZcut{"vertexZcut", 10.0f, "vertex Z cut (cm)"};

  // For one-dimensional plots, where binning is no issue
  ConfigurableAxis axisMultUltraFineFV0A{"axisMultUltraFineFV0A", {60000, 0, 60000}, "FV0A amplitude"};
  ConfigurableAxis axisMultUltraFineFT0A{"axisMultUltraFineFT0A", {60000, 0, 60000}, "FT0A amplitude"};
  ConfigurableAxis axisMultUltraFineFT0C{"axisMultUltraFineFT0C", {60000, 0, 60000}, "FT0C amplitude"};
  ConfigurableAxis axisMultUltraFineFT0M{"axisMultUltraFineFT0M", {50000, 0, 200000}, "FT0M amplitude"};
  ConfigurableAxis axisMultUltraFineFDDA{"axisMultUltraFineFDDA", {60000, 0, 60000}, "FDDA amplitude"};
  ConfigurableAxis axisMultUltraFineFDDC{"axisMultUltraFineFDDC", {60000, 0, 60000}, "FDDC amplitude"};
  ConfigurableAxis axisMultUltraFinePVContributors{"axisMultUltraFinePVContributors", {10000, 0, 10000}, "Number of PV Contributors"};
  ConfigurableAxis axisMultUltraFineGlobalTracks{"axisMultUltraFineGlobalTracks", {5000, 0, 5000}, "Number of global tracks"};
  ConfigurableAxis axisMultUltraFineMFTTracks{"axisMultUltraFineMFTTracks", {5000, 0, 5000}, "Number of MFT tracks"};
  // For profile Z
  ConfigurableAxis axisPVz{"axisPVz", {400, -20.0f, +20.0f}, "PVz (cm)"};
  ConfigurableAxis axisZN{"axisZN", {1100, -50.0f, +500.0f}, "ZN"};
  // For centrality QA
  ConfigurableAxis axisCentrality{"axisCentrality", {10000, 0, 100}, "centrality percentile"};

  // ccdb matters
  Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> pathGRPECSObject{"pathGRPECSObject", "GLO/Config/GRPECS", "Path to GRPECS object"};
  Configurable<std::string> pathVertexZ{"pathVertexZ", "Users/d/ddobrigk/Centrality/Calibration", "Path to vertexZ profiles"};
  Configurable<std::string> pathCentrality{"pathCentrality", "Users/d/ddobrigk/Centrality/Estimators", "Path to centrality calibration"};

  void init(InitContext&)
  {
    hCalibObjects = nullptr;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZNTracks = nullptr;
    hVtxZNGlobals = nullptr;
    hVtxZMFT = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;

    hCentralityFV0A = nullptr;
    hCentralityFT0A = nullptr;
    hCentralityFT0C = nullptr;
    hCentralityFT0M = nullptr;
    hCentralityFDDM = nullptr;
    hCentralityNTPV = nullptr;
    hCentralityNGlo = nullptr;
    hCentralityMFT = nullptr;

    const AxisSpec axisCollisions{100, -0.5f, 99.5f, "Number of collisions"};
    histos.add("hCollisionSelection", "hCollisionSelection", kTH1D, {{20, -0.5f, +19.5f}});
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(11, "Neighbour rejection");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(12, "no ITS in-ROF pileup (standard)");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(13, "no ITS in-ROF pileup (strict)");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(14, "is UPC event");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(15, "rejectCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(16, "em/upc rejection");
    histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(17, "isFlangeEvent");

    histos.add("hFT0A_Collisions", "hFT0A_Collisions", kTH1D, {axisMultUltraFineFT0A});
    histos.add("hFT0C_Collisions", "hFT0C_Collisions", kTH1D, {axisMultUltraFineFT0C});
    histos.add("hFT0M_Collisions", "hFT0M_Collisions", kTH1D, {axisMultUltraFineFT0M});

    histos.add("hFDDA_Collisions", "hFDDA_Collisions", kTH1D, {axisMultUltraFineFDDA});
    histos.add("hFDDC_Collisions", "hFDDC_Collisions", kTH1D, {axisMultUltraFineFDDC});

    histos.add("hFV0A_Collisions", "hFV0A_Collisions", kTH1D, {axisMultUltraFineFV0A});
    histos.add("hNGlobalTracks", "hNGlobalTracks", kTH1D, {axisMultUltraFineGlobalTracks});
    histos.add("hNMFTTracks", "hNMFTTracks", kTH1D, {axisMultUltraFineMFTTracks});
    histos.add("hNPVContributors", "hNPVContributors", kTH1D, {axisMultUltraFinePVContributors});

    histos.add("hFT0AvsPVz_Collisions", "hFT0AvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hFT0CvsPVz_Collisions", "hFT0CvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hFDDAvsPVz_Collisions", "hFDDAvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hFDDCvsPVz_Collisions", "hFDDCvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hFV0AvsPVz_Collisions", "hFV0AvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hNGlobalTracksvsPVz_Collisions", "hNGlobalTracksvsPVz_Collisions", kTProfile, {axisPVz});
    histos.add("hNMFTTracksvsPVz_Collisions", "hNMFTTracksvsPVz_Collisions", kTProfile, {axisPVz});
  }

  template <typename TCollision>
  void initRun(const TCollision& collision)
  {
    if (mRunNumber == collision.multRunNumber()) {
      return;
    }

    mRunNumber = collision.multRunNumber();
    LOGF(info, "Setting up for run: %i", mRunNumber);

    if (applyVertexZEqualization.value) {
      // acquire vertex-Z equalization histograms if requested
      LOGF(info, "Acquiring vertex-Z profiles for run %i", mRunNumber);
      hCalibObjects = ccdb->getForRun<TList>(pathVertexZ, mRunNumber);

      hVtxZFV0A = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFV0A"));
      hVtxZFT0A = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFT0A"));
      hVtxZFT0C = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFT0C"));
      hVtxZFDDA = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFDDA"));
      hVtxZFDDC = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFDDC"));
      hVtxZNTracks = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZNTracksPV"));
      hVtxZNGlobals = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZNGlobals"));
      hVtxZMFT = dynamic_cast<TProfile*>(hCalibObjects->FindObject("hVtxZMFT"));

      // Capture error
      if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZFDDA || !hVtxZFDDC || !hVtxZNTracks || !hVtxZNGlobals || !hVtxZMFT) {
        LOGF(error, "Problem loading CCDB objects! Please check");
      }
    }

    if (doCentralityQA.value) {
      // acquire vertex-Z equalization histograms if requested
      LOGF(info, "Acquiring vertex-Z profiles for run %i", mRunNumber);
      TList* hCentralityObjects = nullptr;
      hCentralityObjects = ccdb->getForRun<TList>(pathCentrality, mRunNumber);

      hCentralityFV0A = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqFV0"));
      hCentralityFT0A = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqFT0A"));
      hCentralityFT0C = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqFT0C"));
      hCentralityFT0M = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqFT0"));
      hCentralityFDDM = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqFDD"));
      hCentralityNTPV = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqNTracksPV"));
      hCentralityNGlo = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqNGlobal"));
      hCentralityMFT = dynamic_cast<TProfile*>(hCentralityObjects->FindObject("hCalibZeqMFT"));

      // won't capture null pointers -> explicitly check for those when attempting to evaluate
    }

    histPath = std::format("Run_{}/", mRunNumber);

    histPointers.insert({histPath + "hCollisionSelection", histos.add((histPath + "hCollisionSelection").c_str(), "hCollisionSelection", {kTH1D, {{20, -0.5f, +19.5f}}})});
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(1, "All collisions");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(2, "sel8 cut");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(3, "posZ cut");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
    getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(11, "Pass INEL > 0");

    histPointers.insert({histPath + "hFT0A_Collisions", histos.add((histPath + "hFT0A_Collisions").c_str(), "hFT0A_Collisions", {kTH1D, {{axisMultUltraFineFT0A}}})});
    histPointers.insert({histPath + "hFT0C_Collisions", histos.add((histPath + "hFT0C_Collisions").c_str(), "hFT0C_Collisions", {kTH1D, {{axisMultUltraFineFT0C}}})});
    histPointers.insert({histPath + "hFT0M_Collisions", histos.add((histPath + "hFT0M_Collisions").c_str(), "hFT0M_Collisions", {kTH1D, {{axisMultUltraFineFT0M}}})});

    histPointers.insert({histPath + "hFDDA_Collisions", histos.add((histPath + "hFDDA_Collisions").c_str(), "hFDDA_Collisions", {kTH1D, {{axisMultUltraFineFDDA}}})});
    histPointers.insert({histPath + "hFDDC_Collisions", histos.add((histPath + "hFDDC_Collisions").c_str(), "hFDDC_Collisions", {kTH1D, {{axisMultUltraFineFDDC}}})});

    histPointers.insert({histPath + "hFV0A_Collisions", histos.add((histPath + "hFV0A_Collisions").c_str(), "hFV0A_Collisions", {kTH1D, {{axisMultUltraFineFV0A}}})});
    histPointers.insert({histPath + "hNGlobalTracks", histos.add((histPath + "hNGlobalTracks").c_str(), "hNGlobalTracks", {kTH1D, {{axisMultUltraFineGlobalTracks}}})});
    histPointers.insert({histPath + "hNMFTTracks", histos.add((histPath + "hNMFTTracks").c_str(), "hNMFTTracks", {kTH1D, {{axisMultUltraFineMFTTracks}}})});
    histPointers.insert({histPath + "hNPVContributors", histos.add((histPath + "hNPVContributors").c_str(), "hNPVContributors", {kTH1D, {{axisMultUltraFinePVContributors}}})});

    if (applyVertexZEqualization && saveUnequalized) {
      histPointers.insert({histPath + "hFT0A_Collisions_Unequalized", histos.add((histPath + "hFT0A_Collisions_Unequalized").c_str(), "hFT0A_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFT0A}}})});
      histPointers.insert({histPath + "hFT0C_Collisions_Unequalized", histos.add((histPath + "hFT0C_Collisions_Unequalized").c_str(), "hFT0C_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFT0C}}})});
      histPointers.insert({histPath + "hFT0M_Collisions_Unequalized", histos.add((histPath + "hFT0M_Collisions_Unequalized").c_str(), "hFT0M_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFT0M}}})});

      histPointers.insert({histPath + "hFDDA_Collisions_Unequalized", histos.add((histPath + "hFDDA_Collisions_Unequalized").c_str(), "hFDDA_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFDDA}}})});
      histPointers.insert({histPath + "hFDDC_Collisions_Unequalized", histos.add((histPath + "hFDDC_Collisions_Unequalized").c_str(), "hFDDC_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFDDC}}})});

      histPointers.insert({histPath + "hFV0A_Collisions_Unequalized", histos.add((histPath + "hFV0A_Collisions_Unequalized").c_str(), "hFV0A_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFV0A}}})});
      histPointers.insert({histPath + "hNGlobalTracks_Unequalized", histos.add((histPath + "hNGlobalTracks_Unequalized").c_str(), "hNGlobalTracks_Unequalized", {kTH1D, {{axisMultUltraFineGlobalTracks}}})});
      histPointers.insert({histPath + "hNMFTTracks_Unequalized", histos.add((histPath + "hNMFTTracks_Unequalized").c_str(), "hNMFTTracks_Unequalized", {kTH1D, {{axisMultUltraFineMFTTracks}}})});
      histPointers.insert({histPath + "hNPVContributors_Unequalized", histos.add((histPath + "hNPVContributors_Unequalized").c_str(), "hNPVContributors_Unequalized", {kTH1D, {{axisMultUltraFinePVContributors}}})});
    }

    if (doCentralityQA.value) {
      histPointers.insert({histPath + "hCentralityDistributionFV0A", histos.add((histPath + "hCentralityDistributionFV0A").c_str(), "hCentralityDistributionFV0A", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionFT0A", histos.add((histPath + "hCentralityDistributionFT0A").c_str(), "hCentralityDistributionFT0A", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionFT0C", histos.add((histPath + "hCentralityDistributionFT0C").c_str(), "hCentralityDistributionFT0C", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionFT0M", histos.add((histPath + "hCentralityDistributionFT0M").c_str(), "hCentralityDistributionFT0M", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionFDDM", histos.add((histPath + "hCentralityDistributionFDDM").c_str(), "hCentralityDistributionFDDM", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionNPTV", histos.add((histPath + "hCentralityDistributionNPTV").c_str(), "hCentralityDistributionNPTV", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionNGlobal", histos.add((histPath + "hCentralityDistributionNGlobal").c_str(), "hCentralityDistributionNGlobal", {kTH1D, {{axisCentrality}}})});
      histPointers.insert({histPath + "hCentralityDistributionMFT", histos.add((histPath + "hCentralityDistributionMFT").c_str(), "hCentralityDistributionMFT", {kTH1D, {{axisCentrality}}})});
    }

    histPointers.insert({histPath + "hFT0AvsPVz_Collisions", histos.add((histPath + "hFT0AvsPVz_Collisions").c_str(), "hFT0AvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hFT0CvsPVz_Collisions", histos.add((histPath + "hFT0CvsPVz_Collisions").c_str(), "hFT0CvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hFDDAvsPVz_Collisions", histos.add((histPath + "hFDDAvsPVz_Collisions").c_str(), "hFDDAvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hFDDCvsPVz_Collisions", histos.add((histPath + "hFDDCvsPVz_Collisions").c_str(), "hFDDCvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hFV0AvsPVz_Collisions", histos.add((histPath + "hFV0AvsPVz_Collisions").c_str(), "hFV0AvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hNGlobalTracksvsPVz_Collisions", histos.add((histPath + "hNGlobalTracksvsPVz_Collisions").c_str(), "hNGlobalTracksvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hNMFTTracksvsPVz_Collisions", histos.add((histPath + "hNMFTTracksvsPVz_Collisions").c_str(), "hNMFTTracksvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    histPointers.insert({histPath + "hNTPVvsPVz_Collisions", histos.add((histPath + "hNTPVvsPVz_Collisions").c_str(), "hNTPVvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
  }

  template <typename TCollision>
  void genericProcessCollision(const TCollision& collision)
  // process this collisions
  {
    initRun(collision);
    histos.fill(HIST("hCollisionSelection"), 0); // all collisions
    getHist(TH1, histPath + "hCollisionSelection")->Fill(0);

    if (applySel8 && !collision.multSel8()) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 1);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(1);

    // calculate vertex-Z-equalized quantities if desired
    float multFV0A = collision.multFV0A();
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFDDA = collision.multFDDA();
    float multFDDC = collision.multFDDC();
    float multNTracksGlobal = collision.multNTracksGlobal();
    float mftNtracks = collision.mftNtracks();
    float multNTracksPV = collision.multNTracksPV();
    if (applyVertexZEqualization) {
      float epsilon = 1e-2; // average value after which this collision will be disregarded
      multFV0A = -1.0f;
      multFT0A = -1.0f;
      multFT0C = -1.0f;
      multFDDA = -1.0f;
      multFDDC = -1.0f;
      multNTracksGlobal = -1.0f;
      mftNtracks = -1.0f;
      multNTracksPV = -1.0f;

      if (hVtxZFV0A->Interpolate(collision.multPVz()) > epsilon) {
        multFV0A = hVtxZFV0A->Interpolate(0.0) * collision.multFV0A() / hVtxZFV0A->Interpolate(collision.multPVz());
      }
      if (hVtxZFT0A->Interpolate(collision.multPVz()) > epsilon) {
        multFT0A = hVtxZFT0A->Interpolate(0.0) * collision.multFT0A() / hVtxZFT0A->Interpolate(collision.multPVz());
      }
      if (hVtxZFT0C->Interpolate(collision.multPVz()) > epsilon) {
        multFT0C = hVtxZFT0C->Interpolate(0.0) * collision.multFT0C() / hVtxZFT0C->Interpolate(collision.multPVz());
      }
      if (hVtxZFDDA->Interpolate(collision.multPVz()) > epsilon) {
        multFDDA = hVtxZFDDA->Interpolate(0.0) * collision.multFDDA() / hVtxZFDDA->Interpolate(collision.multPVz());
      }
      if (hVtxZFDDC->Interpolate(collision.multPVz()) > epsilon) {
        multFDDC = hVtxZFDDC->Interpolate(0.0) * collision.multFDDC() / hVtxZFDDC->Interpolate(collision.multPVz());
      }
      if (hVtxZNGlobals->Interpolate(collision.multPVz()) > epsilon) {
        multNTracksGlobal = hVtxZNGlobals->Interpolate(0.0) * collision.multNTracksGlobal() / hVtxZNGlobals->Interpolate(collision.multPVz());
      }
      if (hVtxZMFT->Interpolate(collision.multPVz()) > epsilon) {
        mftNtracks = hVtxZMFT->Interpolate(0.0) * collision.mftNtracks() / hVtxZMFT->Interpolate(collision.multPVz());
      }
      if (hVtxZNTracks->Interpolate(collision.multPVz()) > epsilon) {
        multNTracksPV = hVtxZNTracks->Interpolate(0.0) * collision.multNTracksPV() / hVtxZNTracks->Interpolate(collision.multPVz());
      }
    }

    bool passRejectITSROFBorder = !(rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder));
    bool passRejectTFBorder = !(rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder));
    bool passRequireIsVertexITSTPC = !(requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC));
    bool passRequireIsGoodZvtxFT0VsPV = !(requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV));
    bool passRequireIsVertexTOFmatched = !(requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched));
    bool passRequireIsVertexTRDmatched = !(requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched));
    bool passRejectSameBunchPileup = !(rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    bool passINELgtZERO = !(requireINELgtZERO && !collision.isInelGt0());

    // _______________________________________________________
    // sidestep vertex-Z rejection for vertex-Z profile histograms
    if (passRejectITSROFBorder && passRejectTFBorder && passRequireIsVertexITSTPC && passRequireIsGoodZvtxFT0VsPV &&
        passRequireIsVertexTOFmatched && passRequireIsVertexTRDmatched && passRejectSameBunchPileup && passINELgtZERO) {
      // all runs
      histos.fill(HIST("hFT0AvsPVz_Collisions"), collision.multPVz(), collision.multFT0A());
      histos.fill(HIST("hFT0CvsPVz_Collisions"), collision.multPVz(), collision.multFT0C());
      histos.fill(HIST("hFDDAvsPVz_Collisions"), collision.multPVz(), collision.multFDDA());
      histos.fill(HIST("hFDDCvsPVz_Collisions"), collision.multPVz(), collision.multFDDC());
      histos.fill(HIST("hFV0AvsPVz_Collisions"), collision.multPVz(), collision.multFV0A());
      histos.fill(HIST("hNGlobalTracksvsPVz_Collisions"), collision.multPVz(), collision.multNTracksGlobal());
      histos.fill(HIST("hNMFTTracksvsPVz_Collisions"), collision.multPVz(), collision.mftNtracks());

      // per run
      getHist(TProfile, histPath + "hFT0CvsPVz_Collisions")->Fill(collision.multPVz(), multFT0C);
      getHist(TProfile, histPath + "hFT0AvsPVz_Collisions")->Fill(collision.multPVz(), multFT0A);
      getHist(TProfile, histPath + "hFDDAvsPVz_Collisions")->Fill(collision.multPVz(), multFDDA);
      getHist(TProfile, histPath + "hFDDCvsPVz_Collisions")->Fill(collision.multPVz(), multFDDC);
      getHist(TProfile, histPath + "hFV0AvsPVz_Collisions")->Fill(collision.multPVz(), multFV0A);
      getHist(TProfile, histPath + "hNGlobalTracksvsPVz_Collisions")->Fill(collision.multPVz(), multNTracksGlobal);
      getHist(TProfile, histPath + "hNMFTTracksvsPVz_Collisions")->Fill(collision.multPVz(), mftNtracks);
      getHist(TProfile, histPath + "hNTPVvsPVz_Collisions")->Fill(collision.multPVz(), multNTracksPV);
    }

    // _______________________________________________________

    if (applyVtxZ && TMath::Abs(collision.multPVz()) > vertexZcut.value) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 2);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(2);

    // _______________________________________________________
    // Extra event selections start here
    if (!passRejectITSROFBorder) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 3 /* Not at ITS ROF border */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(3);

    if (!passRejectTFBorder) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 4 /* Not at TF border */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(4);

    if (!passRequireIsVertexITSTPC) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 5 /* Contains at least one ITS-TPC track */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(5);

    if (!passRequireIsGoodZvtxFT0VsPV) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 6 /* PV position consistency check */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(6);

    if (!passRequireIsVertexTOFmatched) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 7 /* PV with at least one contributor matched with TOF */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(7);

    if (!passRequireIsVertexTRDmatched) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 8 /* PV with at least one contributor matched with TRD */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(8);

    if (!passRejectSameBunchPileup) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 9 /* Not at same bunch pile-up */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(9);

    if (!passINELgtZERO) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 11 /* is INEL > 0 */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(11);

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hNPVContributors"), collision.multNTracksPV());
    histos.fill(HIST("hFT0A_Collisions"), collision.multFT0A());
    histos.fill(HIST("hFT0C_Collisions"), collision.multFT0C());
    histos.fill(HIST("hFT0M_Collisions"), (collision.multFT0A() + collision.multFT0C()));
    histos.fill(HIST("hFDDA_Collisions"), collision.multFDDA());
    histos.fill(HIST("hFDDC_Collisions"), collision.multFDDC());
    histos.fill(HIST("hFV0A_Collisions"), collision.multFV0A());
    histos.fill(HIST("hNGlobalTracks"), collision.multNTracksGlobal());
    histos.fill(HIST("hNMFTTracks"), collision.mftNtracks());

    // save vertex-Z equalized
    getHist(TH1, histPath + "hNPVContributors")->Fill(multNTracksPV);
    getHist(TH1, histPath + "hFT0A_Collisions")->Fill(multFT0A);
    getHist(TH1, histPath + "hFT0C_Collisions")->Fill(multFT0C);
    getHist(TH1, histPath + "hFT0M_Collisions")->Fill((multFT0A + multFT0C));
    getHist(TH1, histPath + "hFV0A_Collisions")->Fill(multFV0A);
    getHist(TH1, histPath + "hNGlobalTracks")->Fill(multNTracksGlobal);
    getHist(TH1, histPath + "hNMFTTracks")->Fill(mftNtracks);

    if (applyVertexZEqualization.value && saveUnequalized.value) {
      // save unequalized for cross-checks
      getHist(TH1, histPath + "hNPVContributors_Unequalized")->Fill(collision.multNTracksPV());
      getHist(TH1, histPath + "hFT0A_Collisions_Unequalized")->Fill(collision.multFT0A());
      getHist(TH1, histPath + "hFT0C_Collisions_Unequalized")->Fill(collision.multFT0C());
      getHist(TH1, histPath + "hFT0M_Collisions_Unequalized")->Fill((collision.multFT0A() + collision.multFT0C()));
      getHist(TH1, histPath + "hFDDA_Collisions_Unequalized")->Fill(collision.multFDDA());
      getHist(TH1, histPath + "hFDDC_Collisions_Unequalized")->Fill(collision.multFDDC());
      getHist(TH1, histPath + "hFV0A_Collisions_Unequalized")->Fill(collision.multFV0A());
      getHist(TH1, histPath + "hNGlobalTracks_Unequalized")->Fill(collision.multNTracksGlobal());
      getHist(TH1, histPath + "hNMFTTracks_Unequalized")->Fill(collision.mftNtracks());
    }

    if (doCentralityQA.value) {
      // generate centralities on the spot in case the centrality histograms are in memory
      if (hCentralityFV0A) {
        getHist(TH1, histPath + "hCentralityDistributionFV0A")->Fill(hCentralityFV0A->GetBinContent(hCentralityFV0A->FindBin(multFV0A)));
      }
      if (hCentralityFT0A) {
        getHist(TH1, histPath + "hCentralityDistributionFT0A")->Fill(hCentralityFT0A->GetBinContent(hCentralityFT0A->FindBin(multFT0A)));
      }
      if (hCentralityFT0C) {
        getHist(TH1, histPath + "hCentralityDistributionFT0C")->Fill(hCentralityFT0C->GetBinContent(hCentralityFT0C->FindBin(multFT0C)));
      }
      if (hCentralityFT0M) {
        getHist(TH1, histPath + "hCentralityDistributionFT0M")->Fill(hCentralityFT0M->GetBinContent(hCentralityFT0M->FindBin(multFT0A + multFT0C)));
      }
      if (hCentralityFDDM) {
        getHist(TH1, histPath + "hCentralityDistributionFDDM")->Fill(hCentralityFDDM->GetBinContent(hCentralityFDDM->FindBin(multFDDA + multFDDC)));
      }
      if (hCentralityNTPV) {
        // note: not vertex-Z-equalized to match what is done in central framework for this estimator
        getHist(TH1, histPath + "hCentralityDistributionNTPV")->Fill(hCentralityNTPV->GetBinContent(hCentralityNTPV->FindBin(collision.multNTracksPV())));
      }
      if (hCentralityNGlo) {
        getHist(TH1, histPath + "hCentralityDistributionNGlobal")->Fill(hCentralityNGlo->GetBinContent(hCentralityNGlo->FindBin(multNTracksGlobal)));
      }
      if (hCentralityMFT) {
        getHist(TH1, histPath + "hCentralityDistributionMFT")->Fill(hCentralityMFT->GetBinContent(hCentralityMFT->FindBin(mftNtracks)));
      }
    }
  }

  void process(soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultsGlobal, aod::MultSelections>::iterator const& collision)
  {
    genericProcessCollision(collision);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<centralityStudypp>(cfgc)};
}
