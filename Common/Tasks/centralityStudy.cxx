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
// This task does dedicated centrality studies for understanding the
// Run 3 Pb-Pb centrality selections in 2023 data. It is compatible with
// derived data.

#include "EventSelectionParams.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPECSObject.h>
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

struct centralityStudy {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;
  std::string histPath;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher mRateFetcher;
  int mRunNumber;
  uint64_t startOfRunTimestamp;

  // vertex Z equalization
  TList* hCalibObjects;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZNTracks;
  TProfile* hVtxZNGlobals;
  TProfile* hVtxZMFT;
  TProfile* hVtxZFDDA;
  TProfile* hVtxZFDDC;

  // Configurables
  Configurable<bool> applyVertexZEqualization{"applyVertexZEqualization", false, "0 - no, 1 - yes"};

  Configurable<bool> do2DPlots{"do2DPlots", true, "0 - no, 1 - yes"};
  Configurable<bool> doOccupancyStudyVsCentrality2d{"doOccupancyStudyVsCentrality2d", true, "0 - no, 1 - yes"};
  Configurable<bool> doOccupancyStudyVsRawValues2d{"doOccupancyStudyVsRawValues2d", true, "0 - no, 1 - yes"};
  Configurable<bool> doOccupancyStudyVsCentrality3d{"doOccupancyStudyVsCentrality3d", false, "0 - no, 1 - yes"};
  Configurable<bool> doOccupancyStudyVsRawValues3d{"doOccupancyStudyVsRawValues3d", false, "0 - no, 1 - yes"};
  Configurable<bool> doTimeStudies{"doTimeStudies", true, "0 - no, 1 - yes"};
  Configurable<bool> doTimeStudyFV0AOuterVsFT0A3d{"doTimeStudyFV0AOuterVsFT0A3d", false, "0 - no, 1 - yes"};
  Configurable<bool> doNGlobalTracksVsRawSignals{"doNGlobalTracksVsRawSignals", true, "0 - no, 1 - yes"};
  Configurable<bool> applySel8{"applySel8", true, "0 - no, 1 - yes"};
  Configurable<bool> applyVtxZ{"applyVtxZ", true, "0 - no, 1 - yes"};

  // Apply extra event selections
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", false, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};

  Configurable<bool> rejectITSinROFpileupStandard{"rejectITSinROFpileupStandard", false, "reject collisions in case of in-ROF ITS pileup (standard)"};
  Configurable<bool> rejectITSinROFpileupStrict{"rejectITSinROFpileupStrict", false, "reject collisions in case of in-ROF ITS pileup (strict)"};
  Configurable<bool> rejectCollInTimeRangeNarrow{"rejectCollInTimeRangeNarrow", false, "reject if extra colls in time range (narrow)"};

  Configurable<bool> selectUPCcollisions{"selectUPCcollisions", false, "select collisions tagged with UPC flag"};

  Configurable<bool> selectCollidingBCs{"selectCollidingBCs", true, "BC analysis: select colliding BCs"};
  Configurable<bool> selectTVX{"selectTVX", true, "BC analysis: select TVX"};
  Configurable<bool> selectFV0OrA{"selectFV0OrA", true, "BC analysis: select FV0OrA"};
  Configurable<float> vertexZwithT0{"vertexZwithT0", 1000.0f, "require a certain vertex-Z in BC analysis"};

  Configurable<float> minTimeDelta{"minTimeDelta", -1.0f, "reject collision if another collision is this close or less in time"};
  Configurable<float> minFT0CforVertexZ{"minFT0CforVertexZ", -1.0f, "minimum FT0C for vertex-Z profile calculation"};

  Configurable<float> scaleSignalFT0A{"scaleSignalFT0A", 1.00f, "scale FT0A signal for convenience"};
  Configurable<float> scaleSignalFT0C{"scaleSignalFT0C", 1.00f, "scale FT0C signal for convenience"};
  Configurable<float> scaleSignalFT0M{"scaleSignalFT0M", 1.00f, "scale FT0M signal for convenience"};
  Configurable<float> scaleSignalFV0A{"scaleSignalFV0A", 1.00f, "scale FV0A signal for convenience"};

  Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> pathGRPECSObject{"pathGRPECSObject", "GLO/Config/GRPECS", "Path to GRPECS object"};
  Configurable<std::string> pathVertexZ{"pathVertexZ", "Users/d/ddobrigk/Centrality/Calibration", "Path to vertexZ profiles"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Source of the interaction rate: (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  Configurable<bool> irCrashOnNull{"irCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<bool> irDoRateVsTime{"irDoRateVsTime", true, "Do IR plots"};

  // _______________________________________
  // upc rejection criteria
  // reject low zna/c
  struct : ConfigurableGroup {
    Configurable<float> minZNACsignal{"minZNACsignal", -999999.0f, "min zna/c signal"};
    Configurable<float> maxFT0CforZNACselection{"maxFT0CforZNACselection", -99999.0f, "max ft0c signal for minZNACsignal to work"};

    Configurable<float> minFV0Asignal{"minFV0Asignal", -999999.0f, "min fv0a signal"};
    Configurable<float> maxFT0CforFV0Aselection{"maxFT0CforFV0Aselection", -99999.0f, "max ft0c signal for minFV0Asignal to work"};

    Configurable<float> minFDDAsignal{"minFDDAsignal", -999999.0f, "min fdda signal"};
    Configurable<float> maxFT0CforFDDAselection{"maxFT0CforFDDAselection", -99999.0f, "max ft0c signal for minFDDAsignal to work"};
  } upcRejection;

  // Configurable Axes for 2d plots, etc
  ConfigurableAxis axisMultFV0A{"axisMultFV0A", {1000, 0, 100000}, "FV0A amplitude"};
  ConfigurableAxis axisMultFT0A{"axisMultFT0A", {1000, 0, 100000}, "FT0A amplitude"};
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {1000, 0, 100000}, "FT0C amplitude"};
  ConfigurableAxis axisMultFT0M{"axisMultFT0M", {1000, 0, 100000}, "FT0M amplitude"};
  ConfigurableAxis axisMultFDDA{"axisMultFDDA", {1000, 0, 100000}, "FDDA amplitude"};
  ConfigurableAxis axisMultFDDC{"axisMultFDDC", {1000, 0, 100000}, "FDDC amplitude"};
  ConfigurableAxis axisMultPVContributors{"axisMultPVContributors", {200, 0, 6000}, "Number of PV Contributors"};
  ConfigurableAxis axisMultGlobalTracks{"axisMultGlobalTracks", {500, 0, 5000}, "Number of global tracks"};
  ConfigurableAxis axisMultMFTTracks{"axisMultMFTTracks", {500, 0, 5000}, "Number of MFT tracks"};
  ConfigurableAxis axisMultMCCounts{"axisMultMCCounts", {1000, 0, 5000}, "N_{ch}"};

  ConfigurableAxis axisTrackOccupancy{"axisTrackOccupancy", {50, 0, 5000}, "Track occupancy"};
  ConfigurableAxis axisFT0COccupancy{"axisFT0COccupancy", {50, 0, 80000}, "FT0C occupancy"};

  // For one-dimensional plots, where binning is no issue
  ConfigurableAxis axisMultUltraFineFV0A{"axisMultUltraFineFV0A", {60000, 0, 60000}, "FV0A amplitude"};
  ConfigurableAxis axisMultUltraFineFT0M{"axisMultUltraFineFT0M", {50000, 0, 200000}, "FT0M amplitude"};
  ConfigurableAxis axisMultUltraFineFT0C{"axisMultUltraFineFT0C", {60000, 0, 60000}, "FT0C amplitude"};
  ConfigurableAxis axisMultUltraFineFT0A{"axisMultUltraFineFT0A", {60000, 0, 60000}, "FT0C amplitude"};
  ConfigurableAxis axisMultUltraFinePVContributors{"axisMultUltraFinePVContributors", {10000, 0, 10000}, "Number of PV Contributors"};
  ConfigurableAxis axisMultUltraFineGlobalTracks{"axisMultUltraFineGlobalTracks", {5000, 0, 5000}, "Number of global tracks"};
  ConfigurableAxis axisMultUltraFineMFTTracks{"axisMultUltraFineMFTTracks", {5000, 0, 5000}, "Number of MFT tracks"};

  ConfigurableAxis axisMultITSOnly{"axisMultITSOnly", {200, 0, 6000}, "Number of ITS only tracks"};
  ConfigurableAxis axisMultITSTPC{"axisMultITSTPC", {200, 0, 6000}, "Number of ITSTPC matched tracks"};

  // For centrality studies if requested
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0, 100}, "FT0C percentile"};
  ConfigurableAxis axisImpactParameter{"axisImpactParameter", {200, 0.0f, 20.0f}, "b (fm)"};
  ConfigurableAxis axisPVChi2{"axisPVChi2", {300, 0, 30}, "FT0C percentile"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {300, 0, 300}, "#Delta time"};

  ConfigurableAxis axisDeltaTimestamp{"axisDeltaTimestamp", {1440, 0, 24}, "#Delta timestamp - sor (hours)"};
  ConfigurableAxis axisInteractionRate{"axisInteractionRate", {500, 0, 100}, "Binning for the interaction rate (kHz)"};
  ConfigurableAxis axisMultCoarseFV0A{"axisMultCoarseFV0A", {350, 0, 70000}, "FV0A amplitude"};

  // For profile Z
  ConfigurableAxis axisPVz{"axisPVz", {400, -20.0f, +20.0f}, "PVz (cm)"};
  ConfigurableAxis axisZN{"axisZN", {1100, -50.0f, +500.0f}, "ZN"};

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

    if (doprocessCollisions || doprocessCollisionsWithCentrality) {
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

      histos.add("hFT0A_Collisions", "hFT0C_Collisions", kTH1D, {axisMultUltraFineFT0A});
      histos.add("hFT0C_Collisions", "hFT0C_Collisions", kTH1D, {axisMultUltraFineFT0C});
      histos.add("hFT0M_Collisions", "hFT0M_Collisions", kTH1D, {axisMultUltraFineFT0M});
      histos.add("hFV0A_Collisions", "hFV0A_Collisions", kTH1D, {axisMultUltraFineFV0A});
      histos.add("hNGlobalTracks", "hNGlobalTracks", kTH1D, {axisMultUltraFineGlobalTracks});
      histos.add("hNMFTTracks", "hNMFTTracks", kTH1D, {axisMultUltraFineMFTTracks});
      histos.add("hNPVContributors", "hNPVContributors", kTH1D, {axisMultUltraFinePVContributors});

      histos.add("hFT0CvsPVz_Collisions_All", "hFT0CvsPVz_Collisions_All", kTProfile, {axisPVz});
      histos.add("hFT0CvsPVz_Collisions", "hFT0CvsPVz_Collisions", kTProfile, {axisPVz});
      histos.add("hFV0AvsPVz_Collisions", "hFV0AvsPVz_Collisions", kTProfile, {axisPVz});
      histos.add("hNGlobalTracksvsPVz_Collisions", "hNGlobalTracksvsPVz_Collisions", kTProfile, {axisPVz});
      histos.add("hNMFTTracksvsPVz_Collisions", "hNMFTTracksvsPVz_Collisions", kTProfile, {axisPVz});
    }

    if (doprocessBCs) {
      histos.add("hBCSelection", "hBCSelection", kTH1D, {{20, -0.5, 19.5f}});
      histos.add("hFT0C_BCs", "hFT0C_BCs", kTH1D, {axisMultUltraFineFT0C});
      histos.add("hFT0M_BCs", "hFT0M_BCs", kTH1D, {axisMultUltraFineFT0M});
      histos.add("hFV0A_BCs", "hFV0A_BCs", kTH1D, {axisMultUltraFineFV0A});
      histos.add("hFT0CvsPVz_BCs_All", "hFT0CvsPVz_BCs_All", kTProfile, {axisPVz});
      histos.add("hFT0CvsPVz_BCs", "hFT0CvsPVz_BCs", kTProfile, {axisPVz});
      histos.add("hVertexZ_BCvsCO", "hVertexZ_BCvsCO", kTH2F, {axisPVz, axisPVz});
      histos.add("hZNAvsFT0C_BCs", "hZNAvsFT0C_BCs", kTH2F, {axisMultFT0C, axisZN});
      histos.add("hZNCvsFT0C_BCs", "hZNCvsFT0C_BCs", kTH2F, {axisMultFT0C, axisZN});
    }

    if (do2DPlots) {
      histos.add("hNContribsVsFT0C", "hNContribsVsFT0C", kTH2F, {axisMultFT0C, axisMultPVContributors});
      histos.add("hNContribsVsFV0A", "hNContribsVsFV0A", kTH2F, {axisMultFV0A, axisMultPVContributors});
      histos.add("hMatchedVsITSOnly", "hMatchedVsITSOnly", kTH2F, {axisMultITSOnly, axisMultITSTPC});

      // 2d correlation of fit signals
      histos.add("hFT0AVsFT0C", "hFT0AVsFT0C", kTH2F, {axisMultFT0C, axisMultFT0A});
      histos.add("hFV0AVsFT0C", "hFV0AVsFT0C", kTH2F, {axisMultFT0C, axisMultFV0A});
      histos.add("hFDDAVsFT0C", "hFDDAVsFT0C", kTH2F, {axisMultFT0C, axisMultFDDA});
      histos.add("hFDDCVsFT0C", "hFDDCVsFT0C", kTH2F, {axisMultFT0C, axisMultFDDC});
    }

    if (doNGlobalTracksVsRawSignals) {
      histos.add("hNGlobalTracksVsFT0A", "hNGlobalTracksVsFT0A", kTH2F, {axisMultFT0A, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsFT0C", "hNGlobalTracksVsFT0C", kTH2F, {axisMultFT0C, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsFT0M", "hNGlobalTracksVsFT0M", kTH2F, {axisMultFT0M, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsFV0A", "hNGlobalTracksVsFV0A", kTH2F, {axisMultFV0A, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsFDDA", "hNGlobalTracksVsFDDA", kTH2F, {axisMultFDDA, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsFDDC", "hNGlobalTracksVsFDDC", kTH2F, {axisMultFDDC, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsZNA", "hNGlobalTracksVsZNA", kTH2F, {axisZN, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsZNC", "hNGlobalTracksVsZNC", kTH2F, {axisZN, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsNMFTTracks", "hNGlobalTracksVsNMFTTracks", kTH2F, {axisMultMFTTracks, axisMultGlobalTracks});
      histos.add("hNGlobalTracksVsNTPV", "hNGlobalTracksVsNTPV", kTH2F, {axisMultPVContributors, axisMultGlobalTracks});
    }

    if (doprocessCollisionsWithResolutionStudy) {
      // histograms with detector signals
      histos.add("hImpactParameterVsFT0A", "hImpactParameterVsFT0A", kTH2F, {axisMultFT0A, axisImpactParameter});
      histos.add("hImpactParameterVsFT0C", "hImpactParameterVsFT0C", kTH2F, {axisMultFT0C, axisImpactParameter});
      histos.add("hImpactParameterVsFT0M", "hImpactParameterVsFT0M", kTH2F, {axisMultFT0M, axisImpactParameter});
      histos.add("hImpactParameterVsFV0A", "hImpactParameterVsFV0A", kTH2F, {axisMultFV0A, axisImpactParameter});
      histos.add("hImpactParameterVsNMFTTracks", "hImpactParameterVsNMFTTracks", kTH2F, {axisMultMFTTracks, axisImpactParameter});
      histos.add("hImpactParameterVsNTPV", "hImpactParameterVsNTPV", kTH2F, {axisMultPVContributors, axisImpactParameter});

      // histograms with actual MC counts in each region
      histos.add("hImpactParameterVsMCFT0A", "hImpactParameterVsMCFT0A", kTH2F, {axisMultMCCounts, axisImpactParameter});
      histos.add("hImpactParameterVsMCFT0C", "hImpactParameterVsMCFT0C", kTH2F, {axisMultMCCounts, axisImpactParameter});
      histos.add("hImpactParameterVsMCFT0M", "hImpactParameterVsMCFT0M", kTH2F, {axisMultMCCounts, axisImpactParameter});
      histos.add("hImpactParameterVsMCFV0A", "hImpactParameterVsMCFV0A", kTH2F, {axisMultMCCounts, axisImpactParameter});
    }

    if (doOccupancyStudyVsRawValues2d) {
      histos.add("hNcontribsProfileVsTrackOccupancyVsFT0C", "hNcontribsProfileVsTrackOccupancyVsFT0C", kTProfile2D, {axisTrackOccupancy, axisMultFT0C});
      histos.add("hNGlobalTracksProfileVsTrackOccupancyVsFT0C", "hNGlobalTracksProfileVsTrackOccupancyVsFT0C", kTProfile2D, {axisTrackOccupancy, axisMultFT0C});
      histos.add("hNcontribsProfileVsFT0COccupancyVsFT0C", "hNcontribsProfileVsFT0COccupancyVsFT0C", kTProfile2D, {axisFT0COccupancy, axisMultFT0C});
      histos.add("hNGlobalTracksProfileVsFT0COccupancyVsFT0C", "hNGlobalTracksProfileVsFT0COccupancyVsFT0C", kTProfile2D, {axisFT0COccupancy, axisMultFT0C});
    }

    if (doOccupancyStudyVsRawValues3d) {
      histos.add("hTrackOccupancyVsNContribsVsFT0C", "hTrackOccupancyVsNContribsVsFT0C", kTH3F, {axisTrackOccupancy, axisMultPVContributors, axisMultFT0C});
      histos.add("hTrackOccupancyVsNGlobalTracksVsFT0C", "hTrackOccupancyVsNGlobalTracksVsFT0C", kTH3F, {axisTrackOccupancy, axisMultGlobalTracks, axisMultFT0C});
      histos.add("hFT0COccupancyVsNContribsVsFT0C", "hFT0COccupancyVsNContribsVsFT0C", kTH3F, {axisFT0COccupancy, axisMultPVContributors, axisMultFT0C});
      histos.add("hFT0COccupancyVsNGlobalTracksVsFT0C", "hFT0COccupancyVsNGlobalTracksVsFT0C", kTH3F, {axisFT0COccupancy, axisMultGlobalTracks, axisMultFT0C});
    }

    if (doprocessCollisionsWithCentrality) {
      // in case requested: do vs centrality debugging
      histos.add("hCentrality", "hCentrality", kTH1F, {axisCentrality});
      histos.add("hNContribsVsCentrality", "hNContribsVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNITSTPCTracksVsCentrality", "hNITSTPCTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNITSOnlyTracksVsCentrality", "hNITSOnlyTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNGlobalTracksVsCentrality", "hNGlobalTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNMFTTracksVsCentrality", "hNMFTTracksVsCentrality", kTH2F, {axisCentrality, axisMultMFTTracks});
      histos.add("hPVChi2VsCentrality", "hPVChi2VsCentrality", kTH2F, {axisCentrality, axisPVChi2});
      histos.add("hDeltaTimeVsCentrality", "hDeltaTimeVsCentrality", kTH2F, {axisCentrality, axisDeltaTime});

      if (doOccupancyStudyVsCentrality2d) {
        histos.add("hNcontribsProfileVsTrackOccupancyVsCentrality", "hNcontribsProfileVsTrackOccupancyVsCentrality", kTProfile2D, {axisTrackOccupancy, axisCentrality});
        histos.add("hNGlobalTracksProfileVsTrackOccupancyVsCentrality", "hNGlobalTracksProfileVsTrackOccupancyVsCentrality", kTProfile2D, {axisTrackOccupancy, axisCentrality});
        histos.add("hNcontribsProfileVsFT0COccupancyVsCentrality", "hNcontribsProfileVsFT0COccupancyVsCentrality", kTProfile2D, {axisFT0COccupancy, axisCentrality});
        histos.add("hNGlobalTracksProfileVsFT0COccupancyVsCentrality", "hNGlobalTracksProfileVsFT0COccupancyVsCentrality", kTProfile2D, {axisFT0COccupancy, axisCentrality});
      }

      if (doOccupancyStudyVsCentrality3d) {
        histos.add("hTrackOccupancyVsNContribsVsCentrality", "hTrackOccupancyVsNContribsVsCentrality", kTH3F, {axisTrackOccupancy, axisMultPVContributors, axisCentrality});
        histos.add("hTrackOccupancyVsNGlobalTracksVsCentrality", "hTrackOccupancyVsNGlobalTracksVsCentrality", kTH3F, {axisTrackOccupancy, axisMultGlobalTracks, axisCentrality});
        histos.add("hFT0COccupancyVsNContribsVsCentrality", "hFT0COccupancyVsNContribsVsCentrality", kTH3F, {axisFT0COccupancy, axisMultPVContributors, axisCentrality});
        histos.add("hFT0COccupancyVsNGlobalTracksVsCentrality", "hFT0COccupancyVsNGlobalTracksVsCentrality", kTH3F, {axisFT0COccupancy, axisMultGlobalTracks, axisCentrality});
      }
    }

    if (doTimeStudies) {
      ccdb->setURL(ccdbURL);
      // ccdb->setCaching(true);
      // ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      if (doTimeStudyFV0AOuterVsFT0A3d) {
        histos.add((histPath + "h3dFV0AVsTime").c_str(), "", {kTH3F, {{axisDeltaTimestamp, axisMultCoarseFV0A, axisMultCoarseFV0A}}});
      }
    }
  }

  template <typename TCollision>
  void initRun(const TCollision& collision)
  {
    if (mRunNumber == collision.multRunNumber()) {
      return;
    }

    mRunNumber = collision.multRunNumber();
    LOGF(info, "Setting up for run: %i", mRunNumber);

    // only get object when switching runs
    o2::parameters::GRPECSObject* grpo = ccdb->getForRun<o2::parameters::GRPECSObject>(pathGRPECSObject, mRunNumber);
    startOfRunTimestamp = grpo->getTimeStart();

    if (applyVertexZEqualization.value) {
      // acquire vertex-Z equalization histograms if requested
      LOGF(info, "Acquiring vertex-Z profiles for run %i", mRunNumber);
      hCalibObjects = ccdb->getForRun<TList>(pathVertexZ, mRunNumber);

      hVtxZFV0A = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFV0A"));
      hVtxZFT0A = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFT0A"));
      hVtxZFT0C = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFT0C"));
      // hVtxZFDDA = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFDDA"));
      // hVtxZFDDC = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZFDDC"));
      hVtxZNTracks = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZNTracksPV"));
      hVtxZNGlobals = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZNGlobals"));
      hVtxZMFT = static_cast<TProfile*>(hCalibObjects->FindObject("hVtxZMFT"));

      // Capture error
      if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZNTracks || !hVtxZNGlobals || !hVtxZMFT) {
        LOGF(error, "Problem loading CCDB objects! Please check");
      }
    }

    histPath = std::format("Run_{}/", mRunNumber);

    if (doprocessCollisions || doprocessCollisionsWithCentrality) {
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
      getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(11, "Neighbour rejection");
      getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(12, "no ITS in-ROF pileup (standard)");
      getHist(TH1, histPath + "hCollisionSelection")->GetXaxis()->SetBinLabel(13, "no ITS in-ROF pileup (strict)");

      histPointers.insert({histPath + "hFT0C_Collisions", histos.add((histPath + "hFT0C_Collisions").c_str(), "hFT0C_Collisions", {kTH1D, {{axisMultUltraFineFT0C}}})});
      histPointers.insert({histPath + "hFT0A_Collisions", histos.add((histPath + "hFT0A_Collisions").c_str(), "hFT0A_Collisions", {kTH1D, {{axisMultUltraFineFT0A}}})});
      histPointers.insert({histPath + "hFT0M_Collisions", histos.add((histPath + "hFT0M_Collisions").c_str(), "hFT0M_Collisions", {kTH1D, {{axisMultUltraFineFT0M}}})});
      histPointers.insert({histPath + "hFV0A_Collisions", histos.add((histPath + "hFV0A_Collisions").c_str(), "hFV0A_Collisions", {kTH1D, {{axisMultUltraFineFV0A}}})});
      histPointers.insert({histPath + "hNGlobalTracks", histos.add((histPath + "hNGlobalTracks").c_str(), "hNGlobalTracks", {kTH1D, {{axisMultUltraFineGlobalTracks}}})});
      histPointers.insert({histPath + "hNMFTTracks", histos.add((histPath + "hNMFTTracks").c_str(), "hNMFTTracks", {kTH1D, {{axisMultUltraFineMFTTracks}}})});
      histPointers.insert({histPath + "hNPVContributors", histos.add((histPath + "hNPVContributors").c_str(), "hNPVContributors", {kTH1D, {{axisMultUltraFinePVContributors}}})});

      if (applyVertexZEqualization) {
        histPointers.insert({histPath + "hFT0C_Collisions_Unequalized", histos.add((histPath + "hFT0C_Collisions_Unequalized").c_str(), "hFT0C_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFT0C}}})});
        histPointers.insert({histPath + "hFT0M_Collisions_Unequalized", histos.add((histPath + "hFT0M_Collisions_Unequalized").c_str(), "hFT0M_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFT0M}}})});
        histPointers.insert({histPath + "hFV0A_Collisions_Unequalized", histos.add((histPath + "hFV0A_Collisions_Unequalized").c_str(), "hFV0A_Collisions_Unequalized", {kTH1D, {{axisMultUltraFineFV0A}}})});
        histPointers.insert({histPath + "hNGlobalTracks_Unequalized", histos.add((histPath + "hNGlobalTracks_Unequalized").c_str(), "hNGlobalTracks_Unequalized", {kTH1D, {{axisMultUltraFineGlobalTracks}}})});
        histPointers.insert({histPath + "hNMFTTracks_Unequalized", histos.add((histPath + "hNMFTTracks_Unequalized").c_str(), "hNMFTTracks_Unequalized", {kTH1D, {{axisMultUltraFineMFTTracks}}})});
        histPointers.insert({histPath + "hNPVContributors_Unequalized", histos.add((histPath + "hNPVContributors_Unequalized").c_str(), "hNPVContributors_Unequalized", {kTH1D, {{axisMultUltraFinePVContributors}}})});
      }

      histPointers.insert({histPath + "hFT0CvsPVz_Collisions_All", histos.add((histPath + "hFT0CvsPVz_Collisions_All").c_str(), "hFT0CvsPVz_Collisions_All", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hFT0AvsPVz_Collisions", histos.add((histPath + "hFT0AvsPVz_Collisions").c_str(), "hFT0AvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hFT0CvsPVz_Collisions", histos.add((histPath + "hFT0CvsPVz_Collisions").c_str(), "hFT0CvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hFV0AvsPVz_Collisions", histos.add((histPath + "hFV0AvsPVz_Collisions").c_str(), "hFV0AvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hNGlobalTracksvsPVz_Collisions", histos.add((histPath + "hNGlobalTracksvsPVz_Collisions").c_str(), "hNGlobalTracksvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hNMFTTracksvsPVz_Collisions", histos.add((histPath + "hNMFTTracksvsPVz_Collisions").c_str(), "hNMFTTracksvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
      histPointers.insert({histPath + "hNTPVvsPVz_Collisions", histos.add((histPath + "hNTPVvsPVz_Collisions").c_str(), "hNTPVvsPVz_Collisions", {kTProfile, {{axisPVz}}})});
    }

    if (do2DPlots) {
      histPointers.insert({histPath + "hNContribsVsFT0C", histos.add((histPath + "hNContribsVsFT0C").c_str(), "hNContribsVsFT0C", {kTH2F, {{axisMultFT0C, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hNContribsVsFV0A", histos.add((histPath + "hNContribsVsFV0A").c_str(), "hNContribsVsFV0A", {kTH2F, {{axisMultFV0A, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hMatchedVsITSOnly", histos.add((histPath + "hMatchedVsITSOnly").c_str(), "hMatchedVsITSOnly", {kTH2F, {{axisMultITSOnly, axisMultITSTPC}}})});

      // 2d correlation of fit signals
      histPointers.insert({histPath + "hFT0AVsFT0C", histos.add((histPath + "hFT0AVsFT0C").c_str(), "hFT0AVsFT0C", {kTH2F, {{axisMultFT0C, axisMultFT0A}}})});
      histPointers.insert({histPath + "hFV0AVsFT0C", histos.add((histPath + "hFV0AVsFT0C").c_str(), "hFV0AVsFT0C", {kTH2F, {{axisMultFT0C, axisMultFV0A}}})});
      histPointers.insert({histPath + "hFDDAVsFT0C", histos.add((histPath + "hFDDAVsFT0C").c_str(), "hFDDAVsFT0C", {kTH2F, {{axisMultFT0C, axisMultFDDA}}})});
      histPointers.insert({histPath + "hFDDCVsFT0C", histos.add((histPath + "hFDDCVsFT0C").c_str(), "hFDDCVsFT0C", {kTH2F, {{axisMultFT0C, axisMultFDDC}}})});
    }

    if (doprocessCollisionsWithCentrality) {
      // in case requested: do vs centrality debugging
      histPointers.insert({histPath + "hCentrality", histos.add((histPath + "hCentrality").c_str(), "hCentrality", {kTH1F, {{axisCentrality}}})});
      histPointers.insert({histPath + "hNContribsVsCentrality", histos.add((histPath + "hNContribsVsCentrality").c_str(), "hNContribsVsCentrality", {kTH2F, {{axisCentrality, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hNITSTPCTracksVsCentrality", histos.add((histPath + "hNITSTPCTracksVsCentrality").c_str(), "hNITSTPCTracksVsCentrality", {kTH2F, {{axisCentrality, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hNITSOnlyTracksVsCentrality", histos.add((histPath + "hNITSOnlyTracksVsCentrality").c_str(), "hNITSOnlyTracksVsCentrality", {kTH2F, {{axisCentrality, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsCentrality", histos.add((histPath + "hNGlobalTracksVsCentrality").c_str(), "hNGlobalTracksVsCentrality", {kTH2F, {{axisCentrality, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hNMFTTracksVsCentrality", histos.add((histPath + "hNMFTTracksVsCentrality").c_str(), "hNMFTTracksVsCentrality", {kTH2F, {{axisCentrality, axisMultMFTTracks}}})});
      histPointers.insert({histPath + "hPVChi2VsCentrality", histos.add((histPath + "hPVChi2VsCentrality").c_str(), "hPVChi2VsCentrality", {kTH2F, {{axisCentrality, axisPVChi2}}})});
      histPointers.insert({histPath + "hDeltaTimeVsCentrality", histos.add((histPath + "hDeltaTimeVsCentrality").c_str(), "hDeltaTimeVsCentrality", {kTH2F, {{axisCentrality, axisDeltaTime}}})});
    }

    if (doNGlobalTracksVsRawSignals) {
      histPointers.insert({histPath + "hNGlobalTracksVsFT0A", histos.add((histPath + "hNGlobalTracksVsFT0A").c_str(), "hNGlobalTracksVsFT0A", {kTH2F, {{axisMultFT0A, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsFT0C", histos.add((histPath + "hNGlobalTracksVsFT0C").c_str(), "hNGlobalTracksVsFT0C", {kTH2F, {{axisMultFT0C, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsFT0M", histos.add((histPath + "hNGlobalTracksVsFT0M").c_str(), "hNGlobalTracksVsFT0M", {kTH2F, {{axisMultFT0M, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsFV0A", histos.add((histPath + "hNGlobalTracksVsFV0A").c_str(), "hNGlobalTracksVsFV0A", {kTH2F, {{axisMultFV0A, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsNMFTTracks", histos.add((histPath + "hNGlobalTracksVsNMFTTracks").c_str(), "hNGlobalTracksVsNMFTTracks", {kTH2F, {{axisMultMFTTracks, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNGlobalTracksVsNTPV", histos.add((histPath + "hNGlobalTracksVsNTPV").c_str(), "hNGlobalTracksVsNTPV", {kTH2F, {{axisMultPVContributors, axisMultGlobalTracks}}})});
    }

    if (doTimeStudies) {
      histPointers.insert({histPath + "hFT0AVsTime", histos.add((histPath + "hFT0AVsTime").c_str(), "hFT0AVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultFT0A}}})});
      histPointers.insert({histPath + "hFT0CVsTime", histos.add((histPath + "hFT0CVsTime").c_str(), "hFT0CVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultFT0C}}})});
      histPointers.insert({histPath + "hFT0MVsTime", histos.add((histPath + "hFT0MVsTime").c_str(), "hFT0MVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultFT0M}}})});
      histPointers.insert({histPath + "hFV0AVsTime", histos.add((histPath + "hFV0AVsTime").c_str(), "hFV0AVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultFV0A}}})});
      histPointers.insert({histPath + "hFV0AOuterVsTime", histos.add((histPath + "hFV0AOuterVsTime").c_str(), "hFV0AOuterVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultFV0A}}})});
      histPointers.insert({histPath + "hMFTTracksVsTime", histos.add((histPath + "hMFTTracksVsTime").c_str(), "hMFTTracksVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultMFTTracks}}})});
      histPointers.insert({histPath + "hNGlobalVsTime", histos.add((histPath + "hNGlobalVsTime").c_str(), "hNGlobalVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultGlobalTracks}}})});
      histPointers.insert({histPath + "hNTPVContributorsVsTime", histos.add((histPath + "hNTPVContributorsVsTime").c_str(), "hNTPVContributorsVsTime", {kTH2F, {{axisDeltaTimestamp, axisMultPVContributors}}})});
      histPointers.insert({histPath + "hPVzProfileCoVsTime", histos.add((histPath + "hPVzProfileCoVsTime").c_str(), "hPVzProfileCoVsTime", {kTProfile, {{axisDeltaTimestamp}}})});
      histPointers.insert({histPath + "hPVzProfileBcVsTime", histos.add((histPath + "hPVzProfileBcVsTime").c_str(), "hPVzProfileBcVsTime", {kTProfile, {{axisDeltaTimestamp}}})});
      if (irDoRateVsTime) {
        histPointers.insert({histPath + "hIRProfileVsTime", histos.add((histPath + "hIRProfileVsTime").c_str(), "hIRProfileVsTime", {kTProfile, {{axisDeltaTimestamp}}})});
      }
    }
  }

  template <typename TCollision>
  void genericProcessCollision(const TCollision& collision)
  // process this collisions
  {
    initRun(collision);
    histos.fill(HIST("hCollisionSelection"), 0); // all collisions
    getHist(TH1, histPath + "hCollisionSelection")->Fill(0);

    if (applySel8 && !collision.multSel8())
      return;
    histos.fill(HIST("hCollisionSelection"), 1);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(1);

    // calculate vertex-Z-equalized quantities if desired
    float multFV0A = collision.multFV0A();
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multNTracksGlobal = collision.multNTracksGlobal();
    float mftNtracks = collision.mftNtracks();
    float multNTracksPV = collision.multNTracksPV();
    if (applyVertexZEqualization) {
      float epsilon = 1e-2; // average value after which this collision will be disregarded
      multFV0A = -1.0f;
      multFT0A = -1.0f;
      multFT0C = -1.0f;
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
    bool passRejectITSinROFpileupStandard = !(rejectITSinROFpileupStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard));
    bool passRejectITSinROFpileupStrict = !(rejectITSinROFpileupStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict));
    bool passSelectUPCcollisions = !(selectUPCcollisions && collision.flags() < 1);
    bool passRejectCollInTimeRangeNarrow = !(rejectCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow));
    // _______________________________________________________
    // sidestep vertex-Z rejection for vertex-Z profile histograms
    if (passRejectITSROFBorder && passRejectTFBorder && passRequireIsVertexITSTPC && passRequireIsGoodZvtxFT0VsPV &&
        passRequireIsVertexTOFmatched && passRequireIsVertexTRDmatched && passRejectSameBunchPileup && passRejectITSinROFpileupStandard && passRejectITSinROFpileupStrict &&
        passSelectUPCcollisions && passRejectCollInTimeRangeNarrow) {
      getHist(TProfile, histPath + "hFT0CvsPVz_Collisions_All")->Fill(collision.multPVz(), multFT0C * scaleSignalFT0C);
      getHist(TProfile, histPath + "hFT0CvsPVz_Collisions")->Fill(collision.multPVz(), multFT0C * scaleSignalFT0C);
      getHist(TProfile, histPath + "hFT0AvsPVz_Collisions")->Fill(collision.multPVz(), multFT0A * scaleSignalFT0C);
      getHist(TProfile, histPath + "hFV0AvsPVz_Collisions")->Fill(collision.multPVz(), multFV0A * scaleSignalFV0A);
      getHist(TProfile, histPath + "hNGlobalTracksvsPVz_Collisions")->Fill(collision.multPVz(), multNTracksGlobal);
      getHist(TProfile, histPath + "hNMFTTracksvsPVz_Collisions")->Fill(collision.multPVz(), mftNtracks);
      getHist(TProfile, histPath + "hNTPVvsPVz_Collisions")->Fill(collision.multPVz(), multNTracksPV);
    }

    // _______________________________________________________

    if (applyVtxZ && TMath::Abs(collision.multPVz()) > 10)
      return;
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

    // do this only if information is available
    if constexpr (requires { collision.timeToNext(); }) {
      float timeToNeighbour = TMath::Min(
        std::abs(collision.timeToNext()),
        std::abs(collision.timeToPrevious()));
      histos.fill(HIST("hDeltaTimeVsCentrality"), collision.centFT0C(), timeToNeighbour);
      if (timeToNeighbour < minTimeDelta) {
        return;
      }
      histos.fill(HIST("hCollisionSelection"), 10 /* has suspicious neighbour */);
      getHist(TH1, histPath + "hCollisionSelection")->Fill(10);
    }

    if (!passRejectITSinROFpileupStandard) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 11 /* Not ITS ROF pileup (standard) */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(11);

    if (!passRejectITSinROFpileupStrict) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 12 /* Not ITS ROF pileup (strict) */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(12);

    if (!passSelectUPCcollisions) { // if zero then NOT upc, otherwise UPC
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 13 /* is UPC event */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(13);

    if (!passRejectCollInTimeRangeNarrow) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 14 /* Not ITS ROF pileup (strict) */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(14);

    if (collision.multFT0C() < upcRejection.maxFT0CforZNACselection &&
        collision.multZNA() < upcRejection.minZNACsignal &&
        collision.multZNC() < upcRejection.minZNACsignal) {
      return;
    }
    if (collision.multFT0C() < upcRejection.maxFT0CforFV0Aselection &&
        collision.multFV0A() < upcRejection.minFV0Asignal) {
      return;
    }
    if (collision.multFT0C() < upcRejection.maxFT0CforFDDAselection &&
        collision.multFDDA() < upcRejection.minFDDAsignal) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 15 /* pass em/upc rejection */);
    getHist(TH1, histPath + "hCollisionSelection")->Fill(15);

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hNPVContributors"), collision.multNTracksPV());
    histos.fill(HIST("hFT0A_Collisions"), collision.multFT0A() * scaleSignalFT0C);
    histos.fill(HIST("hFT0C_Collisions"), collision.multFT0C() * scaleSignalFT0C);
    histos.fill(HIST("hFT0M_Collisions"), (collision.multFT0A() + collision.multFT0C()) * scaleSignalFT0M);
    histos.fill(HIST("hFV0A_Collisions"), collision.multFV0A() * scaleSignalFV0A);
    histos.fill(HIST("hNGlobalTracks"), collision.multNTracksGlobal());
    histos.fill(HIST("hNMFTTracks"), collision.mftNtracks());
    histos.fill(HIST("hFT0CvsPVz_Collisions_All"), collision.multPVz(), collision.multFT0C() * scaleSignalFT0C);
    histos.fill(HIST("hFV0AvsPVz_Collisions"), collision.multPVz(), collision.multFV0A() * scaleSignalFV0A);
    histos.fill(HIST("hNGlobalTracksvsPVz_Collisions"), collision.multPVz(), collision.multNTracksGlobal());
    histos.fill(HIST("hNMFTTracksvsPVz_Collisions"), collision.multPVz(), collision.mftNtracks());

    // save vertex-Z equalized
    getHist(TH1, histPath + "hNPVContributors")->Fill(multNTracksPV);
    getHist(TH1, histPath + "hFT0A_Collisions")->Fill(multFT0A * scaleSignalFT0A);
    getHist(TH1, histPath + "hFT0C_Collisions")->Fill(multFT0C * scaleSignalFT0C);
    getHist(TH1, histPath + "hFT0M_Collisions")->Fill((multFT0A + multFT0C) * scaleSignalFT0M);
    getHist(TH1, histPath + "hFV0A_Collisions")->Fill(multFV0A * scaleSignalFV0A);
    getHist(TH1, histPath + "hNGlobalTracks")->Fill(multNTracksGlobal);
    getHist(TH1, histPath + "hNMFTTracks")->Fill(mftNtracks);

    if (applyVertexZEqualization.value) {
      // save unequalized for cross-checks
      getHist(TH1, histPath + "hNPVContributors_Unequalized")->Fill(collision.multNTracksPV());
      getHist(TH1, histPath + "hFT0C_Collisions_Unequalized")->Fill(collision.multFT0C() * scaleSignalFT0C);
      getHist(TH1, histPath + "hFT0M_Collisions_Unequalized")->Fill((collision.multFT0A() + collision.multFT0C()) * scaleSignalFT0M);
      getHist(TH1, histPath + "hFV0A_Collisions_Unequalized")->Fill(collision.multFV0A() * scaleSignalFV0A);
      getHist(TH1, histPath + "hNGlobalTracks_Unequalized")->Fill(collision.multNTracksGlobal());
      getHist(TH1, histPath + "hNMFTTracks_Unequalized")->Fill(collision.mftNtracks());
    }

    if (do2DPlots) {
      histos.fill(HIST("hNContribsVsFT0C"), collision.multFT0C() * scaleSignalFT0C, collision.multPVTotalContributors());
      histos.fill(HIST("hNContribsVsFV0A"), collision.multFV0A() * scaleSignalFV0A, collision.multPVTotalContributors());
      histos.fill(HIST("hMatchedVsITSOnly"), collision.multNTracksITSOnly(), collision.multNTracksITSTPC());
      getHist(TH2, histPath + "hNContribsVsFT0C")->Fill(collision.multFT0C() * scaleSignalFT0C, collision.multPVTotalContributors());
      getHist(TH2, histPath + "hNContribsVsFV0A")->Fill(collision.multFV0A() * scaleSignalFV0A, collision.multPVTotalContributors());
      getHist(TH2, histPath + "hMatchedVsITSOnly")->Fill(collision.multNTracksITSOnly(), collision.multNTracksITSTPC());

      // correlate also FIT detector signals
      histos.fill(HIST("hFT0AVsFT0C"), collision.multFT0C() * scaleSignalFT0C, collision.multFT0A());
      histos.fill(HIST("hFV0AVsFT0C"), collision.multFT0C() * scaleSignalFT0C, collision.multFV0A());
      histos.fill(HIST("hFDDAVsFT0C"), collision.multFT0C() * scaleSignalFT0C, collision.multFDDA());
      histos.fill(HIST("hFDDCVsFT0C"), collision.multFT0C() * scaleSignalFT0C, collision.multFDDC());
      getHist(TH2, histPath + "hFT0AVsFT0C")->Fill(collision.multFT0C() * scaleSignalFT0C, collision.multFT0A());
      getHist(TH2, histPath + "hFV0AVsFT0C")->Fill(collision.multFT0C() * scaleSignalFT0C, collision.multFV0A());
      getHist(TH2, histPath + "hFDDAVsFT0C")->Fill(collision.multFT0C() * scaleSignalFT0C, collision.multFDDA());
      getHist(TH2, histPath + "hFDDCVsFT0C")->Fill(collision.multFT0C() * scaleSignalFT0C, collision.multFDDC());
    }

    if (doOccupancyStudyVsCentrality2d) {
      histos.fill(HIST("hNcontribsProfileVsTrackOccupancyVsFT0C"), collision.trackOccupancyInTimeRange(), collision.multFT0C(), collision.multPVTotalContributors());
      histos.fill(HIST("hNGlobalTracksProfileVsTrackOccupancyVsFT0C"), collision.trackOccupancyInTimeRange(), collision.multFT0C(), collision.multNTracksGlobal());
      histos.fill(HIST("hNcontribsProfileVsFT0COccupancyVsFT0C"), collision.ft0cOccupancyInTimeRange(), collision.multFT0C(), collision.multPVTotalContributors());
      histos.fill(HIST("hNGlobalTracksProfileVsFT0COccupancyVsFT0C"), collision.ft0cOccupancyInTimeRange(), collision.multFT0C(), collision.multNTracksGlobal());
    }

    if (doOccupancyStudyVsRawValues3d) {
      histos.fill(HIST("hTrackOccupancyVsNContribsVsFT0C"), collision.trackOccupancyInTimeRange(), collision.multPVTotalContributors(), collision.multFT0C());
      histos.fill(HIST("hTrackOccupancyVsNGlobalTracksVsFT0C"), collision.trackOccupancyInTimeRange(), collision.multNTracksGlobal(), collision.multFT0C());
      histos.fill(HIST("hFT0COccupancyVsNContribsVsFT0C"), collision.ft0cOccupancyInTimeRange(), collision.multPVTotalContributors(), collision.multFT0C());
      histos.fill(HIST("hFT0COccupancyVsNGlobalTracksVsFT0C"), collision.ft0cOccupancyInTimeRange(), collision.multNTracksGlobal(), collision.multFT0C());
    }

    if (doNGlobalTracksVsRawSignals) {
      histos.fill(HIST("hNGlobalTracksVsFT0A"), multFT0A, multNTracksGlobal);
      histos.fill(HIST("hNGlobalTracksVsFT0C"), multFT0C, multNTracksGlobal);
      histos.fill(HIST("hNGlobalTracksVsFT0M"), (multFT0A + multFT0C), multNTracksGlobal);
      histos.fill(HIST("hNGlobalTracksVsFV0A"), multFV0A, multNTracksGlobal);
      histos.fill(HIST("hNGlobalTracksVsNMFTTracks"), mftNtracks, multNTracksGlobal);
      histos.fill(HIST("hNGlobalTracksVsNTPV"), multNTracksPV, multNTracksGlobal);

      // per run
      getHist(TH2, histPath + "hNGlobalTracksVsFT0A")->Fill(multFT0A, multNTracksGlobal);
      getHist(TH2, histPath + "hNGlobalTracksVsFT0C")->Fill(multFT0C, multNTracksGlobal);
      getHist(TH2, histPath + "hNGlobalTracksVsFT0M")->Fill(multFT0A + multFT0C, multNTracksGlobal);
      getHist(TH2, histPath + "hNGlobalTracksVsFV0A")->Fill(multFV0A, multNTracksGlobal);
      getHist(TH2, histPath + "hNGlobalTracksVsNMFTTracks")->Fill(mftNtracks, multNTracksGlobal);
      getHist(TH2, histPath + "hNGlobalTracksVsNTPV")->Fill(multNTracksPV, multNTracksGlobal);
    }

    if constexpr (requires { collision.multMCExtraId(); }) {
      // requires monte carlo information
      if (collision.multMCExtraId() > -1) {
        auto mcCollision = collision.template multMCExtra_as<soa::Join<aod::MultMCExtras, aod::MultHepMCHIs>>();
        histos.fill(HIST("hImpactParameterVsFT0A"), multFT0A, mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsFT0C"), multFT0C, mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsFT0M"), (multFT0A + multFT0C), mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsFV0A"), multFV0A, mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsNMFTTracks"), mftNtracks, mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsNTPV"), multNTracksPV, mcCollision.impactParameter());

        histos.fill(HIST("hImpactParameterVsMCFT0A"), mcCollision.multMCFT0A(), mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsMCFT0C"), mcCollision.multMCFT0C(), mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsMCFT0M"), (mcCollision.multMCFT0A() + mcCollision.multMCFT0C()), mcCollision.impactParameter());
        histos.fill(HIST("hImpactParameterVsMCFV0A"), mcCollision.multMCFV0A(), mcCollision.impactParameter());
      }
    }

    // if the table has centrality information
    if constexpr (requires { collision.centFT0C(); }) {
      // process FT0C centrality plots
      histos.fill(HIST("hCentrality"), collision.centFT0C());
      histos.fill(HIST("hNContribsVsCentrality"), collision.centFT0C(), collision.multPVTotalContributors());
      histos.fill(HIST("hNITSTPCTracksVsCentrality"), collision.centFT0C(), collision.multNTracksITSTPC());
      histos.fill(HIST("hNITSOnlyTracksVsCentrality"), collision.centFT0C(), collision.multNTracksITSOnly());
      histos.fill(HIST("hNGlobalTracksVsCentrality"), collision.centFT0C(), collision.multNTracksGlobal());
      histos.fill(HIST("hNMFTTracksVsCentrality"), collision.centFT0C(), collision.mftNtracks());
      histos.fill(HIST("hPVChi2VsCentrality"), collision.centFT0C(), collision.multPVChi2());
      getHist(TH1, histPath + "hCentrality")->Fill(collision.centFT0C());
      getHist(TH2, histPath + "hNContribsVsCentrality")->Fill(collision.centFT0C(), collision.multPVTotalContributors());
      getHist(TH2, histPath + "hNITSTPCTracksVsCentrality")->Fill(collision.centFT0C(), collision.multNTracksITSTPC());
      getHist(TH2, histPath + "hNITSOnlyTracksVsCentrality")->Fill(collision.centFT0C(), collision.multNTracksITSOnly());
      getHist(TH2, histPath + "hNGlobalTracksVsCentrality")->Fill(collision.centFT0C(), collision.multNTracksGlobal());
      getHist(TH2, histPath + "hNMFTTracksVsCentrality")->Fill(collision.centFT0C(), collision.mftNtracks());
      getHist(TH2, histPath + "hPVChi2VsCentrality")->Fill(collision.centFT0C(), collision.multPVChi2());

      if (doOccupancyStudyVsCentrality2d) {
        histos.fill(HIST("hNcontribsProfileVsTrackOccupancyVsCentrality"), collision.trackOccupancyInTimeRange(), collision.centFT0C(), collision.multPVTotalContributors());
        histos.fill(HIST("hNGlobalTracksProfileVsTrackOccupancyVsCentrality"), collision.trackOccupancyInTimeRange(), collision.centFT0C(), collision.multNTracksGlobal());
        histos.fill(HIST("hNcontribsProfileVsFT0COccupancyVsCentrality"), collision.ft0cOccupancyInTimeRange(), collision.centFT0C(), collision.multPVTotalContributors());
        histos.fill(HIST("hNGlobalTracksProfileVsFT0COccupancyVsCentrality"), collision.ft0cOccupancyInTimeRange(), collision.centFT0C(), collision.multNTracksGlobal());
      }

      if (doOccupancyStudyVsCentrality3d) {
        histos.fill(HIST("hTrackOccupancyVsNContribsVsCentrality"), collision.trackOccupancyInTimeRange(), collision.multPVTotalContributors(), collision.centFT0C());
        histos.fill(HIST("hTrackOccupancyVsNGlobalTracksVsCentrality"), collision.trackOccupancyInTimeRange(), collision.multNTracksGlobal(), collision.centFT0C());
        histos.fill(HIST("hFT0COccupancyVsNContribsVsCentrality"), collision.ft0cOccupancyInTimeRange(), collision.multPVTotalContributors(), collision.centFT0C());
        histos.fill(HIST("hFT0COccupancyVsNGlobalTracksVsCentrality"), collision.ft0cOccupancyInTimeRange(), collision.multNTracksGlobal(), collision.centFT0C());
      }
    }

    if constexpr (requires { collision.has_multBC(); }) {
      if (doTimeStudies && collision.has_multBC()) {
        initRun(collision);
        auto multbc = collision.template multBC_as<aod::MultBCs>();
        uint64_t bcTimestamp = multbc.timestamp();
        float hoursAfterStartOfRun = static_cast<float>(bcTimestamp - startOfRunTimestamp) / 3600000.0;

        getHist(TH2, histPath + "hFT0AVsTime")->Fill(hoursAfterStartOfRun, collision.multFT0A());
        getHist(TH2, histPath + "hFT0CVsTime")->Fill(hoursAfterStartOfRun, collision.multFT0C());
        getHist(TH2, histPath + "hFT0MVsTime")->Fill(hoursAfterStartOfRun, collision.multFT0M());
        getHist(TH2, histPath + "hFV0AVsTime")->Fill(hoursAfterStartOfRun, collision.multFV0A());
        getHist(TH2, histPath + "hFV0AOuterVsTime")->Fill(hoursAfterStartOfRun, collision.multFV0AOuter());
        getHist(TH2, histPath + "hMFTTracksVsTime")->Fill(hoursAfterStartOfRun, collision.mftNtracks());
        getHist(TH2, histPath + "hNGlobalVsTime")->Fill(hoursAfterStartOfRun, collision.multNTracksGlobal());
        getHist(TH2, histPath + "hNTPVContributorsVsTime")->Fill(hoursAfterStartOfRun, collision.multPVTotalContributors());
        getHist(TProfile, histPath + "hPVzProfileCoVsTime")->Fill(hoursAfterStartOfRun, collision.multPVz());
        getHist(TProfile, histPath + "hPVzProfileBcVsTime")->Fill(hoursAfterStartOfRun, multbc.multFT0PosZ());
        if (doTimeStudyFV0AOuterVsFT0A3d) {
          histos.fill(HIST("h3dFV0AVsTime"), hoursAfterStartOfRun, collision.multFV0A(), collision.multFV0AOuter());
        }

        if (irDoRateVsTime) {
          float interactionRate = mRateFetcher.fetch(ccdb.service, bcTimestamp, mRunNumber, irSource.value, irCrashOnNull) / 1000.; // kHz
          getHist(TProfile, histPath + "hIRProfileVsTime")->Fill(hoursAfterStartOfRun, interactionRate);
        }
      }
    }
  }

  void processCollisions(soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultsGlobal, aod::MultSelections, aod::Mults2BC, aod::FV0AOuterMults>::iterator const& collision, aod::MultBCs const&)
  {
    genericProcessCollision(collision);
  }

  void processCollisionsWithResolutionStudy(soa::Join<aod::MultsRun3, aod::MFTMults, aod::Mult2MCExtras, aod::MultsExtra, aod::MultsGlobal, aod::MultSelections, aod::Mults2BC, aod::FV0AOuterMults>::iterator const& collision, soa::Join<aod::MultMCExtras, aod::MultHepMCHIs> const&)
  {
    genericProcessCollision(collision);
  }

  void processCollisionsWithCentrality(soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal, aod::Mults2BC, aod::FV0AOuterMults>::iterator const& collision, aod::MultBCs const&)
  {
    genericProcessCollision(collision);
  }

  void processCollisionsWithCentralityWithNeighbours(soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal, aod::MultNeighs, aod::FV0AOuterMults>::iterator const& collision)
  {
    genericProcessCollision(collision);
  }

  void processBCs(soa::Join<aod::BC2Mults, aod::MultBCs>::iterator const& multbc, soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal> const&)
  {
    // process BCs, calculate FT0C distribution
    // conditionals suggested by FIT team (Jacek O. et al)
    histos.fill(HIST("hBCSelection"), 0); // all BCs
    if (selectCollidingBCs && !multbc.multCollidingBC())
      return;
    histos.fill(HIST("hBCSelection"), 1); // colliding
    if (selectTVX && !multbc.multTVX())
      return;
    histos.fill(HIST("hBCSelection"), 2); // TVX
    if (selectFV0OrA && !multbc.multFV0OrA())
      return;
    histos.fill(HIST("hBCSelection"), 3); // FV0OrA
    if (vertexZwithT0 < 100.0f) {
      if (!multbc.multFT0PosZValid())
        return;
      if (TMath::Abs(multbc.multFT0PosZ()) > vertexZwithT0)
        return;
    }
    histos.fill(HIST("hBCSelection"), 4); // FV0OrA

    if (multbc.multFT0C() < upcRejection.maxFT0CforZNACselection &&
        multbc.multZNA() < upcRejection.minZNACsignal &&
        multbc.multZNC() < upcRejection.minZNACsignal) {
      return;
    }
    if (multbc.multFT0C() < upcRejection.maxFT0CforFV0Aselection &&
        multbc.multFV0A() < upcRejection.minFV0Asignal) {
      return;
    }
    if (multbc.multFT0C() < upcRejection.maxFT0CforFDDAselection &&
        multbc.multFDDA() < upcRejection.minFDDAsignal) {
      return;
    }

    histos.fill(HIST("hBCSelection"), 5); // znac

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hFT0C_BCs"), multbc.multFT0C() * scaleSignalFT0C);

    // ZN signals
    histos.fill(HIST("hZNAvsFT0C_BCs"), multbc.multFT0C() * scaleSignalFT0C, multbc.multZNA());
    histos.fill(HIST("hZNCvsFT0C_BCs"), multbc.multFT0C() * scaleSignalFT0C, multbc.multZNC());

    histos.fill(HIST("hFT0M_BCs"), (multbc.multFT0A() + multbc.multFT0C()) * scaleSignalFT0M);
    histos.fill(HIST("hFV0A_BCs"), multbc.multFV0A() * scaleSignalFV0A);
    if (multbc.multFT0PosZValid()) {
      histos.fill(HIST("hFT0CvsPVz_BCs_All"), multbc.multFT0PosZ(), multbc.multFT0C() * scaleSignalFT0C);
      if (multbc.multFT0C() > minFT0CforVertexZ) {
        histos.fill(HIST("hFT0CvsPVz_BCs"), multbc.multFT0PosZ(), multbc.multFT0C() * scaleSignalFT0C);
      }
    }

    if (multbc.has_ft0Mult()) {
      auto multco = multbc.ft0Mult_as<soa::Join<aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal>>();
      if (multbc.multFT0PosZValid()) {
        histos.fill(HIST("hVertexZ_BCvsCO"), multco.multPVz(), multbc.multFT0PosZ());
      }
    }
  }

  PROCESS_SWITCH(centralityStudy, processCollisions, "per-collision analysis", false);
  PROCESS_SWITCH(centralityStudy, processCollisionsWithResolutionStudy, "per-collision analysis, with reso study", false);
  PROCESS_SWITCH(centralityStudy, processCollisionsWithCentrality, "per-collision analysis", true);
  PROCESS_SWITCH(centralityStudy, processCollisionsWithCentralityWithNeighbours, "per-collision analysis", false);
  PROCESS_SWITCH(centralityStudy, processBCs, "per-BC analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<centralityStudy>(cfgc)};
}
