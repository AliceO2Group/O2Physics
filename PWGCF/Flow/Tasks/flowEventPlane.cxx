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

/// \file flowEventPlane.cxx
/// \brief Flow calculation using event plane.
/// \author Yash Patley <yash.patley@cern.ch>

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <array>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

enum GainClibCorr {
  kGainCalibA = 0,
  kGainCalibC,
  kNGainCalib
};

enum CorrectionType {
  kFineCorr = 0,
  kCoarseCorr,
  kNCorr
};

enum CollisionParameterType {
  kCent = 0,
  kVx,
  kVy,
  kVz
};

enum ZDCXYType {
  kXa = 0,
  kYa,
  kXc,
  kYc,
  kXYAC
};

enum ParticleType {
  kPi = 0,
  kKa,
  kPr,
  kNPart
};

struct FlowEventPlane {
  // Configurables
  // Collisions
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinCent{"cMinCent", 0., "Minumum Centrality"};
  Configurable<float> cMaxCent{"cMaxCent", 100.0, "Maximum Centrality"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<float> cMinOccupancy{"cMinOccupancy", 0, "Minimum FT0C Occupancy"};
  Configurable<float> cMaxOccupancy{"cMaxOccupancy", 1e6, "Maximum FT0C Occupancy"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.1, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 10.0, "p_{T} maximum"};
  Configurable<int> cNEtaBins{"cNEtaBins", 7, "# of eta bins"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<bool> cTrackGlobal{"cTrackGlobal", true, "Global Track"};
  Configurable<float> cTrackDcaXYCut{"cTrackDcaXYCut", 0.1, "DcaXY Cut"};
  Configurable<float> cTrackDcaZCut{"cTrackDcaZCut", 1., "DcaXY Cut"};
  Configurable<int> cNRapBins{"cNRapBins", 5, "# of y bins"};
  Configurable<int> cNInvMassBins{"cNInvMassBins", 500, "# of m bins"};

  // Track PID
  Configurable<float> cTpcNSigmaCut{"cTpcNSigmaCut", 2, "TPC NSigma Cut"};
  Configurable<float> cTpcRejCut{"cTpcRejCut", 3, "TPC Rej Cut"};
  Configurable<float> cTofNSigmaCut{"cTofNSigmaCut", 2, "TOF NSigma Cut"};
  Configurable<float> cTofRejCut{"cTofRejCut", 3, "TOF Rej Cut"};
  Configurable<float> cPionPtCut{"cPionPtCut", 0.6, "Pion TPC pT cutoff"};
  Configurable<float> cKaonPtCut{"cKaonPtCut", 0.6, "Kaon TPC pT cutoff"};
  Configurable<float> cProtonPtCut{"cProtonPtCut", 1.1, "Proton TPC pT cutoff"};

  // Resonance
  Configurable<float> cResRapCut{"cResRapCut", 0.5, "Resonance rapidity cut"};

  // Gain calibration
  Configurable<bool> cDoGainCalib{"cDoGainCalib", false, "Gain Calib Flag"};
  Configurable<bool> cUseAlphaZDC{"cUseAlphaZDC", true, "Use Alpha ZDC"};

  // Coarse binning factor
  Configurable<int> cAxisCBF{"cAxisCBF", 5, "Coarse Bin Factor"};

  // Cent Vx Vy Vz Bins
  Configurable<int> cAxisCentBins{"cAxisCentBins", 20, "NBins Centrality"};
  Configurable<int> cAxisVxyBins{"cAxisVxyBins", 20, "NBins Vx Vy"};
  Configurable<int> cAxisVzBins{"cAxisVzBins", 20, "NBins Vz"};
  Configurable<float> cAxisVxMin{"cAxisVxMin", -0.06, "Vx Min"};
  Configurable<float> cAxisVxMax{"cAxisVxMax", -0.02, "Vx Max"};
  Configurable<float> cAxisVyMin{"cAxisVyMin", -0.01, "Vy Min"};
  Configurable<float> cAxisVyMax{"cAxisVyMax", 0.006, "Vy Max"};

  // Corrections
  Configurable<bool> cApplyRecentCorr{"cApplyRecentCorr", false, "Apply recentering"};
  Configurable<std::vector<int>> cCorrFlagVector{"cCorrFlagVector", {0, 0, 0, 0, 0, 0}, "Correction Flag"};

  // CCDB
  Configurable<std::string> cCcdbUrl{"cCcdbUrl", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cCcdbPath{"cCcdbPath", "Users/y/ypatley/DFOO", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdbService;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Global objects
  const float zdcDenThrs = 1e-4;
  float cent = 0., mult = 0.;
  float posX = 0., posY = 0., posZ = 0.;
  std::vector<std::vector<std::string>> vCoarseCorrHistNames = {
    {"hXZNAVsCentVxVyVz"},
    {"hYZNAVsCentVxVyVz"},
    {"hXZNCVsCentVxVyVz"},
    {"hYZNCVsCentVxVyVz"}};
  std::vector<std::vector<std::string>> vFineCorrHistNames = {
    {"hXZNAVsCent", "hXZNAVsVx", "hXZNAVsVy", "hXZNAVsVz"},
    {"hYZNAVsCent", "hYZNAVsVx", "hYZNAVsVy", "hYZNAVsVz"},
    {"hXZNCVsCent", "hXZNCVsVx", "hXZNCVsVy", "hXZNCVsVz"},
    {"hYZNCVsCent", "hYZNCVsVx", "hYZNCVsVy", "hYZNCVsVz"}};
  std::map<CorrectionType, std::vector<std::vector<std::string>>> corrTypeHistNameMap = {{kFineCorr, vFineCorrHistNames}, {kCoarseCorr, vCoarseCorrHistNames}};

  // Container for histograms
  struct CorrectionHistContainer {
    std::array<TH2F*, 2> hGainCalib;
    std::array<std::array<std::array<THnSparseF*, 1>, 4>, 6> vCoarseCorrHist;
    std::array<std::array<std::array<TProfile*, 4>, 4>, 6> vFineCorrHist;
  } CorrectionHistContainer;

  // Run number
  int cRunNum = 0, lRunNum = 0;

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdbService->setURL(cCcdbUrl.value);
    ccdbService->setCaching(true);

    // Define axes
    const AxisSpec axisZDCEnergy{500, 0, 500, "ZD[AC] Signal"};

    const AxisSpec axisCent{100, 0., 100, "FT0C%"};
    const AxisSpec axisVx{cAxisVxyBins, cAxisVxMin, cAxisVxMax, "V_{X}(cm)"};
    const AxisSpec axisVy{cAxisVxyBins, cAxisVyMin, cAxisVyMax, "V_{Y}(cm)"};
    const AxisSpec axisVz{cAxisVzBins, cMinZVtx, cMaxZVtx, "V_{Z}(cm)"};

    const AxisSpec axisCoarseCent{cAxisCentBins / cAxisCBF, cMinCent, cMaxCent, "FT0C%"};
    const AxisSpec axisCoarseVx{cAxisVxyBins / cAxisCBF, cAxisVxMin, cAxisVxMax, "V_{x}"};
    const AxisSpec axisCoarseVy{cAxisVxyBins / cAxisCBF, cAxisVyMin, cAxisVyMax, "V_{y}"};
    const AxisSpec axisCoarseVz{cAxisVzBins / cAxisCBF, cMinZVtx, cMaxZVtx, "V_{z}"};

    const AxisSpec axisFineCent{cAxisCentBins, cMinCent, cMaxCent, "FT0C%"};
    const AxisSpec axisFineVx{cAxisVxyBins, cAxisVxMin, cAxisVxMax, "V_{x}"};
    const AxisSpec axisFineVy{cAxisVxyBins, cAxisVyMin, cAxisVyMax, "V_{x}"};
    const AxisSpec axisFineVz{cAxisVzBins, cMinZVtx, cMaxZVtx, "V_{z}"};

    const AxisSpec axisXa{40, -1, 1, "X^{ZNA}_{1}"};
    const AxisSpec axisYa{40, -1, 1, "Y^{ZNA}_{1}"};
    const AxisSpec axisXc{40, -1, 1, "X^{ZNC}_{1}"};
    const AxisSpec axisYc{40, -1, 1, "Y^{ZNC}_{1}"};

    const AxisSpec axisPsi{18, -PIHalf, PIHalf, "#Psi_{SP}"};

    const AxisSpec axisXYac{600, -6, 6, "Q^{t}Q^{p}"};
    const AxisSpec axisV1{400, -4, 4, "v_{1}"};

    const AxisSpec axisTrackPt{100, 0., 10., "p_{T} (GeV/#it{c})"};
    const AxisSpec axisTrackEta{cNEtaBins, -0.8, 0.8, "#eta"};
    const AxisSpec axisTrackRap{cNRapBins, -0.5, 0.5, "y"};
    const AxisSpec axisInvMass{cNInvMassBins, 0.87, 1.12, "M_{KK} (GeV/#it{c}^{2}"};

    const AxisSpec axisTrackDcaXY{60, -0.15, 0.15, "DCA_{XY}"};
    const AxisSpec axisTrackDcaZ{230, -1.15, 1.15, "DCA_{XY}"};
    const AxisSpec axisTrackdEdx{360, 20, 200, "#frac{dE}{dx}"};
    const AxisSpec axisTrackNSigma{161, -4.025, 4.025, {"n#sigma"}};

    // Create histograms
    // Event
    histos.add("Event/hCent", "FT0C%", kTH1F, {axisCent});
    histos.add("Event/hVx", "V_{x}", kTH1F, {axisVx});
    histos.add("Event/hVy", "V_{y}", kTH1F, {axisVy});
    histos.add("Event/hVz", "V_{z}", kTH1F, {axisVz});

    // Gain calib
    histos.add("QA/GainCalib/hZNASignal", "ZNA Signal", kTH2F, {{4, 0, 4}, {axisZDCEnergy}});
    histos.add("QA/GainCalib/hZNCSignal", "ZNC Signal", kTH2F, {{4, 0, 4}, {axisZDCEnergy}});
    histos.add("QA/hZNASignal", "ZNA Signal", kTProfile2D, {{4, 0, 4}, {axisVz}});
    histos.add("QA/hZNCSignal", "ZNC Signal", kTProfile2D, {{4, 0, 4}, {axisVz}});
    histos.add("QA/hZNAEnergyCommon", "ZNA Energy Common", kTProfile, {axisVz});
    histos.add("QA/hZNCEnergyCommon", "ZNC Energy Common", kTProfile, {axisVz});

    // Corrections
    histos.add("CorrHist/hWtXZNA", "X^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtYZNA", "Y^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtXZNC", "X^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtYZNC", "Y^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtXZNA", "X^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtYZNA", "Y^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtXZNC", "X^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtYZNC", "Y^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hXZNAVsCent", "X^{ZNA}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hXZNAVsVx", "X^{ZNA}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hXZNAVsVy", "X^{ZNA}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hXZNAVsVz", "X^{ZNA}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hYZNAVsCent", "Y^{ZNA}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hYZNAVsVx", "Y^{ZNA}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hYZNAVsVy", "Y^{ZNA}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hYZNAVsVz", "Y^{ZNA}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hXZNCVsCent", "X^{ZNC}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hXZNCVsVx", "X^{ZNC}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hXZNCVsVy", "X^{ZNC}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hXZNCVsVz", "X^{ZNC}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hYZNCVsCent", "Y^{ZNC}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hYZNCVsVx", "Y^{ZNC}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hYZNCVsVy", "Y^{ZNC}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hYZNCVsVz", "Y^{ZNC}_{1} Vs V_{z}", kTProfile, {axisFineVz});

    // Checks
    histos.add("Checks/hPsiSPA", "#Psi_{SP}^{A} distribution", kTH2F, {axisCent, axisPsi});
    histos.add("Checks/hPsiSPC", "#Psi_{SP}^{C} distribution", kTH2F, {axisCent, axisPsi});
    histos.add("Checks/hCosPsiSPAC", "Cos(#Psi_{SP}^{A} #minus #Psi_{SP}^{C}) distribution", kTProfile, {axisCent});
    histos.add("Checks/hSinPsiSPAC", "Sin(#Psi_{SP}^{A} #minus #Psi_{SP}^{C}) distribution", kTProfile, {axisCent});
    histos.add("Checks/hXaXc", "X^{A}_{1}X^{C}_{1}", kTProfile, {axisCent});
    histos.add("Checks/hYaYc", "Y^{A}_{1}Y^{C}_{1}", kTProfile, {axisCent});
    histos.add("Checks/hXaYc", "X^{A}_{1}Y^{C}_{1}", kTProfile, {axisCent});
    histos.add("Checks/hYaXc", "Y^{A}_{1}X^{C}_{1}", kTProfile, {axisCent});

    // Track QA
    histos.add("TrackQA/hPtDcaXY", "DCA_{XY} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaXY});
    histos.add("TrackQA/hPtDcaZ", "DCA_{Z} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaZ});
    histos.add("TrackQA/hTrackTPCdEdX", "hTrackTPCdEdX", kTH2F, {axisTrackPt, axisTrackdEdx});

    // Charged particle directed flow
    histos.add("DF/hQaQc", "X^{A}_{1}X^{C}_{1} + Y^{A}_{1}Y^{C}_{1}", kTProfile, {axisCent});
    histos.add("DF/hAQu", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("DF/hCQu", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("DF/hAQuPos", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("DF/hCQuPos", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("DF/hAQuNeg", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("DF/hCQuNeg", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTProfile2D, {axisCent, axisTrackEta});

    // Identified particle
    histos.add("PartId/Pion/hdEdX", "PartId/Pion/hdEdX", kTH2F, {axisTrackPt, axisTrackdEdx});
    histos.add("PartId/Pion/hTPCNSigma", "PartId/Pion/hTPCNSigma", kTH2F, {axisTrackPt, axisTrackNSigma});
    histos.add("PartId/Pion/hTOFNSigma", "PartId/Pion/hTOFNSigma", kTH2F, {axisTrackPt, axisTrackNSigma});
    histos.add("PartId/Pion/hAQuPos", "PartId/Pion/hAQuPos", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("PartId/Pion/hAQuNeg", "PartId/Pion/hAQuNeg", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("PartId/Pion/hCQuPos", "PartId/Pion/hCQuPos", kTProfile2D, {axisCent, axisTrackEta});
    histos.add("PartId/Pion/hCQuNeg", "PartId/Pion/hCQuNeg", kTProfile2D, {axisCent, axisTrackEta});
    histos.addClone("PartId/Pion/", "PartId/Kaon/");
    histos.addClone("PartId/Pion/", "PartId/Proton/");

    // Resonance
    histos.add("Reso/Phi/hSigCentPtInvMass", "hUSCentPtInvMass", kTH3F, {axisCent, axisTrackPt, axisInvMass});
    histos.add("Reso/Phi/hBkgCentPtInvMass", "hLSCentPtInvMass", kTH3F, {axisCent, axisTrackPt, axisInvMass});
    histos.add("Reso/Phi/Sig/hPhiQuA", "hPhiQuA", kTProfile3D, {axisCent, axisTrackRap, axisInvMass});
    histos.add("Reso/Phi/Sig/hPhiQuC", "hPhiQuC", kTProfile3D, {axisCent, axisTrackRap, axisInvMass});
    histos.add("Reso/Phi/Bkg/hPhiQuA", "hPhiQuA", kTProfile3D, {axisCent, axisTrackRap, axisInvMass});
    histos.add("Reso/Phi/Bkg/hPhiQuC", "hPhiQuC", kTProfile3D, {axisCent, axisTrackRap, axisInvMass});
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) { // VtxZ selection
      return false;
    }

    if (cSel8Trig && !col.sel8()) { // Sel8 selection
      return false;
    }

    cent = col.centFT0C();
    if (cent <= cMinCent || cent >= cMaxCent) { // Centrality selection
      return false;
    }

    if (col.ft0cOccupancyInTimeRange() < cMinOccupancy || col.ft0cOccupancyInTimeRange() > cMaxOccupancy) { // Occupancy cut
      return false;
    }

    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX)) { // Time and Vertex trigger
      return false;
    }

    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) { // Time frame border
      return false;
    }

    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) { // ITS Readout frame border
      return false;
    }

    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) { // ITS+TPC Vertex
      return false;
    }

    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) { // Pile-up rejection
      return false;
    }

    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) { // ZvtxFT0 vs PV
      return false;
    }

    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) { // All ITS layer active
      return false;
    }

    // Set Multiplicity
    mult = col.multNTracksHasTPC();

    return true;
  }

  // Track Selection
  template <typename T>
  bool selectTrack(T const& track)
  {
    if (track.pt() <= cTrackMinPt || track.pt() >= cTrackMaxPt || std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    if (cTrackGlobal && !track.isGlobalTrackWoDCA()) {
      return false;
    }

    if (std::abs(track.dcaXY()) > cTrackDcaXYCut || std::abs(track.dcaZ()) > cTrackDcaZCut) {
      return false;
    }

    return true;
  }

  template <ParticleType part1, ParticleType part2, ParticleType part3>
  bool checkTrackPid(float const& ptCut, float const& trackPt, std::vector<float> const& vTpcNsig, std::vector<float> const& vTofNsig, bool const& tofFlag)
  {
    bool retFlag = false;
    if (tofFlag) {
      if (vTofNsig[part1] < cTofNSigmaCut && vTofNsig[part2] > cTofRejCut && vTofNsig[part3] > cTofRejCut && vTpcNsig[part1] < cTpcNSigmaCut) {
        retFlag = true;
      }
    } else {
      if (trackPt < ptCut && vTpcNsig[part1] < cTpcNSigmaCut && vTpcNsig[part2] > cTpcRejCut && vTpcNsig[part3] > cTpcRejCut) {
        retFlag = true;
      }
    }
    return retFlag;
  }

  template <ParticleType partType, typename T>
  bool identifyTrack(T const& track)
  {
    std::vector<float> vPtCut = {cPionPtCut, cKaonPtCut, cProtonPtCut};
    std::vector<float> vTpcNsig = {std::abs(track.tpcNSigmaPi()), std::abs(track.tpcNSigmaKa()), std::abs(track.tpcNSigmaPr())};
    std::vector<float> vTofNsig = {std::abs(track.tofNSigmaPi()), std::abs(track.tofNSigmaKa()), std::abs(track.tofNSigmaPr())};
    bool retFlag = false;

    if (partType == kPi && checkTrackPid<kPi, kKa, kPr>(vPtCut[kPi], track.pt(), vTpcNsig, vTofNsig, track.hasTOF())) {
      retFlag = true;
    } else if (partType == kKa && checkTrackPid<kKa, kPi, kPr>(vPtCut[kKa], track.pt(), vTpcNsig, vTofNsig, track.hasTOF())) {
      retFlag = true;
    } else if (partType == kPr && checkTrackPid<kPr, kPi, kKa>(vPtCut[kPr], track.pt(), vTpcNsig, vTofNsig, track.hasTOF())) {
      retFlag = true;
    } else {
      return false;
    }

    return retFlag;
  }

  void gainCalib(bool const& loadGainCalib, float const& vz, std::array<float, 4>& eA, std::array<float, 4>& eC)
  {
    // Store gain calibration histograms per run number
    if (loadGainCalib) {
      std::string ccdbPath = static_cast<std::string>(cCcdbPath) + "/GainCalib" + "/Run" + std::to_string(cRunNum);
      auto ccdbObj = ccdbService->getForTimeStamp<TList>(ccdbPath, -1);
      CorrectionHistContainer.hGainCalib[0] = reinterpret_cast<TH2F*>(ccdbObj->FindObject("hZNASignal"));
      CorrectionHistContainer.hGainCalib[1] = reinterpret_cast<TH2F*>(ccdbObj->FindObject("hZNCSignal"));
    }

    // Apply gain calibration
    float vA = 0., vC = 0.;
    for (int i = 0; i < static_cast<int>(eA.size()); ++i) {
      vA = CorrectionHistContainer.hGainCalib[0]->GetBinContent(CorrectionHistContainer.hGainCalib[0]->FindBin(i + 0.5, vz + 0.00001));
      vC = CorrectionHistContainer.hGainCalib[1]->GetBinContent(CorrectionHistContainer.hGainCalib[1]->FindBin(i + 0.5, vz + 0.00001));
      eA[i] *= vA;
      eC[i] *= vC;
    }
  }

  std::vector<float> getAvgCorrFactors(int const& itr, CorrectionType const& corrType, std::array<float, 4> const& vCollParam)
  {
    std::vector<float> vAvgOutput = {0., 0., 0., 0.};
    int binarray[4];
    if (corrType == kCoarseCorr) {
      int cntrx = 0;
      for (auto const& v : CorrectionHistContainer.vCoarseCorrHist[itr]) {
        for (auto const& h : v) {
          binarray[kCent] = h->GetAxis(kCent)->FindBin(vCollParam[kCent] + 0.0001);
          binarray[kVx] = h->GetAxis(kVx)->FindBin(vCollParam[kVx] + 0.0001);
          binarray[kVy] = h->GetAxis(kVy)->FindBin(vCollParam[kVy] + 0.0001);
          binarray[kVz] = h->GetAxis(kVz)->FindBin(vCollParam[kVz] + 0.0001);
          vAvgOutput[cntrx] += h->GetBinContent(h->GetBin(binarray));
        }
        ++cntrx;
      }
    } else {
      int cntrx = 0;
      for (auto const& v : CorrectionHistContainer.vFineCorrHist[itr]) {
        int cntry = 0;
        for (auto const& h : v) {
          vAvgOutput[cntrx] += h->GetBinContent(h->GetXaxis()->FindBin(vCollParam[cntry] + 0.0001));
          ++cntry;
        }
        ++cntrx;
      }
    }

    return vAvgOutput;
  }

  void applyCorrection(bool const& loadShiftCorr, std::array<float, 4> const& inputParam, std::array<float, 4>& outputParam)
  {
    std::vector<int> vCorrFlags = static_cast<std::vector<int>>(cCorrFlagVector);
    int nitr = vCorrFlags.size();
    CorrectionType corrType = kFineCorr;
    std::string ccdbPath;

    // Correction iterations
    for (int i = 0; i < nitr; ++i) {
      // Don't correct if corrFlag != 1
      if (vCorrFlags[i] != 1) {
        continue;
      }

      // Set correction type
      if (i % kNCorr == 0) {
        corrType = kCoarseCorr;
      } else {
        corrType = kFineCorr;
      }

      // Check current and last run number, fetch ccdb object and store corrections in container
      if (loadShiftCorr) {
        // Set ccdb path
        ccdbPath = static_cast<std::string>(cCcdbPath) + "/CorrItr_" + std::to_string(i + 1) + "/Run" + std::to_string(cRunNum);

        // Get object from CCDB
        auto ccdbObject = ccdbService->getForTimeStamp<TList>(ccdbPath, -1);

        // Check CCDB Object
        if (!ccdbObject) {
          LOGF(warning, "CCDB OBJECT NOT FOUND");
          return;
        }

        // Store histograms in Hist Container
        std::vector<std::vector<std::string>> vHistNames = corrTypeHistNameMap.at(corrType);
        int cntrx = 0;
        for (auto const& x : vHistNames) {
          int cntry = 0;
          for (auto const& y : x) {
            if (corrType == kFineCorr) {
              CorrectionHistContainer.vFineCorrHist[i][cntrx][cntry] = reinterpret_cast<TProfile*>(ccdbObject->FindObject(y.c_str()));
            } else {
              CorrectionHistContainer.vCoarseCorrHist[i][cntrx][cntry] = reinterpret_cast<THnSparseF*>(ccdbObject->FindObject(y.c_str()));
            }
            ++cntry;
          }
          ++cntrx;
        }
      }

      // Get averages
      std::vector<float> vAvg = getAvgCorrFactors(i, corrType, inputParam);

      // Apply correction
      outputParam[kXa] -= vAvg[kXa];
      outputParam[kYa] -= vAvg[kYa];
      outputParam[kXc] -= vAvg[kXc];
      outputParam[kYc] -= vAvg[kYc];
    }
  }

  template <typename T>
  void getResoFlow(T const& tracks, std::array<float, 4> const& vSP)
  {
    float ux = 0., uy = 0., v1a = 0., v1c = 0.;
    for (auto const& [track1, track2] : soa::combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      // Discard same track
      if (track1.index() == track2.index()) {
        continue;
      }

      // Discard same charge track
      if (track1.sign() == track2.sign()) {
        continue;
      }

      // Select track
      if (!selectTrack(track1) || !selectTrack(track2)) {
        continue;
      }

      // Identify track
      if (!identifyTrack<kKa>(track1) || !identifyTrack<kKa>(track2)) {
        continue;
      }

      // Apply rapidity acceptance
      std::array<float, 3> v = {track1.px() + track2.px(), track1.py() + track2.py(), track1.pz() + track2.pz()};
      if (RecoDecay::y(v, MassPhi) >= cResRapCut) {
        continue;
      }

      // Reconstruct invariant mass
      float p = RecoDecay::p((track1.px() + track2.px()), (track1.py() + track2.py()), (track1.pz() + track2.pz()));
      float e = RecoDecay::e(track1.px(), track1.py(), track1.pz(), MassKaonCharged) + RecoDecay::e(track2.px(), track2.py(), track2.pz(), MassKaonCharged);
      float m = std::sqrt(RecoDecay::m2(p, e));

      // Get directed flow
      ux = std::cos(RecoDecay::phi(v));
      uy = std::sin(RecoDecay::phi(v));
      v1a = ux * vSP[kXa] + uy * vSP[kYa];
      v1c = ux * vSP[kXc] + uy * vSP[kYc];

      // Fill signal histogram
      histos.fill(HIST("Reso/Phi/hSigCentPtInvMass"), cent, RecoDecay::pt(v), m);
      histos.fill(HIST("Reso/Phi/Sig/hPhiQuA"), cent, RecoDecay::y(v, MassPhi), m, v1a);
      histos.fill(HIST("Reso/Phi/Sig/hPhiQuC"), cent, RecoDecay::y(v, MassPhi), m, v1c);

      // Get background
      p = RecoDecay::p((track1.px() - track2.px()), (track1.py() - track2.py()), (track1.pz() - track2.pz()));
      m = std::sqrt(RecoDecay::m2(p, e));
      v[0] = track1.px() - track2.px();
      v[1] = track1.py() - track2.py();
      v[2] = track1.pz() - track2.pz();
      ux = std::cos(RecoDecay::phi(v));
      uy = std::sin(RecoDecay::phi(v));
      v1a = ux * vSP[kXa] + uy * vSP[kYa];
      v1c = ux * vSP[kXc] + uy * vSP[kYc];

      // Fill bkg histogram
      histos.fill(HIST("Reso/Phi/hBkgCentPtInvMass"), cent, RecoDecay::pt(v), m);
      histos.fill(HIST("Reso/Phi/Bkg/hPhiQuA"), cent, RecoDecay::y(v, MassPhi), m, v1a);
      histos.fill(HIST("Reso/Phi/Bkg/hPhiQuC"), cent, RecoDecay::y(v, MassPhi), m, v1c);
    }
  }

  template <ParticleType part, typename T>
  void getIdHadronFlow(float const& cent, T const& track, float const& v1a, float const& v1c)
  {
    static constexpr std::string_view SubDir[] = {"Pion/", "Kaon/", "Proton/"};
    float tpcNsigma = 0., tofNsigma = 0.;
    if (part == kPi) {
      tpcNsigma = track.tpcNSigmaPi();
      tofNsigma = track.tofNSigmaPi();
    } else if (part == kKa) {
      tpcNsigma = track.tpcNSigmaKa();
      tofNsigma = track.tofNSigmaKa();
    } else if (part == kPr) {
      tpcNsigma = track.tpcNSigmaPr();
      tofNsigma = track.tofNSigmaPr();
    } else {
      return;
    }
    histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hdEdX"), track.pt(), track.tpcSignal());
    histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hTPCNSigma"), track.pt(), tpcNsigma);
    if (track.hasTOF()) {
      histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hTOFNSigma"), track.pt(), tofNsigma);
    }
    if (track.sign() > 0) {
      histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hAQuPos"), cent, track.eta(), v1a);
      histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hCQuPos"), cent, track.eta(), v1c);
    } else {
      histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hAQuNeg"), cent, track.eta(), v1a);
      histos.fill(HIST("PartId/") + HIST(SubDir[part]) + HIST("hCQuNeg"), cent, track.eta(), v1c);
    }
  }

  void fillCorrHist(std::array<float, 4> const& vCollParam, std::array<float, 4> const& vSP)
  {
    histos.fill(HIST("CorrHist/hWtXZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kXa]);
    histos.fill(HIST("CorrHist/hWtYZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kYa]);
    histos.fill(HIST("CorrHist/hWtXZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kXc]);
    histos.fill(HIST("CorrHist/hWtYZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kYc]);
    histos.fill(HIST("CorrHist/hUWtXZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtYZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtXZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtYZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hXZNAVsCent"), vCollParam[kCent], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVx"), vCollParam[kVx], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVy"), vCollParam[kVy], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVz"), vCollParam[kVz], vSP[kXa]);
    histos.fill(HIST("CorrHist/hYZNAVsCent"), vCollParam[kCent], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVx"), vCollParam[kVx], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVy"), vCollParam[kVy], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVz"), vCollParam[kVz], vSP[kYa]);
    histos.fill(HIST("CorrHist/hXZNCVsCent"), vCollParam[kCent], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVx"), vCollParam[kVx], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVy"), vCollParam[kVy], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVz"), vCollParam[kVz], vSP[kXc]);
    histos.fill(HIST("CorrHist/hYZNCVsCent"), vCollParam[kCent], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVx"), vCollParam[kVx], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVy"), vCollParam[kVy], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVz"), vCollParam[kVz], vSP[kYc]);
  }

  template <typename T>
  void fillTrackHist(T const& track)
  {
    histos.fill(HIST("TrackQA/hPtDcaZ"), track.pt(), track.dcaZ());
    histos.fill(HIST("TrackQA/hPtDcaXY"), track.pt(), track.dcaXY());
  }

  template <typename C>
  bool analyzeCollision(C const& collision, std::array<float, 4>& vSP)
  {
    // Event selection
    if (!selCollision(collision)) {
      return false;
    }
    posX = collision.posX();
    posY = collision.posY();
    posZ = collision.posZ();
    std::array<float, 4> vCollParam = {cent, posX, posY, posZ};

    // Fill event QA
    histos.fill(HIST("Event/hCent"), cent);
    histos.fill(HIST("Event/hVx"), posX);
    histos.fill(HIST("Event/hVy"), posY);
    histos.fill(HIST("Event/hVz"), posZ);

    // Get bunch crossing
    auto bc = collision.template foundBC_as<BCsRun3>();
    cRunNum = collision.template foundBC_as<BCsRun3>().runNumber();

    // Load calibration flags
    bool loadGainCalib = false, loadShiftCorr = false;
    if (cRunNum != lRunNum) {
      loadGainCalib = true;
      loadShiftCorr = true;
    }

    // check zdc
    if (!bc.has_zdc()) {
      lRunNum = cRunNum;
      return false;
    }

    auto zdc = bc.zdc();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergyCommon = zdc.energyCommonZNA();
    auto zncEnergyCommon = zdc.energyCommonZNC();

    // check energy deposits
    if (znaEnergyCommon <= 0 || zncEnergyCommon <= 0 || znaEnergy[0] <= 0 || znaEnergy[1] <= 0 || znaEnergy[2] <= 0 || znaEnergy[3] <= 0 || zncEnergy[0] <= 0 || zncEnergy[1] <= 0 || zncEnergy[2] <= 0 || zncEnergy[3] <= 0) {
      lRunNum = cRunNum;
      return false;
    }

    // Fill gain calib histograms
    for (int iCh = 0; iCh < kXYAC; ++iCh) {
      histos.fill(HIST("QA/hZNASignal"), iCh + 0.5, vCollParam[kVz], znaEnergy[iCh]);
      histos.fill(HIST("QA/hZNCSignal"), iCh + 0.5, vCollParam[kVz], zncEnergy[iCh]);
      histos.fill(HIST("QA/hZNAEnergyCommon"), vCollParam[kVz], znaEnergyCommon);
      histos.fill(HIST("QA/hZNCEnergyCommon"), vCollParam[kVz], zncEnergyCommon);
    }

    // Do gain calibration
    if (cDoGainCalib) {
      gainCalib(loadGainCalib, vCollParam[kVz], znaEnergy, zncEnergy);
    }

    // Fill zdc signal
    for (int iCh = 0; iCh < kXYAC; ++iCh) {
      histos.fill(HIST("QA/GainCalib/hZNASignal"), iCh + 0.5, znaEnergy[iCh]);
      histos.fill(HIST("QA/GainCalib/hZNCSignal"), iCh + 0.5, zncEnergy[iCh]);
    }

    auto alphaZDC = 0.395;
    const double x[4] = {-1.75, 1.75, -1.75, 1.75};
    const double y[4] = {-1.75, -1.75, 1.75, 1.75};

    // Calculate X and Y
    float znaXNum = 0., znaYNum = 0., zncXNum = 0., zncYNum = 0.;
    float znaDen = 0., zncDen = 0.;
    float znaWt = 0., zncWt = 0.;

    // Loop over zdc sectors
    for (int i = 0; i < kXYAC; ++i) {
      if (cUseAlphaZDC) {
        znaWt = std::pow(znaEnergy[i], alphaZDC);
        zncWt = std::pow(zncEnergy[i], alphaZDC);
      } else {
        znaWt = znaEnergy[i];
        zncWt = zncEnergy[i];
      }
      znaXNum -= znaWt * x[i];
      znaYNum += znaWt * y[i];
      zncXNum += zncWt * x[i];
      zncYNum += zncWt * y[i];
      znaDen += znaWt;
      zncDen += zncWt;
    }

    if (znaDen < zdcDenThrs || zncDen < zdcDenThrs) {
      lRunNum = cRunNum;
      return false;
    }

    // Get X and Y for A and C side ZNA
    vSP[kXa] = znaXNum / znaDen;
    vSP[kYa] = znaYNum / znaDen;
    vSP[kXc] = zncXNum / zncDen;
    vSP[kYc] = zncYNum / zncDen;

    // Do corrections
    if (cApplyRecentCorr) {
      applyCorrection(loadShiftCorr, vCollParam, vSP);
    }

    // Fill X and Y histograms for corrections after each iteration
    fillCorrHist(vCollParam, vSP);
    return true;
  }

  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::MultsExtra>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::TOFSignal, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPr, aod::pidTOFPr, aod::TrackCompColls>;

  void process(CollisionsRun3::iterator const& collision, BCsRun3 const& /* bcs*/, aod::Zdcs const&, Tracks const& tracks)
  {
    // Analyze collison and get SP vector
    std::array<float, 4> vSP = {0., 0., 0., 0.};
    if (!analyzeCollision(collision, vSP)) {
      lRunNum = cRunNum;
      return;
    }

    // Evaluate spectator plane angle and [X,Y] correlations
    float psiA = std::atan2(vSP[kYa], vSP[kXa]);
    float psiC = std::atan2(vSP[kYc], vSP[kXc]);
    histos.fill(HIST("Checks/hPsiSPA"), cent, psiA);
    histos.fill(HIST("Checks/hPsiSPC"), cent, psiC);
    histos.fill(HIST("Checks/hCosPsiSPAC"), cent, std::cos(psiA - psiC));
    histos.fill(HIST("Checks/hSinPsiSPAC"), cent, std::sin(psiA - psiC));
    histos.fill(HIST("Checks/hXaXc"), cent, (vSP[kXa] * vSP[kXc]));
    histos.fill(HIST("Checks/hYaYc"), cent, (vSP[kYa] * vSP[kYc]));
    histos.fill(HIST("Checks/hXaYc"), cent, (vSP[kXa] * vSP[kYc]));
    histos.fill(HIST("Checks/hYaXc"), cent, (vSP[kYa] * vSP[kXc]));

    // Directed flow
    float qac = (vSP[kXa] * vSP[kXc]) + (vSP[kYa] * vSP[kYc]);
    histos.fill(HIST("DF/hQaQc"), cent, qac);

    // Loop over tracks
    float ux = 0., uy = 0., v1a = 0., v1c = 0.;
    for (auto const& track : tracks) {
      // Select track
      if (!selectTrack(track)) {
        continue;
      }

      // Fill track QA
      fillTrackHist(track);

      // Get directed flow
      ux = std::cos(track.phi());
      uy = std::sin(track.phi());
      v1a = ux * vSP[kXa] + uy * vSP[kYa];
      v1c = ux * vSP[kXc] + uy * vSP[kYc];

      // Charged particle directed flow
      histos.fill(HIST("DF/hAQu"), cent, track.eta(), v1a);
      histos.fill(HIST("DF/hCQu"), cent, track.eta(), v1c);
      if (track.sign() > 0) {
        histos.fill(HIST("DF/hAQuPos"), cent, track.eta(), v1a);
        histos.fill(HIST("DF/hCQuPos"), cent, track.eta(), v1c);
      } else {
        histos.fill(HIST("DF/hAQuNeg"), cent, track.eta(), v1a);
        histos.fill(HIST("DF/hCQuNeg"), cent, track.eta(), v1c);
      }

      // Identified directed flow
      if (identifyTrack<kPi>(track)) {
        getIdHadronFlow<kPi>(cent, track, v1a, v1c);
      } else if (identifyTrack<kKa>(track)) {
        getIdHadronFlow<kKa>(cent, track, v1a, v1c);
      } else if (identifyTrack<kPr>(track)) {
        getIdHadronFlow<kPr>(cent, track, v1a, v1c);
      }
    }

    // Get resonance flow
    getResoFlow(tracks, vSP);

    // Update runnumber
    lRunNum = cRunNum;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEventPlane>(cfgc)};
}
